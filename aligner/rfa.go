package aligner

import (
	"bufio"
	"fmt"
	"os"
	"sync"

	"arachne/fastqreader"
	"arachne/gobwa"
	"arachne/optimizer"

	"github.com/biogo/hts/sam"
)

// Holds configuration parameters for the optimizer. These to be constant after set in main.
type RFAConfig struct {
	improper_penalty float64
}

// Holds statistics and coordination points for RFA
type RFAStats struct {
	total                 int64
	correct               int64
	correct_mapq10        int64
	total_mapq10          int64
	total_improper        int64
	total_improper_before int64
	mapq                  int64

	// This lock synchronizes access to the "mapq.csv" file
	file *bufio.Writer
	lock sync.Mutex
}

type Optimizer struct {
	candidate_molecules       []*CandidateMolecule
	alignments                [][]*Alignment
	total_alignment_score     float64
	currentMoleculeMoveSource int
	currentScore              float64
	log_unpaired_probability  float64
	barcode                   string
}

func DoRFAForOneBarcode(work *WorkUnit,
	out chan *sam.Record,
	ref *gobwa.GoBwaReference,
	settings *gobwa.GoBwaSettings,
	config *RFAConfig,
	stats *RFAStats,
	contigs map[string]*sam.Reference,
	debugtags *bool,
	reads []fastqreader.FastQRecord) {
	var worthRFA bool
	//TODO STATS ARENT USED ANYWHERE?
	stats.total = 0
	stats.mapq = 0
	//barcode_num := work.barcodenum
	barcode_reads := work.reads
	arena := gobwa.NewArena()
	if !work.reads[0].Valid {
		worthRFA = worthRunningRFA(barcode_reads, work.unique_barcode)
	}
	barcode_chains, barcode := GetChains(ref, settings, barcode_reads, arena, 25)
	alignments, stashed_alignments := GetAlignments(ref, settings, barcode_chains, 17, arena)
	//stashed_alignments = StashAlignments(alignments);

	//	positions := tagBestAlignments(alignments, -17)
	positions := tagBestAlignments(alignments)

	if len(barcode_reads) > 2 {
		if *verbose {
			fmt.Fprintf(os.Stderr, "working on barcode %s  num reads: %d  doing RFA: %v  unique_barcode %v\n",
				string(barcode_reads[0].Barcode),
				len(barcode_reads),
				worthRFA,
				work.unique_barcode)
		}
	}

	if !worthRFA {
		//estimateMapQualities(-1, alignments, nil, config.improper_penalty, stats)
		estimateMapQualities(alignments, nil, config.improper_penalty)
		markDuplicates(alignments)
		CheckSplitReads(stashed_alignments, centromeres)
		flushToChannel(alignments, out, contigs, debugtags)
		ReturnBuffer(reads) // was inside BamThread after DoDumpToBam, move here
		arena.Free()
		return
	}

	candidate_molecules := inferMolecules(positions)
	markBestAlignmentForReadInMolecule(candidate_molecules)
	candidate_molecules = scrapMolecules(candidate_molecules)

	setMoleculeDifferences(candidate_molecules, false)

	optimizer_obj := &Optimizer{
		candidate_molecules:       candidate_molecules,
		alignments:                alignments,
		currentMoleculeMoveSource: 0,
		log_unpaired_probability:  config.improper_penalty,
		barcode:                   barcode,
	}

	optimized := optimizer.Optimize(optimizer.Optimizable(*optimizer_obj), 1, 2, 4*len(candidate_molecules)).(Optimizer)

	//estimateMapQualities(barcode_num, optimized.alignments, optimized.candidate_molecules, optimized.log_unpaired_probability, stats)
	estimateMapQualities(optimized.alignments, optimized.candidate_molecules, optimized.log_unpaired_probability)
	markDuplicates(alignments)
	CheckSplitReads(stashed_alignments, centromeres)
	flushToChannel(alignments, out, contigs, debugtags)
	ReturnBuffer(reads) // was inside BamThread after DoDumpToBam, move here
	arena.Free()
}

// Was the deconv format super necessary?

// Determine if there are enough reads (5) to run RFA
func worthRunningRFA(barcode_reads []fastqreader.FastQRecord, uniqueBarcode bool) bool {
	if len(barcode_reads) == 0 || !uniqueBarcode {
		return false
	}
	//bcParts := strings.Split(string(barcode_reads[0].Barcode), "-")
	//if len(bcParts) < 2 {
	//	return false
	//}
	if len(barcode_reads) < 5 {
		return false
	}
	return true
}

func (o Optimizer) fastScore(sourceMolecule, sinkMolecule *CandidateMolecule) (float64, Move) {
	change, move := fastScore(sourceMolecule, sinkMolecule, o.log_unpaired_probability)
	return change, move
}

func acceptMove(move Move) {
	toDelete := move.toDelete
	toSet := move.toSet
	if *debugPrintMove {
		fmt.Println("Accepting move from ", move.source.start, " to ", move.sink.start)
	}
	for i := range toDelete {
		read_id := toDelete[i]
		sinkAlignment := toSet[i]
		sourceAlignment := move.source.active_alignments.Get(read_id)
		for _, mismatchLoc := range sourceAlignment.mismatchLocs {
			num, has := move.source.mismatchLocs[mismatchLoc]
			if !has || num == 0 {
				//there is a problem
				panic("source molecule should have this entry")
			}
			if *debugPrintMove {
				fmt.Println("removing mismatchLoc", mismatchLoc, move.source.mismatchLocs[mismatchLoc])
			}
			move.source.mismatchLocs[mismatchLoc]--
			if *debugPrintMove {
				fmt.Println(move.source.mismatchLocs[mismatchLoc])
			}
		}
		for _, mismatchLoc := range sinkAlignment.mismatchLocs {
			_, has := move.sink.mismatchLocs[mismatchLoc]
			if has {
				move.sink.mismatchLocs[mismatchLoc]++
			} else {
				move.sink.mismatchLocs[mismatchLoc] = 1
			}
		}
		move.source.active_alignments.Delete(read_id)
		move.sink.active_alignments.Set(read_id, sinkAlignment)
		sourceAlignment.active = false
		sinkAlignment.active = true
	}
}
