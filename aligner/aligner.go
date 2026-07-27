// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package aligner

import (
	"C"
	"bufio"
	"crypto/md5"
	"encoding/binary"
	"fmt"
	"math"
	"math/rand"
	"os"
	"runtime"
	"sort"
	"strconv"
	"sync"
	"unsafe"

	"arachne/fastqreader"
	"arachne/gobwa"
	"arachne/optimizer"

	"github.com/biogo/hts/sam"
)

var VERSION string

type ArachneArgs struct {
	R1                    *string
	R2                    *string
	Improper_pair_penalty *float64
	//Read_groups           *string
	Sample_id         *string
	Threads           *int
	DEBUG             *bool
	PositionChunkSize *int
	InferDistance     *int64
	DebugTags         *bool
	DebugPrintMove    *bool
	Reference         *string
	Centromeres       *string
	Verbose           *bool
	Comments          *bool
}

type ChainedHit struct {
	contig    string
	pos       int64
	aend      int64
	read_id   int
	mate_id   int
	hit_id    int
	secondary bool
	read1     bool
	score     int
	chain     unsafe.Pointer //*C.mem_alnreg_t
	read      *[]byte
	fastq     *fastqreader.FastQRecord
	aln       *gobwa.EasyAlignment
	// trim_seq  *[]byte
	// trim_qual *[]byte
}

func FindRead(alignments [][]*Alignment, molecules []*CandidateMolecule, qname string) {
	for _, alignmentArray := range alignments {
		for _, alignment := range alignmentArray {

			if *alignment.read_name == qname && alignment.active {
				fmt.Println("printing", alignment.read1)
				alignment.Print()
				fmt.Println("and its mate")
				if alignment.mate_alignment != nil {
					alignment.mate_alignment.Print()
				}
				if alignment.molecule_id != -1 {
					fmt.Println(molecules[alignment.molecule_id].active_alignments.Len())
				}
			}
		}
	}
}

// types and functions to be able to sort a list of aligntments by position,
// must already be of the same contig otherwise it doesnt make sense
type ByPosition []*Alignment

func (a ByPosition) Len() int           { return len(a) }
func (a ByPosition) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a ByPosition) Less(i, j int) bool { return a[i].pos < a[j].pos }

type CandidateMolecule struct {
	id                      int
	chrom                   string
	start                   int64
	stop                    int64
	alignments              *OrderedMap          //map[int]map[int]*Alignment read id to map of alignment id to alignment
	best_alignment_for_read *OrderedAlignmentMap //map[int]*Alignment read id to alignment
	active_alignments       *OrderedAlignmentMap //map[int]*Alignment read id to alignment
	log_probability         float64
	true_molecule           bool //only for simulation
	active_molecule         bool
	molecule_confidence     float64
	differences             float64
	soft_clipped            int
	mismatchLocs            map[int]int
}

// Holds a single unit of "work" to be processed by an RFA thread
type WorkUnit struct {
	reads          []fastqreader.FastQRecord
	barcode_num    int
	unique_barcode bool
	poolBuf        *[]fastqreader.FastQRecord
}

var barcode_reads [][]fastqreader.FastQRecord
var barcode_reads_lock sync.Mutex

func ReturnBuffer(reads []fastqreader.FastQRecord) {
	barcode_reads_lock.Lock()
	barcode_reads = append(barcode_reads, reads)
	barcode_reads_lock.Unlock()
}

// Command line arguments
var r1 *string
var r2 *string
var improper_pair_penalty *float64

// var read_groups *string
var sample_id *string
var threads *int
var debugTags *bool
var reference *string
var DEBUG *bool
var debugPrintMove *bool
var centromeres map[string]Region
var verbose *bool
var inferDistance *int64
var AddComments *bool

// this is the actual arachne program
func Arachne(args ArachneArgs) {
	r1 = args.R1
	r2 = args.R2
	improper_pair_penalty = args.Improper_pair_penalty
	//read_groups = args.Read_groups
	sample_id = args.Sample_id
	threads = args.Threads
	debugTags = args.DebugTags
	reference = args.Reference
	centromeres = loadCentromeres(args.Centromeres)
	inferDistance = args.InferDistance
	AddComments = args.Comments
	verbose = args.Verbose
	// unused
	DEBUG = args.DEBUG
	debugPrintMove = args.DebugPrintMove
	runtime.GOMAXPROCS(*threads + 2)
	fastq, err := fastqreader.OpenFastQPair(*r1, *r2)
	if err != nil {
		panic(err)
	}
	defer fastq.Close()

	// ── process reference ─────────────────────────────────────────────────────────
	if *verbose {
		fmt.Fprintf(os.Stderr, "Loading reference: %s\n", *reference)
	}
	ref := gobwa.GoBwaLoadReference(*reference)
	if *verbose {
		fmt.Fprint(os.Stderr, "Reference loaded\n")
	}
	settings := gobwa.GoBwaAllocSettings()
	config := &RFAConfig{*improper_pair_penalty}
	config.improper_penalty = float64(*improper_pair_penalty)

	// ------- SAM output writer -------------------------
	// chanCap sized to absorb worker bursts
	header, contigs := buildHeader(ref, *sample_id, VERSION)
	writeChannel, doneChan := NewSamWriterChannel(header, *threads*500, 2<<20, *threads)

	// ── buffer pool ─────────────────────────────────────────────────────────
	// One extra buffer beyond *threads so the main loop always has something
	// to fill while all workers are busy with their own buffers.
	bufChan := make(chan *[]fastqreader.FastQRecord, *threads+1)
	for range *threads + 1 {
		s := make([]fastqreader.FastQRecord, 0, 50000)
		bufChan <- &s
	}

	var w *bufio.Writer
	stats := &RFAStats{}
	stats.file = w

	// ── workers ─────────────────────────────────────────────────────────────
	work_to_do := make(chan *WorkUnit, 2)
	var wg sync.WaitGroup
	for range *threads {
		wg.Go(
			func() {
				WorkerThread(work_to_do, writeChannel, bufChan, ref, settings, config, stats, contigs, debugTags)
			})
	}
	//finished := make (chan bool);
	barcode_reads = [][]fastqreader.FastQRecord{}

	// ── main loop ───────────────────────────────────────────────────────────
	for barcode_num := 1; ; barcode_num++ {
		buf := <-bufChan
		*buf = (*buf)[:0]

		records, err, full_barcode := fastq.ReadBarcodeSet(buf)
		if err != nil {
			bufChan <- buf // return unused buffer
			break
		}
		if !records[0].Valid {
			if *verbose {
				fmt.Printf("Detected invalid barcode (%v), skipping.\n", records[0].Barcode)
			}
			bufChan <- buf // return unused buffer on skip too
			continue
		}
		work_to_do <- &WorkUnit{records, barcode_num, full_barcode, buf}
	}

	// ── teardown ────────────────────────────────────────────────────────────
	close(work_to_do)
	wg.Wait()
	close(writeChannel)
	<-doneChan
}

// A single "worker" thread. It tries to grab work units until it gets nil, then it shuts down.
func WorkerThread(
	input chan *WorkUnit,
	out chan *sam.Record,
	bufChan chan *[]fastqreader.FastQRecord,
	ref *gobwa.GoBwaReference,
	settings *gobwa.GoBwaSettings,
	config *RFAConfig,
	stats *RFAStats,
	contigs map[string]*sam.Reference,
	debugtags *bool) {
	for work := range input {
		DoRFAForOneBarcode(work, out, ref, settings, config, stats, contigs, debugtags, work.reads)
		bufChan <- work.poolBuf // return buffer for reuse
	}
}

func StashAlignments(a *OrderedMap) *OrderedMap {
	res := NewOrderedMap()
	for id, read := range a.Iter() {
		a_for_read := NewOrderedMap()
		res.Set(id, a_for_read)
		alignments := FixGetForTypeOrderedMap(read)
		for id2, a_ptr := range alignments.Iter() {
			aln := FixGetForTypeAlignment(a_ptr)
			a_for_read.Set(id2, aln)
		}
	}
	return res
}

func DeAlignCrappyReads(reads [][]*Alignment) {
	for _, readArray := range reads {
		for _, aln := range readArray {
			if !aln.is_proper && aln.score-17 < 19 {
				aln.pos = -1
			}
		}
	}
}

/*
TODO FLAG FOR REMOVAL IT'S UNUSED

	func unbarcodeAlignments(alignments [][]*Alignment) {
		for _, alignmentList := range alignments {
			for _, alignment := range alignmentList {
				if alignment.active {
					barcode := []byte(strings.Split(string(*alignment.barcode), "-")[0])
					alignment.barcode = &barcode
				}
			}
		}
	}
*/
func setMoleculeDifferences(candidate_molecules []*CandidateMolecule, setBad bool) {
	for i := range candidate_molecules {
		differences := 0
		for _, alignment := range candidate_molecules[i].active_alignments.Iter() {
			differences += alignment.mismatches
		}
		if setBad {
			if float64(differences)/float64(candidate_molecules[i].active_alignments.Len()) > 1.5 {
				for _, alignment := range candidate_molecules[i].active_alignments.Iter() {
					alignment.bad_molecule = true
				}
			}
		}
		candidate_molecules[i].differences = float64(differences) / float64(candidate_molecules[i].active_alignments.Len())
		for _, alignment := range candidate_molecules[i].active_alignments.Iter() {
			alignment.molecule_difference = candidate_molecules[i].differences
		}
	}
}

func psuedoCountAlignmentScore(aln *Alignment, log_molecule_penalty float64) float64 {
	psuedoAlignmentLength := 25.0
	//score := 0.0
	//score -= 10.0                                                        //maximum soft clipping penalty
	score := -10.0
	score -= (float64(len(*aln.read_seq)) - psuedoAlignmentLength) * 0.5 // soft clipping length penalty for 25bp alignment
	score += log_molecule_penalty
	return score
}

func scoreAlignment(aln *Alignment, mate *Alignment, log_molecule_penalty float64) float64 {
	score := 0.0
	if aln != nil {
		score += float64(aln.mismatches*-2.0 + aln.indels*-3.0)
		if aln.soft_clipped > 0 {
			score -= 5.0 * float64(aln.soft_clipped)
			score -= float64(aln.soft_clipped_length) * 0.5
		}
	}
	if mate != nil {
		score += float64(mate.mismatches*-2.0 + mate.indels*-3.0)
		if mate.soft_clipped > 0 {
			score -= 5.0 * float64(mate.soft_clipped)
			score -= float64(mate.soft_clipped_length) * 0.5
		}
	}
	if mate == nil || aln == nil || !isPair(aln, mate) {
		score += *improper_pair_penalty
	}
	if aln != nil {
		if !aln.active_molecule {
			score += log_molecule_penalty
		}
	}
	return score
}

func updateAlignmentsMoleculeStatus(alignments [][]*Alignment, candidate_molecules []*CandidateMolecule, read_copies_in_active_molecule, read_copies_not_in_active_molecule map[int]int, unique_molecules_active map[int]map[int]bool) {
	if candidate_molecules != nil {
		setMoleculeConfidences(candidate_molecules)
		setMoleculeDifferences(candidate_molecules, false)

		// #1 update alignment probabilities based on whether they are in an active molecule or not
		// here we take 2 sources of power. 1 is taht singletons are rare
		for read_id, alignmentArray := range alignments {
			for _, alignment := range alignmentArray {
				is_molecule_active := false
				if alignment.molecule_id != -1 {
					molecule := candidate_molecules[alignment.molecule_id]
					is_molecule_active = molecule.active_alignments.Len()-molecule.soft_clipped > 4 && molecule.molecule_confidence > 0.1
					alignment.active_molecule = is_molecule_active
				}
				if is_molecule_active {
					candidate_molecules[alignment.molecule_id].active_molecule = true
					read_copies_in_active_molecule[read_id]++ //TODO remove, book keeping
					_, has := unique_molecules_active[read_id]
					if !has {
						unique_molecules_active[read_id] = map[int]bool{}
					}
					unique_molecules_active[read_id][alignment.molecule_id] = true
				} else {
					read_copies_not_in_active_molecule[read_id]++ //TODO remove, book keeping
				}
				if alignment.molecule_id != -1 {
					alignment.mapq_data.reads_in_molecule = candidate_molecules[alignment.molecule_id].active_alignments.Len()
				}
			}
		}
	}
}

func appendPsuedocountAlignmentScore(scores []float64, alignmentArray []*Alignment, alignments [][]*Alignment, log_molecule_penalty float64) []float64 {
	// get psuedocount alignment by pairing psuedocount alignment with best single mate
	if len(alignmentArray) > 0 {
		mateArray := alignments[alignmentArray[0].mate_id]
		bestSingleMateScore := -math.MaxFloat64
		for _, mateAlignment := range mateArray {
			scoreMateSingle := scoreAlignment(nil, mateAlignment, log_molecule_penalty)
			if scoreMateSingle > bestSingleMateScore {
				bestSingleMateScore = scoreMateSingle
			}
		}
		if len(mateArray) > 0 {
			scores = append(scores, bestSingleMateScore+psuedoCountAlignmentScore(alignmentArray[0], log_molecule_penalty))
		} else {
			scores = append(scores, psuedoCountAlignmentScore(alignmentArray[0], log_molecule_penalty))
		}
	}
	return scores
}

func moleculeMapqProbabilitySums(candidate_molecules []*CandidateMolecule, log_unpaired_probability float64) {
	for mol_id, sourceMolecule := range candidate_molecules {
		for mol2_id, sinkMolecule := range candidate_molecules {
			if mol_id == mol2_id {
				continue
			}
			sourceAlignments := []*Alignment{}
			for _, aln := range sourceMolecule.active_alignments.Iter() {
				alt_aln := sinkMolecule.best_alignment_for_read.Get(aln.read_id)
				if alt_aln != nil {
					sourceAlignments = append(sourceAlignments, aln)
				}
			}
			sourceSinkChange, _ := fastScore(sourceMolecule, sinkMolecule, log_unpaired_probability)
			moleculeMoveProbability := math.Pow(10, sourceSinkChange)
			for _, alignment := range sourceAlignments {
				if !alignment.active {
					panic("setting molecule mapq for non active alignment")
				}
				alignment.sum_move_probability_change += moleculeMoveProbability
			}
		}
	}
}

func calculateLogMoleculePenalty(candidate_molecules []*CandidateMolecule, referenceLength float64) float64 {
	dnaLength := 1000.0
	if len(candidate_molecules) == 0 {
		return 0.0
	}
	for i := range candidate_molecules {
		mol := candidate_molecules[i]
		if mol.active_molecule {
			smallest := int64(math.MaxInt64)
			biggest := int64(-1)
			for _, alignment := range mol.active_alignments.Iter() {
				if alignment.pos > biggest {
					biggest = alignment.pos
				}
				if alignment.pos < smallest {
					smallest = alignment.pos
				}
			}
			if biggest >= smallest {
				dnaLength += float64(biggest-smallest) + 1000.0
			}
		} else {
			for _, alignment := range mol.active_alignments.Iter() {
				dnaLength += float64(alignment.aend-alignment.pos) * 2.0
			}
		}
	}
	// singletonProb := 0.05
	moleculePenalty := math.Log10(dnaLength / referenceLength * 0.05)
	return moleculePenalty
}

func checkMates(alignments [][]*Alignment) {
	for _, alignmentArray := range alignments {
		for _, alignment := range alignmentArray {
			if alignment.active {
				if alignment.mate_alignment == nil {
					continue
				}
				if !alignment.mate_alignment.active {
					var activeMate *Alignment
					for _, mateAln := range alignments[alignment.mate_id] {
						if mateAln.active {
							activeMate = mateAln
						}
					}
					panic(fmt.Sprint("something wrong with mate of ", alignment.id, alignment.read_id, " is ", alignment.mate_alignment.id, alignment.mate_alignment.read_id, "actually active id", activeMate.id, activeMate.read_id))
				}
			}
		}
	}
}

/*
So the basic strategy here is two-fold
 1. is based on calculating a mapq vs being able to move that read to any of the other
    alignments that it has. This includes the probabilities of a read being a singleton vs being in an active molecule.
    Every alignment has some probability calculated from its alignment score/cigar string. Then alignments are
    penalized for improper_pair and for being outside of active molecules. The mapq is then calculated by normalizing the sum of
    alignment probabilities to 1 and then the mapq is -10*log10(1-probability_of_chosen_alignment)
    If there is a problem with this method it is that it dings alignments for not being in an active read, but if a whole
    molecule could have moved to that location with very little probability change making it an active molecule, the read
    shouldn't have gotten that ding. I have tried not dinging reads in candidate molecules that are inactive but have
    greater than some number of alignments but have found it better to just include method #2
 2. is a mapq considering the probability change in the case of a candidate molecule making a sub-move (moving all of the reads
    that can move from this molecule to another candidate molecule) to every candidate molecule that it can make a sub move to.
    So for each active molecule we calculate the probability change of a sub-move to every candidate molecule that it can. For each
    sub-move calculated, every alignment that was involved in that sub move gets that probability change added to a sum.
    The log probability change of making a sub-move from a molecule to itself is 0 so probability 1.0. then we normalize the probability
    of the final selected sub-move which is the molecule moving to itself p = 1.0 vs the sum of all of the probabilities of the moves
    which is p_final = 1.0/sum_i(p(move_i)) and then calculate a mapq from that probability = -10*log10(1-p_final)

Finally we take the min of the two approaches as the final mapq
*/
func estimateMapQualities( //barcode int, //TODO remove, this isn't used
	alignments [][]*Alignment,
	candidate_molecules []*CandidateMolecule,
	log_unpaired_probability float64,
	//stats *RFAStats    //TODO remove, this isn't used
) {
	read_copies_in_active_molecule := map[int]int{}     //TODO remove, book keeping
	read_copies_not_in_active_molecule := map[int]int{} //TODO remove, book keeping
	unique_molecules_active := map[int]map[int]bool{}
	// mapq strategy 2: sum probabilities of full molecule moves
	if *debugPrintMove {
		fmt.Println("NOW TESTING MAPQS")
	}
	moleculeMapqProbabilitySums(candidate_molecules, log_unpaired_probability)

	// Now to update alignment probabilities for being singleton/outside active molecules
	// this part only happens if we ran RFA, bad barcodes etc get no more probability updates
	updateAlignmentsMoleculeStatus(alignments, candidate_molecules, read_copies_in_active_molecule, read_copies_not_in_active_molecule, unique_molecules_active)
	log_molecule_penalty := calculateLogMoleculePenalty(candidate_molecules, 3200000000.0) //hard coding length of human reference
	//now go through every read_id and normalize all alternate alignment probabilities
	for read_id, alignmentArray := range alignments {
		// find best pair for alignments and make list of those alignment pair scores for use of probability normalization to sum to 1.0
		scores := []float64{}
		scores = appendPsuedocountAlignmentScore(scores, alignmentArray, alignments, log_molecule_penalty)
		total_probability := float64(0.0)
		for _, alignment := range alignmentArray {
			mateArray := alignments[alignment.mate_id]
			for _, mateAlignment := range mateArray {
				if alignment.active && mateAlignment.active {
					alignment.mate_alignment = mateAlignment
					mateAlignment.mate_alignment = alignment
				}
			}
		}

		for _, alignment := range alignmentArray {
			mateArray := alignments[alignment.mate_id]
			best_score := -math.MaxFloat64
			for _, mateAlignment := range mateArray {
				score := scoreAlignment(alignment, mateAlignment, log_molecule_penalty)
				if score > best_score {
					best_score = score
				}
			}
			if len(mateArray) == 0 {
				best_score = scoreAlignment(alignment, nil, log_molecule_penalty)
			}
			scores = append(scores, best_score)
		}

		// gather and record info about the second best pair alignment
		second_best_proper_pair := false
		second_best_raw_score := scores[0] //psuedoCountAlignmentScore
		second_best_log_probability := -1000.0
		second_best_molecule_reads := -1
		var second_best_alignment *Alignment
		second_best_molecule_confidence := -1.0
		for _, alignment := range alignmentArray {
			mateArray := alignments[alignment.mate_id]
			for _, mateAlignment := range mateArray {
				score := scoreAlignment(alignment, mateAlignment, log_molecule_penalty)
				if !alignment.active {
					if score > second_best_log_probability {
						second_best_log_probability = score
						second_best_raw_score = scoreAlignment(alignment, mateAlignment, 0.0)
						second_best_alignment = alignment
						alignment.mate_alignment = mateAlignment
						second_best_proper_pair = alignment.is_proper
						if alignment.molecule_id != -1 {
							alt_mol := candidate_molecules[alignment.molecule_id]
							second_best_molecule_confidence = alt_mol.molecule_confidence
							second_best_molecule_reads = alt_mol.active_alignments.Len()
						}
					}
				}
			}
		}
		// store meta data for use in determining why a read got a certain mapq.
		debug_strings := map[int]map[int]string{}
		for _, alignment := range alignmentArray {
			if alignment.active {
				alignment.mapq_data.second_best = second_best_alignment
				alignment.mapq_data.second_best_score = second_best_raw_score
				alignment.mapq_data.second_best_proper_pair = second_best_proper_pair
				alignment.mapq_data.second_best_molecule_confidence = second_best_molecule_confidence
				alignment.mapq_data.second_best_molecule_reads = second_best_molecule_reads
				alignment.mapq_data.copies = len(alignmentArray)
				alignment.mapq_data.second_best_molecule_confidence = second_best_molecule_confidence
				alignment.mapq_data.copies_in_active_molecules = read_copies_in_active_molecule[alignment.read_id]
				alignment.mapq_data.copies_outside_active_molecules = read_copies_not_in_active_molecule[read_id]
				alignment.mapq_data.unique_molecules_active = len(unique_molecules_active[read_id])
				alignment.mapq_data.score = scoreAlignment(alignment, alignment.mate_alignment, 0.0) // for the purposes of the AS bam tag, want pair alignment score without molecule penalties
				debugStrings(alignment, alignments, candidate_molecules, debug_strings, log_unpaired_probability)
			}
		}

		//sort scores and limit analysis to top 15 scoring alignment pairs
		//fmt.Println(scores)
		sort.Float64s(scores)
		for i := len(scores) - 1; i >= 0 && len(scores)-i <= 15; i-- {
			total_probability += math.Pow(10, scores[i])
		}

		// calculate mapq
		for _, alignment := range alignmentArray {

			score := scoreAlignment(alignment, alignment.mate_alignment, log_molecule_penalty)
			mapq := -10.0 * math.Log10(1.0-math.Pow(10.0, score)/total_probability)             // method 1: read probability normalization w/ molecule penalties
			moleculeMapq := -10.0 * math.Log10(1.0-(1.0/alignment.sum_move_probability_change)) // method 2: molecule move probability normalization
			mapq = math.Min(mapq, moleculeMapq)                                                 // take min of both techniques
			mapq = math.Min(float64(60.0), mapq)                                                // cap at q60
			centromereRegion, ok := centromeres[alignment.contig]
			start := -1
			end := -1
			if ok {
				start = centromereRegion.start
				end = centromereRegion.end
			}
			if alignment.pos > int64(start) && alignment.pos <= int64(end) {
				mapq = 0.0
			}
			alignment.mapq = int(mapq)
		}
	}
	checkMates(alignments)
}

func debugStrings(alignment *Alignment, alignments [][]*Alignment, candidate_molecules []*CandidateMolecule, debug_strings map[int]map[int]string, log_unpaired_probability float64) {
	if *DEBUG {
		alt_alignments := alignments[alignment.read_id]
		for _, alignment_alt := range alt_alignments {
			if alignment_alt.molecule_id != -1 {
				chrom := alignment_alt.contig
				start := candidate_molecules[alignment_alt.molecule_id].start
				end := candidate_molecules[alignment_alt.molecule_id].stop
				sinksource := 0
				sourcesink := 0
				molstring := ""
				_, has := debug_strings[alignment.molecule_id]
				if has {
					_, has = debug_strings[alignment.molecule_id][alignment_alt.molecule_id]
				}
				if !has {
					for _, read1 := range candidate_molecules[alignment.molecule_id].active_alignments.Iter() {
						read := candidate_molecules[alignment_alt.molecule_id].best_alignment_for_read.Get(read1.read_id)
						if read != nil {
							sourcesink++
						}
					}
					for _, rid := range candidate_molecules[alignment_alt.molecule_id].active_alignments.Iter() {
						has := FixGetForTypeAlignment(candidate_molecules[alignment.molecule_id].best_alignment_for_read.Get(rid.read_id)) != nil

						if has {
							sinksource++
						}
					}

					ST := strconv.FormatInt(int64(sourcesink), 10)
					TS := strconv.FormatInt(int64(sinksource), 10)
					sourcesinkchange, _ := fastScore(candidate_molecules[alignment.molecule_id], candidate_molecules[alignment_alt.molecule_id], log_unpaired_probability)
					sinksourcechange, _ := fastScore(candidate_molecules[alignment_alt.molecule_id], candidate_molecules[alignment.molecule_id], log_unpaired_probability)
					active := strconv.FormatInt(int64(candidate_molecules[alignment_alt.molecule_id].active_alignments.Len()), 10)
					spots := strconv.FormatInt(int64(candidate_molecules[alignment_alt.molecule_id].best_alignment_for_read.Len()), 10)
					STC := strconv.FormatInt(int64(sourcesinkchange), 10)
					TSC := strconv.FormatInt(int64(sinksourcechange), 10)
					molstring = chrom + ":" + strconv.FormatInt(start, 10) + "-" + strconv.FormatInt(end, 10) + ":alignments:" + active + ":spots:" + spots + ":mv_S->T:" + ST + ":" + STC + ":mv_T->S:" + TS + ":" + TSC + ","
					_, has_key := debug_strings[alignment.molecule_id]
					if !has_key {
						debug_strings[alignment.molecule_id] = map[int]string{}
					}
					debug_strings[alignment.molecule_id][alignment_alt.molecule_id] = molstring
				} else {
					molstring = debug_strings[alignment.molecule_id][alignment_alt.molecule_id]
				}
				alignment.mapq_data.active_alignments_in_molecules += molstring

			}
		}
	}
}

func setMoleculeConfidences(molecules []*CandidateMolecule) {
	for i := range molecules {
		molecules[i].molecule_confidence = moleculeConfidence(molecules[i], molecules[i].active_alignments.Len())
		for _, alignment := range molecules[i].active_alignments.Iter() {

			if alignment.soft_clipped > 0 {
				molecules[i].soft_clipped++
			}
			alignment.molecule_confidence = molecules[i].molecule_confidence
		}
	}
}

func scrapMolecules(candidate_molecules []*CandidateMolecule) []*CandidateMolecule {
	toReturn := []*CandidateMolecule{}
	count := 0
	for i := range candidate_molecules {
		if candidate_molecules[i].active_alignments.Len() > 0 {
			toReturn = append(toReturn, candidate_molecules[i])
			for _, read_id := range candidate_molecules[i].alignments.IterKeys() {
				alignment_map := FixGetForTypeOrderedMap(candidate_molecules[i].alignments.Get(read_id))
				for _, alignment_id := range alignment_map.IterKeys() {
					alignment := FixGetForTypeAlignment(alignment_map.Get(alignment_id))
					alignment.molecule_id = count
				}
			}
			count++
		} else {
			for _, read_id := range candidate_molecules[i].alignments.IterKeys() {
				alignment_map := FixGetForTypeOrderedMap(candidate_molecules[i].alignments.Get(read_id))
				for _, alignment_id := range alignment_map.IterKeys() {
					alignment := FixGetForTypeAlignment(alignment_map.Get(alignment_id))
					alignment.molecule_id = -1
				}
			}
		}
	}
	return toReturn
}

func (o Optimizer) GenerateMove(accept_move func(p_curr float64, p_next float64) bool) optimizer.Optimizable {
	sourceMolecule := o.candidate_molecules[o.currentMoleculeMoveSource]

	if sourceMolecule.active_alignments.Len() == 0 {
		o.currentMoleculeMoveSource = (o.currentMoleculeMoveSource + 1) % len(o.candidate_molecules)
		return o
	}
	var sinkMolecule *CandidateMolecule
	best_move := Move{score_change: -math.MaxFloat64}

	for i := 0; i < len(o.candidate_molecules); i++ {
		if i == o.currentMoleculeMoveSource {
			continue
		}
		sinkMolecule = o.candidate_molecules[i]

		score, move := o.fastScore(sourceMolecule, sinkMolecule)

		if (score > best_move.score_change ||
			(score == best_move.score_change && move.sink.active_alignments.Len() > best_move.sink.active_alignments.Len())) && move.num_moved > 0 {
			best_move = move
		}
	}

	best_score := best_move.score_change

	if best_score > 0 || (best_score == 0 && best_move.sink.active_alignments.Len() > sourceMolecule.active_alignments.Len()) {
		acceptMove(best_move)
	}

	o.currentMoleculeMoveSource = (o.currentMoleculeMoveSource + 1) % len(o.candidate_molecules)
	return o
}

type Move struct {
	score_change     float64
	alignment_change float64
	source           *CandidateMolecule
	sink             *CandidateMolecule
	toDelete         []int
	toSet            []*Alignment
	num_moved        int
}

func fastScore(sourceMolecule, sinkMolecule *CandidateMolecule, log_unpaired_probability float64) (float64, Move) {
	//initialization
	change := float64(0)
	alignment_change := float64(0)
	num := 0
	toDelete := []int{}
	sourceMismatchRemoveCount := map[int]int{}
	sinkMismatchAddCount := map[int]int{}
	toSet := []*Alignment{}
	soft_clipped := 0
	if *debugPrintMove {
		fmt.Println("test move ", sourceMolecule.id, " to ", sinkMolecule.id, sourceMolecule.start, sinkMolecule.start, "current alignments", sourceMolecule.active_alignments.Len(), sinkMolecule.active_alignments.Len())
	}
	if *debugPrintMove {
		fmt.Println("  source mol mismatch locs ", sourceMolecule.mismatchLocs)
	}
	if *debugPrintMove {
		fmt.Println("  sink mol mismatch locs ", sinkMolecule.mismatchLocs)
	}
	for _, sourceAlignment := range sourceMolecule.active_alignments.Iter() {
		if sourceAlignment.soft_clipped > 0 {
			soft_clipped += 1
		}
		read_id := sourceAlignment.read_id
		sinkAlignment := sinkMolecule.best_alignment_for_read.Get(read_id)
		// that have alternate alignments in the sink molecule
		if sinkAlignment != nil {
			//check if the mate pairing info is changing in this move
			mate_id := sourceAlignment.mate_id
			sourceMate := sourceMolecule.active_alignments.Get(mate_id)
			source_has_mate := sourceMate != nil
			source_has_mate_pair := source_has_mate && isPair(sourceAlignment, sourceMate)
			mate := sinkMolecule.best_alignment_for_read.Get(mate_id)
			sink_has_mate_pair := mate != nil && isPair(sinkAlignment, mate) && source_has_mate

			if !source_has_mate_pair || (source_has_mate && sink_has_mate_pair) {
				toDelete = append(toDelete, read_id)
				toSet = append(toSet, sinkAlignment)
			}
			alignment_change += sinkAlignment.log_alignment_probability - sourceAlignment.log_alignment_probability
			if *debugPrintMove {
				fmt.Println("\talignment ", sourceAlignment.pos, " to ", sinkAlignment.pos, " change score ", sinkAlignment.updated_log_alignment_probability-sourceAlignment.updated_log_alignment_probability)
				fmt.Println("\t\tsource mismatches ", sourceAlignment.mismatchLocs)
				fmt.Println("\t\tsink mismatches ", sinkAlignment.mismatchLocs)
			}

			for _, mismatchLoc := range sourceAlignment.mismatchLocs {
				numMismatch, has := sourceMolecule.mismatchLocs[mismatchLoc]
				if !has || numMismatch == 0 {
					//there is a problem
					panic(fmt.Sprint("source molecule should have this entry", has, numMismatch, sourceAlignment.contig, sourceAlignment.pos, sourceAlignment.id, sourceAlignment.read_id, mismatchLoc))
				}
				_, has = sourceMismatchRemoveCount[mismatchLoc]
				sourceMismatchRemoveCount[mismatchLoc]++
				if numMismatch-sourceMismatchRemoveCount[mismatchLoc] == 0 {
					//alignment_change += 2.0
					if *debugPrintMove {
						fmt.Println("\t\tmismatch in source alignment ", mismatchLoc, "being removed")
					}
				}
			}
			for _, mismatchLoc := range sinkAlignment.mismatchLocs {
				numMismatch, _ := sinkMolecule.mismatchLocs[mismatchLoc]
				toAdd, _ := sinkMismatchAddCount[mismatchLoc]
				if numMismatch == 0 && toAdd == 0 {
					//alignment_change -= 2.0
					if *debugPrintMove {
						fmt.Println("\t\tmismatch in sink alignment ", mismatchLoc, "being added")
					}
				}
				sinkMismatchAddCount[mismatchLoc]++
			}
			if source_has_mate_pair && !sink_has_mate_pair && sourceMolecule.id != sinkMolecule.id {

				alignment_change += log_unpaired_probability / 2.0
				if *debugPrintMove {
					fmt.Println("\t\tsource was paired and sink isnt so adding ", log_unpaired_probability/2.0)
				}
			} else if !source_has_mate_pair && sink_has_mate_pair && sourceMolecule.id != sinkMolecule.id {
				alignment_change -= log_unpaired_probability / 2.0
				if *debugPrintMove {
					fmt.Println("\t\tsink is paired and source wasnt so adding ", -log_unpaired_probability/2.0)
				}
			}
			num++
		}
	}

	source_active_before := isActiveMolecule(sourceMolecule, 0)
	source_active_after := isActiveMolecule(sourceMolecule, -num)
	if !source_active_after && source_active_before && sourceMolecule.id != sinkMolecule.id {
		change -= float64(sourceMolecule.best_alignment_for_read.Len()) * -0.5
		if *debugPrintMove {
			fmt.Println(">>> source killed adding ", -float64(sourceMolecule.best_alignment_for_read.Len())*-0.5)
		}
	}
	sink_active_before := isActiveMolecule(sinkMolecule, 0)
	sink_active_after := isActiveMolecule(sinkMolecule, num)
	if sink_active_after && !sink_active_before && sourceMolecule.id != sinkMolecule.id {
		change += float64(sinkMolecule.best_alignment_for_read.Len()) * -0.5
		if *debugPrintMove {
			fmt.Println(">>> sink created adding ", float64(sinkMolecule.best_alignment_for_read.Len())*-0.5)
		}
	}
	if sourceMolecule.active_alignments.Len()-num == 0 && num > 0 && sourceMolecule.id != sinkMolecule.id {
		change -= -3.0
		if *debugPrintMove {
			fmt.Println(">>>>>> adding 3")
		}
	}
	if sinkMolecule.active_alignments.Len() == 0 && num > 0 && sourceMolecule.id != sinkMolecule.id {
		change += -3.0
		if *debugPrintMove {
			fmt.Println(">>>>>> adding -3")
		}
	}
	change += alignment_change
	if *debugPrintMove {
		fmt.Println("\t======= final alignment change ", alignment_change)
	}
	if *debugPrintMove {
		fmt.Println("&&&&&&& final change ", change)
	}
	return change, Move{source: sourceMolecule, sink: sinkMolecule, toDelete: toDelete, toSet: toSet, num_moved: num, score_change: change, alignment_change: alignment_change}
}

func isActiveMolecule(mol *CandidateMolecule, read_change int) bool {
	active := float64(mol.active_alignments.Len() + read_change)
	potential := float64(mol.best_alignment_for_read.Len())
	if active <= 4 {
		return false
	}
	if active/potential < 0.1 {
		return false
	}
	return true
}

func moleculeConfidence(mol *CandidateMolecule, num_active int) float64 {
	density := float64(num_active) / float64(mol.best_alignment_for_read.Len())
	return density
}

func inferMolecules(positions [][]*Alignment) []*CandidateMolecule {
	toReturn := []*CandidateMolecule{}
	molecule_num := 0
	var currentMolecule *CandidateMolecule
	for _, position_list := range positions {
		for i, pos := range position_list {
			if i == 0 || (i > 0 && pos.pos-position_list[i-1].pos > *inferDistance) {
				if i > 0 {
					currentMolecule.stop = position_list[i-1].pos
				}
				currentMolecule = &CandidateMolecule{
					chrom:               pos.contig,
					start:               pos.pos,
					id:                  molecule_num,
					alignments:          NewOrderedMap(),
					molecule_confidence: 1.0,
					mismatchLocs:        map[int]int{},
				}
				aln_map := NewOrderedMap()
				aln_map.Set(pos.id, pos)
				currentMolecule.alignments.Set(pos.read_id, aln_map)
				toReturn = append(toReturn, currentMolecule)
				molecule_num++
			}
			alignment_map := FixGetForTypeOrderedMap(currentMolecule.alignments.Get(pos.read_id))
			if alignment_map != nil {
				alignment_map.Set(pos.id, pos)
			} else {
				aln_map := NewOrderedMap()
				aln_map.Set(pos.id, pos)
				currentMolecule.alignments.Set(pos.read_id, aln_map)
			}
		}
		if len(position_list) > 0 {
			currentMolecule.stop = position_list[len(position_list)-1].pos
		}
	}
	return toReturn
}

func markBestAlignmentForReadInMolecule(molecules []*CandidateMolecule) {
	active_alignment_num := 0
	for _, molecule := range molecules {
		active_alignments := NewOrderedAlignmentMap()
		best_alignment_for_read := NewOrderedAlignmentMap()
		for _, read_id := range molecule.alignments.IterKeys() {
			alignments := FixGetForTypeOrderedMap(molecule.alignments.Get(read_id))
			best_score := -math.MaxFloat64
			var best_alignment *Alignment
			for _, alignment_id := range alignments.IterKeys() {
				alignment := FixGetForTypeAlignment(alignments.Get(alignment_id))
				mate_alignments := FixGetForTypeOrderedMap(molecule.alignments.Get(alignment.mate_id))
				if mate_alignments != nil && mate_alignments.Len() > 0 {
					score := float64(0)
					// loop over mates to get best pair
					for _, mate_alignment_id := range mate_alignments.IterKeys() {
						mate_alignment := FixGetForTypeAlignment(mate_alignments.Get(mate_alignment_id))
						score = scoreAlignment(alignment, mate_alignment, 0.0)
						if score > best_score {
							best_score = score
							best_alignment = alignment
						}
					}
				} else {
					if alignment.log_alignment_probability > best_score {
						best_score = alignment.log_alignment_probability
						best_alignment = alignment
					}
				}
				if alignment.active {
					active_alignments.Set(read_id, alignment)
					active_alignment_num++
				}
			}
			if best_alignment.active {
				active_alignments.Set(read_id, best_alignment)
			}
			best_alignment_for_read.Set(read_id, best_alignment)
		}
		for _, aln := range active_alignments.Iter() {
			for _, mismatchLoc := range aln.mismatchLocs {
				_, has := molecule.mismatchLocs[mismatchLoc]
				if has {
					molecule.mismatchLocs[mismatchLoc]++
				} else {
					molecule.mismatchLocs[mismatchLoc] = 1
				}
			}
		}
		molecule.active_alignments = active_alignments
		molecule.best_alignment_for_read = best_alignment_for_read
	}
}

// alignments sent back sorted by position
// func tagBestAlignments(alignments [][]*Alignment, improper_pair_penalty float64) [][]*Alignment {
func tagBestAlignments(alignments [][]*Alignment) [][]*Alignment {
	positions := [][]*Alignment{}
	contigs := map[string]int{}
	read_ids_touched := make([]bool, len(alignments))
	improper_pair := 0
	proper_pair := 0
	reads := 0

	for read_id := range alignments {
		alignmentArray := alignments[read_id]
		touched := read_ids_touched[read_id]
		bestScore := -math.MaxFloat64
		first := true
		var bestAlignment *Alignment
		var bestMate *Alignment
		seed := int64(1)
		if len(alignmentArray) > 0 {
			md5sum := md5.Sum([]byte(*alignmentArray[0].read_name))
			seed = int64(binary.LittleEndian.Uint64(md5sum[0:8]))
		}
		random := rand.New(rand.NewSource(seed))
		reads++
		for _, alignment := range alignmentArray {
			if first {
				first = false
			}
			if read_id != alignment.read_id {
				panic(fmt.Sprintf("ids not idey%d%d", read_id, alignment.read_id))
			}
			mateAlignments := alignments[alignment.mate_id]
			totalScore := float64(0)
			for _, mateAlignment := range mateAlignments {
				//must consider the option of them separate vs them being together
				totalScore = scoreAlignment(alignment, mateAlignment, 0.0) + (random.Float64() / 2.0)
				if alignment.mate_id != mateAlignment.read_id {
					panic(fmt.Sprintf("mates not matey%d%d", alignment.mate_id, mateAlignment.read_id))
				}
				if totalScore > bestScore {
					bestScore = totalScore
					bestAlignment = alignment
					bestMate = mateAlignment
				}
			}
			if len(mateAlignments) == 0 {
				score := float64(alignment.score) + random.Float64()/2.0
				if score > bestScore {
					bestScore = score
					bestAlignment = alignment
				}
			}
			index, has_chrom := contigs[alignment.contig]

			if has_chrom {
				positions[index] = append(positions[index], alignment)
			} else {
				contigs[alignment.contig] = len(positions)
				positions = append(positions, []*Alignment{alignment})
			}
		}
		if !touched {
			bestAlignment.active = true
			bestAlignment.bwa_pick = true

			if bestMate != nil {
				if !isPair(bestAlignment, bestMate) {
					improper_pair += 2
				} else {
					proper_pair += 2
					bestAlignment.is_proper = true
					bestMate.is_proper = true
				}
				bestMate.active = true
				bestMate.bwa_pick = true
				read_ids_touched[bestMate.read_id] = true
			} else {
				improper_pair += 1
			}
		}
	}
	for _, position := range positions {
		sort.Sort(ByPosition(position))
	}
	return positions
}

// returns a map from read id to a map of
func GetAlignments(ref *gobwa.GoBwaReference, settings *gobwa.GoBwaSettings, barcode_chains [][]ChainedHit, delta int, arena *gobwa.Arena) ([][]*Alignment, [][]*Alignment) {

	toReturn := make([][]*Alignment, len(barcode_chains))
	full := make([][]*Alignment, len(barcode_chains))
	for i := range barcode_chains {
		bestScore := 0
		for _, chain := range barcode_chains[i] {
			if chain.score > bestScore {
				bestScore = chain.score
			}
		}
		// replace barcode_chains[i][j] by calling it as `chain` at the for loop
		for _, chain := range barcode_chains[i] {
			var alignment gobwa.SingleReadAlignment
			if chain.chain != nil {
				alignment = gobwa.GoBwaSmithWaterman(ref, settings, string(*(chain.read)), chain.chain, arena)
			} else {
				alignment = gobwa.SingleReadAlignment{}
			}

			matches := 0
			indels := 0
			indel_length := 0
			soft_clipping := 0
			soft_clipping_num := 0
			soft_clipping_length := 0
			refStart := chain.pos
			refEnd := chain.aend
			if alignment.Reversed {
				refStart = chain.aend + 1
				refEnd = chain.pos + 1
			}
			mismatchLocs := []int{}
			mismatchReadLocs := []int{}
			refSeq := ref.GetSeq(alignment.Chrom, refStart, refEnd, alignment.Reversed)
			refSeqOffset := 0
			readOffset := 0
			readSeq := *chain.read
			cigarStart := 0
			cigarIncrement := 2
			if alignment.Reversed {
				cigarStart = len(alignment.Cigar) - 2
				cigarIncrement = -2
			}
			for k := cigarStart; k < len(alignment.Cigar) && k >= 0; k += cigarIncrement {
				op := alignment.Cigar[k]
				opLen := int(alignment.Cigar[k+1])

				switch op {
				case 0:
					matches += opLen
					for match := range opLen {
						if refSeqOffset+match >= len(refSeq) {
							continue
						}
						if readOffset+match >= len(readSeq) {
							panic(fmt.Sprint("cigar string represents sequence larger than read?", len(readSeq), alignment.Cigar))
						}
						if refSeqOffset+match < len(refSeq) && readOffset+match < len(readSeq) && refSeq[refSeqOffset+match] != readSeq[readOffset+match] {
							if alignment.Reversed {
								mismatchLocs = append(mismatchLocs, int(refEnd)-(refSeqOffset+match))
							} else {
								mismatchLocs = append(mismatchLocs, refSeqOffset+int(refStart)+match)
							}
							mismatchReadLocs = append(mismatchReadLocs, readOffset+match)
						}
					}
					refSeqOffset += opLen
					readOffset += opLen
				case 1:
					indels += 1
					indel_length += opLen
					readOffset += opLen
				case 2:
					indels += 1
					indel_length += opLen
					refSeqOffset += opLen
				case 3:
					soft_clipping += 1
					soft_clipping_num += 1
					soft_clipping_length += opLen
					readOffset += opLen
				}
			}
			mismatches := alignment.EditDistance - indel_length
			matches -= mismatches
			if mismatches < 0 {
				mismatches = 0
			}

			var quals *[]byte
			if chain.read1 {
				quals = &chain.fastq.ReadQual1
			} else {
				quals = &chain.fastq.ReadQual2
			}
			pos := chain.pos
			aend := chain.aend
			if pos != -1 && alignment.Reversed {
				pos = chain.aend + 1
				aend = chain.pos + 1
			}

			//trim_seq := &chain.fastq.TrimBases
			//trim_qual := &chain.fastq.TrimQuals

			full_alignment := Alignment{
				id:                          chain.hit_id,
				comments:                    &chain.fastq.Tags,
				aend:                        aend,
				read_name:                   &chain.fastq.ReadInfo,
				read_seq:                    chain.read,
				read_qual:                   quals,
				matches:                     matches,
				mismatches:                  mismatches,
				mismatchLocs:                mismatchLocs,
				mismatchReadLocs:            mismatchReadLocs,
				indels:                      indels,
				soft_clipped:                soft_clipping,
				soft_clipped_length:         soft_clipping_length,
				read1:                       chain.read1,
				mapq_data:                   &MapQData{active_alignments_in_molecules: ""},
				barcode:                     &chain.fastq.Barcode,
				contig:                      alignment.Chrom,
				pos:                         pos,
				molecule_id:                 -1,
				score:                       chain.score,
				cigar:                       alignment.Cigar,
				read_id:                     chain.read_id,
				mate_id:                     chain.mate_id,
				reversed:                    alignment.Reversed,
				read_group:                  &chain.fastq.ReadGroupId,
				sum_move_probability_change: 1.0,
				molecule_confidence:         0.00001875, //0.00075 * 0.025
				duplicate:                   false,
			}

			full_alignment.log_alignment_probability = scoreAlignment(&full_alignment, nil, 0.0) - *improper_pair_penalty //remove improper pair penalty
			full_alignment.updated_log_alignment_probability = full_alignment.log_alignment_probability + 2.0*float64(len(mismatchLocs))
			if chain.aln != nil {
				full_alignment.readmap_s = chain.aln.ReadS
				full_alignment.readmap_e = chain.aln.ReadE
			}
			full[chain.read_id] = append(full[chain.read_id], &full_alignment)
			if full_alignment.score >= bestScore-delta {
				toReturn[chain.read_id] = append(toReturn[chain.read_id], &full_alignment)
			}
		}
	}
	return toReturn, full
}

func GetChains(ref *gobwa.GoBwaReference, settings *gobwa.GoBwaSettings, reads_for_barcode []fastqreader.FastQRecord, arena *gobwa.Arena, score_delta int) ([][]ChainedHit, string) {
	toReturn := [][]ChainedHit{}
	hit_num := 0
	var barcode string
	for i := range reads_for_barcode {
		read1_chains, read2_chains := gobwa.GoBwaMemMateSW(ref, settings, &reads_for_barcode[i].Read1, &reads_for_barcode[i].Read2, arena, score_delta)
		barcode = string(reads_for_barcode[i].Barcode)
		read1_num := 0
		toReturn = append(toReturn, []ChainedHit{})
		for j := range read1_chains {
			read1_chain_n := ChainedHit{
				contig:    read1_chains[j].Contig,
				pos:       read1_chains[j].Offset,
				aend:      read1_chains[j].Alignment_end,
				read_id:   i * 2,
				mate_id:   i*2 + 1,
				hit_id:    hit_num,
				read1:     true,
				secondary: read1_chains[j].Secondary,
				score:     read1_chains[j].Score,
				chain:     unsafe.Pointer(read1_chains[j].ChainedHit),
				fastq:     &reads_for_barcode[i],
				read:      &reads_for_barcode[i].Read1,
				aln:       &read1_chains[j],
			}
			read1_num++
			toReturn[len(toReturn)-1] = append(toReturn[len(toReturn)-1], read1_chain_n)
			hit_num++
		}
		if read1_num == 0 {
			toReturn[len(toReturn)-1] = append(toReturn[len(toReturn)-1], ChainedHit{
				read_id: i * 2,
				mate_id: i*2 + 1,
				pos:     -1,
				read1:   true,
				chain:   nil,
				fastq:   &reads_for_barcode[i],
				read:    &reads_for_barcode[i].Read1,
			})
			hit_num++
		}
		toReturn = append(toReturn, []ChainedHit{})
		read2_num := 0
		for j := range read2_chains {
			read2_chain_n := ChainedHit{
				contig:    read2_chains[j].Contig,
				pos:       read2_chains[j].Offset,
				aend:      read2_chains[j].Alignment_end,
				read_id:   i*2 + 1,
				mate_id:   i * 2,
				hit_id:    hit_num,
				read1:     false,
				score:     read2_chains[j].Score,
				chain:     unsafe.Pointer(read2_chains[j].ChainedHit),
				secondary: read2_chains[j].Secondary,
				fastq:     &reads_for_barcode[i],
				read:      &reads_for_barcode[i].Read2,
				aln:       &read2_chains[j],
			}
			read2_num++
			toReturn[len(toReturn)-1] = append(toReturn[len(toReturn)-1], read2_chain_n)
			hit_num++
		}
		if read2_num == 0 {
			toReturn[len(toReturn)-1] = append(toReturn[len(toReturn)-1], ChainedHit{
				read_id: i*2 + 1,
				mate_id: i * 2,
				pos:     -1,
				hit_id:  hit_num,
				read1:   false,
				chain:   nil,
				fastq:   &reads_for_barcode[i],
				read:    &reads_for_barcode[i].Read2,
			})
			hit_num++
		}
	}
	return toReturn, barcode
}
