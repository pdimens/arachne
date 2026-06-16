package aligner

import "fmt"

// #TODO ALIGNMENT NEEDS tO GET RID OF 10X STUFF
type Alignment struct {
	//	trim_seq                          *[]byte
	//	trim_qual                         *[]byte
	//	raw_barcode                       *[]byte
	//  barcode_qual                      *[]byte
	//  sample_index                      *[]byte
	//  sample_index_qual                 *[]byte
	id                                int
	read1                             bool
	is_proper                         bool
	soft_clipped                      int
	soft_clipped_length               int
	barcode                           *[]byte
	read_name                         *string
	read_seq                          *[]byte
	read_qual                         *[]byte
	mapq                              int
	molecule_difference               float64
	contig                            string
	pos                               int64
	aend                              int64
	score                             int
	mismatches                        int
	matches                           int
	mismatchLocs                      []int
	mismatchReadLocs                  []int
	indels                            int
	read_id                           int
	bad_molecule                      bool
	correctly_placed                  bool
	mate_id                           int
	mate_alignment                    *Alignment
	reversed                          bool
	molecule_id                       int
	cigar                             []uint32
	read_group                        *string
	active                            bool    // the selected alignment for this read
	log_alignment_probability         float64 // does not include penalty for improperly paired
	updated_log_alignment_probability float64
	bwa_pick                          bool
	mapq_data                         *MapQData
	sum_move_probability_change       float64
	molecule_confidence               float64
	active_molecule                   bool
	readmap_s                         int
	readmap_e                         int
	secondary                         *Alignment
	primary                           *Alignment
	duplicate                         bool
}

func (aln *Alignment) Print() {
	fmt.Println(
		"read ", *aln.read_name,
		"read1", aln.read1,
		"chrom", aln.contig,
		"pos", aln.pos,
		"reverse", aln.reversed,
		"mismatches", aln.mismatches,
		"indels", aln.indels,
		"soft clipped sides", aln.soft_clipped,
		"soft clipping length", aln.soft_clipped_length,
		"active_molecule", aln.active_molecule,
		"barcode", aln.barcode,
		"molecule id", aln.molecule_id,
	)
}

func (aln *Alignment) IsUnmapped() bool {
	if !aln.is_proper && aln.score-17 < 19 {
		return true
	}
	return false
}

type MapQData struct {
	copies                          int
	copies_in_active_molecules      int
	unique_molecules_active         int
	copies_outside_active_molecules int
	reads_in_molecule               int
	active_alignments_in_molecules  string
	second_best                     *Alignment
	second_best_score               float64
	score                           float64
	second_best_proper_pair         bool
	second_best_molecule_reads      int
	second_best_molecule_confidence float64
}

func isPair(read1, read2 *Alignment) bool {
	if read1.reversed == read2.reversed || read1.contig != read2.contig {
		return false
	}
	var forward, reverse *Alignment
	if read1.reversed {
		forward = read2
		reverse = read1
	} else {
		forward = read1
		reverse = read2
	}
	dist := reverse.pos - forward.pos
	//TODO dont delete, trimming code to be turned on at later date // this is for if you have a reverse read exend further left than the forward read starts due to soft clipping and random bases matching the reference by chance.
	// if dist < 0 && dist >= -35 {
	// 	if reverse.cigar[0] == uint32(3) && len(reverse.cigar) > 2 && reverse.cigar[2] == 0 && int64(reverse.cigar[3]) > -dist {
	// 		reverse.cigar[1] += uint32(-dist)
	// 		reverse.cigar[3] -= uint32(-dist)
	// 		reverse.pos += -dist
	// 	}
	// }
	// 	cigar_end := len(forward.cigar)
	// 	overhang := reverse.aend - forward.aend
	// 	if overhang < 0 {
	// 		if forward.cigar[cigar_end-2] == uint32(3) && cigar_end > 2 && int64(forward.cigar[cigar_end-3]) > -overhang {
	// 			forward.cigar[cigar_end-1] += uint32(-overhang)
	// 			forward.cigar[cigar_end-3] -= uint32(-overhang)
	// 			forward.aend -= -overhang
	// 		}
	// 	}
	return dist >= int64(-35) && dist < int64(750)
}
