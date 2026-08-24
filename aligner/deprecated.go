package aligner

import (
	"fmt"
	"strings"
)

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

/* old implementation
func estimateMapQualities(alignments [][]*Alignment, candidate_molecules []*CandidateMolecule, log_unpaired_probability float64) {
	read_copies_in_active_molecule := map[int]int{}     //TODO remove, book keeping
	read_copies_not_in_active_molecule := map[int]int{} //TODO remove, book keeping
	unique_molecules_active := map[int]map[int]bool{}
	// mapq strategy 2: sum probabilities of full molecule moves
	if *debugPrintMove {
		fmt.Fprintln(os.Stderr, "NOW TESTING MAPQS")
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
		//debug_strings := map[int]map[int]string{}
		arraylen := len(alignmentArray)
		uniqmolactive := len(unique_molecules_active[read_id])
		for _, alignment := range alignmentArray {
			if alignment.active {
				alignment.mapq_data.second_best = second_best_alignment
				alignment.mapq_data.second_best_score = second_best_raw_score
				alignment.mapq_data.second_best_proper_pair = second_best_proper_pair
				alignment.mapq_data.second_best_molecule_confidence = second_best_molecule_confidence
				alignment.mapq_data.second_best_molecule_reads = second_best_molecule_reads
				alignment.mapq_data.copies = arraylen
				alignment.mapq_data.second_best_molecule_confidence = second_best_molecule_confidence
				alignment.mapq_data.copies_in_active_molecules = read_copies_in_active_molecule[alignment.read_id]
				alignment.mapq_data.copies_outside_active_molecules = read_copies_not_in_active_molecule[read_id]
				alignment.mapq_data.unique_molecules_active = uniqmolactive
				// for the purposes of the AS bam tag, want pair alignment score without molecule penalties
				alignment.mapq_data.score = scoreAlignment(alignment, alignment.mate_alignment, 0.0)
				//debugStrings(alignment, alignments, candidate_molecules, debug_strings, log_unpaired_probability)
			}
		}

		//sort scores and limit analysis to top 15 scoring alignment pairs
		sort.Float64s(scores)
		n := len(scores)
		start := max(0, n-15)
		for _, s := range scores[start:] {
			total_probability += math.Exp(s * math.Ln10)
		}
		//old
		//for i := len(scores) - 1; i >= 0 && len(scores)-i <= 15; i-- {
		//	total_probability += math.Pow(10, scores[i])
		//}

		// calculate mapq
		for _, alignment := range alignmentArray {
			score := scoreAlignment(alignment, alignment.mate_alignment, log_molecule_penalty)
			mapq := -10.0 * math.Log10(1.0-math.Exp(score*math.Ln10)/total_probability)
			//mapq := -10.0 * math.Log10(1.0-math.Pow(10.0, score)/total_probability)             // method 1: read probability normalization w/ molecule penalties
			moleculeMapq := -10.0 * math.Log10(1.0-(1.0/alignment.sum_move_probability_change)) // method 2: molecule move probability normalization
			mapq = math.Min(math.Min(mapq, moleculeMapq), 60.0)                                 // cap at q60
			// centromeres
			start := -1
			end := -1
			if centromereRegion, ok := centromeres[alignment.contig]; ok {
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
*/
