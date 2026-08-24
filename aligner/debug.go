package aligner

import "strconv"

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
