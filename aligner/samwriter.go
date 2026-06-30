package aligner

import (
	"bytes"
	"log"
	"strconv"

	sam "github.com/biogo/hts/sam"
)

var _TAB = []byte("\t")
var auxSA = []byte("SA")
var auxAS = []byte("AS")
var auxRG = []byte("RG")
var auxBX = []byte("BX")
var auxVX = []byte("VX")
var auxXB = []byte("XB")
var auxXC = []byte("XC")
var auxXM = []byte("XM")
var auxXS = []byte("XS")
var auxAC = []byte("AC")
var auxXT = []byte("XT")
var auxAM = []byte("AM")
var auxDM = []byte("DM")

func fixCigar(in []uint32) []uint32 {
	var out = make([]uint32, len(in))
	for i := 0; i < len(in)/2; i++ {
		idx := i * 2
		if int(in[idx]) > len(cigartable) {
			log.Printf("BAMOP: %v", in[idx])
			panic("ILLEGAL CIGAR OP")
		}
		out[idx] = cigartable[int(in[idx])]
		out[idx+1] = in[idx+1]
	}
	return out
}

func buildRecord(aln, primary *Alignment, debugTags *bool, contigs map[string]*sam.Reference) *sam.Record {
	rec := &sam.Record{}
	// replace b.Contigs[...] with contigs[...]

	ref := contigs[aln.contig]
	var flags int32

	if !aln.is_proper && aln.score-17 < 19 {
		aln.pos = -1
		aln.mapq = 0
	}
	if aln.mate_id >= 0 {
		flags |= 1
		if aln.is_proper {
			if aln == primary {
				flags |= 0x2
			} else {
				if isPair(aln, primary.mate_alignment) {
					flags |= 0x2
				}
			}
		}

		if primary.mate_alignment.pos == -1 || (!primary.is_proper && primary.mate_alignment.score-17 < 19) {
			// Mate is unmapped
			flags |= 0x8
			rec.MatePos = -1
			rec.MateRef = nil
		} else {
			// Mate mapped
			if primary.mate_alignment.reversed {
				flags |= 0x20
			}
			rec.MateRef = contigs[primary.mate_alignment.contig]
			rec.MatePos = int(primary.mate_alignment.pos)
		}

		if aln.read1 {
			flags |= 0x40
		} else {
			flags |= 0x80
		}
		if aln.duplicate {
			flags |= 0x400
		}

		if primary.mate_alignment.pos == -1 {
			rec.MateRef = nil
			rec.TempLen = 0
		} else if aln == primary {
			if aln.contig == aln.mate_alignment.contig && (primary.is_proper || primary.mate_alignment.score-17 >= 19) {
				if aln.reversed {
					rec.TempLen = -int(aln.aend - aln.mate_alignment.pos)
				} else {
					rec.TempLen = int(aln.mate_alignment.aend - aln.pos)
				}
			} else {
				rec.TempLen = 0
			}
		} else {
			rec.TempLen = 0
		}
	} else {
		rec.MatePos = -1
		rec.MateRef = nil
	}

	if aln != primary {
		flags |= 256
	}

	rec.Ref = ref

	rec.MapQ = byte(aln.mapq)
	if aln.pos == -1 {
		flags |= 0x4
		rec.MapQ = byte(0)
		rec.Ref = nil
	}
	if aln.reversed {
		flags |= 0x10
	}
	rec.Name = string(*aln.read_name)
	rec.Flags = sam.Flags(flags)

	seq := *aln.read_seq
	pos := int(aln.pos)
	cigar := fixCigar(aln.cigar)
	qual := *aln.read_qual

	if aln.reversed {
		seq = ReverseComp(seq)
		qual = ReverseQual(qual)
	}

	if primary != aln {
		//var deltapos int // NOTE: never assigned — this block is currently a no-op; check HardClip's contract
		seq, qual, cigar = HardClip(seq, qual, cigar, aln.reversed)
		//if pos > 0 {
		//	pos += deltapos
		//}
	}

	rec.Pos = pos
	rec.Cigar = FixCigar(cigar)
	rec.Seq = sam.NewSeq(seq)
	rec.Qual = FixQual(qual)

	aux := []sam.Aux{}
	as := AuxifyInt(auxAS, aln.score)
	if len(*aln.read_group) > 0 {
		rg := AuxifyString(auxRG, []byte(*aln.read_group))
		aux = append(aux, sam.Aux(rg))
	}
	if aln.mapq_data != nil {
		xs := AuxifyInt(auxXS, int(aln.mapq_data.second_best_score))
		aux = append(aux, sam.Aux(xs))
		as = AuxifyInt(auxAS, int(aln.mapq_data.score))

		var xcBuf []byte
		if aln.mapq_data.second_best != nil {
			mismatchReadLocs := aln.mapq_data.second_best.mismatchReadLocs
			mismatchLocs := aln.mapq_data.second_best.mismatchLocs
			xcBuf = make([]byte, 0, len(mismatchReadLocs)*8)
			for i := range mismatchReadLocs {
				xcBuf = strconv.AppendInt(xcBuf, int64(mismatchLocs[i]), 10)
				xcBuf = append(xcBuf, ',')
				xcBuf = strconv.AppendInt(xcBuf, int64(mismatchReadLocs[i]), 10)
				xcBuf = append(xcBuf, ',', '1', ';')
			}
		}
		xc := AuxifyString(auxXC, xcBuf)
		if len(xcBuf) > 0 {
			aux = append(aux, sam.Aux(xc))
		}

		mismatchReadLocs := aln.mismatchReadLocs
		mismatchLocs := aln.mismatchLocs
		acBuf := make([]byte, 0, len(mismatchReadLocs)*8)
		for i := range mismatchReadLocs {
			acBuf = strconv.AppendInt(acBuf, int64(mismatchLocs[i]), 10)
			acBuf = append(acBuf, ',')
			acBuf = strconv.AppendInt(acBuf, int64(mismatchReadLocs[i]), 10)
			acBuf = append(acBuf, ',', '1', ';')
		}
		ac := AuxifyString(auxAC, acBuf)
		if len(acBuf) > 0 {
			aux = append(aux, sam.Aux(ac))
		}
	}
	aux = append(aux, sam.Aux(as))

	second_best_active_molecule := 0
	if aln.mapq_data != nil && aln.mapq_data.second_best != nil && aln.mapq_data.second_best.active_molecule {
		second_best_active_molecule = 1
	}
	xm := AuxifyString(auxXM, []byte(strconv.FormatInt(int64(second_best_active_molecule), 10)))
	aux = append(aux, sam.Aux(xm))

	active_molecule := "0"
	if aln.active_molecule {
		active_molecule = "1"
	}
	am := AuxifyString(auxAM, []byte(active_molecule))
	aux = append(aux, sam.Aux(am))

	tandem := 0
	if aln.mapq_data != nil && aln.mapq_data.second_best != nil && aln.molecule_id == aln.mapq_data.second_best.molecule_id {
		tandem = 1
	}
	xt := AuxifyInt(auxXT, tandem)
	aux = append(aux, sam.Aux(xt))

	var secondaryAlignment *Alignment
	if aln.secondary != nil {
		secondaryAlignment = aln.secondary
	} else if aln.primary != nil {
		secondaryAlignment = aln.primary
	}
	if secondaryAlignment != nil && secondaryAlignment.pos > -1 {
		cigarBytes := secondaryAlignment.cigar
		strandByte := byte('+')
		if secondaryAlignment.reversed {
			strandByte = '-'
			cigarBytes = reverseCigar(cigarBytes)
		}

		cigarBuf := make([]byte, 0, len(cigarBytes)*3)
		indelLength := 0
		for cig := 0; cig < len(cigarBytes); cig += 2 {
			var cigChar byte
			if cigarBytes[cig] == 3 && aln.secondary != nil {
				cigChar = 'H'
			} else {
				cigChar = cigarCharacter[cigarBytes[cig]][0]
			}
			if cigarBytes[cig] == 1 || cigarBytes[cig] == 2 {
				indelLength += int(cigarBytes[cig+1])
			}
			cigarBuf = strconv.AppendInt(cigarBuf, int64(cigarBytes[cig+1]), 10)
			cigarBuf = append(cigarBuf, cigChar)
		}

		sndBuf := make([]byte, 0, len(secondaryAlignment.contig)+len(cigarBuf)+24)
		sndBuf = append(sndBuf, secondaryAlignment.contig...)
		sndBuf = append(sndBuf, ',')
		sndBuf = strconv.AppendInt(sndBuf, int64(secondaryAlignment.pos), 10)
		sndBuf = append(sndBuf, ',', strandByte, ',')
		sndBuf = append(sndBuf, cigarBuf...)
		sndBuf = append(sndBuf, ',')
		sndBuf = strconv.AppendInt(sndBuf, int64(secondaryAlignment.mapq), 10)
		sndBuf = append(sndBuf, ',')
		sndBuf = strconv.AppendInt(sndBuf, int64(len(secondaryAlignment.mismatchLocs)+indelLength), 10)
		sndBuf = append(sndBuf, ';')

		sa := AuxifyString(auxSA, sndBuf)
		aux = append(aux, sam.Aux(sa))
	}

	if *debugTags && aln.mapq_data != nil {
		addAuxDebug(aln, primary, aux)
	}

	bx := AuxifyString(auxBX, *aln.barcode)
	aux = append(aux, sam.Aux(bx))
	vx := AuxifyInt(auxVX, 1)
	aux = append(aux, sam.Aux(vx))

	if aln.active_molecule {
		md := AuxifyString(auxDM, []byte(strconv.FormatFloat(aln.molecule_difference, 'f', 6, 64)))
		aux = append(aux, sam.Aux(md))
	}

	if *AddComments && len(*aln.comments) > 0 {
		for field := range bytes.SplitSeq(*aln.comments, _TAB) {
			a, err := sam.ParseAux(field)
			if err != nil {
				// ignore malformed field
				continue
			}
			aux = append(aux, sam.Aux(a))
		}
	}
	rec.AuxFields = aux
	return rec
}

func flushToChannel(alignments [][]*Alignment, out chan *sam.Record, contigs map[string]*sam.Reference, debugTags *bool) {
	for _, alignmentArray := range alignments {
		if len(alignmentArray) == 0 {
			panic("not all read_ids are spoken for")
		}
		read_output := false
		for _, alignment := range alignmentArray {
			if alignment.active {
				out <- buildRecord(alignment, alignment, debugTags, contigs)
				if alignment.secondary != nil {
					out <- buildRecord(alignment.secondary, alignment, debugTags, contigs)
				}
				read_output = true
			}
		}
		if !read_output {
			panic("read_id has no active alignment but more than one alignment")
		}
	}
}

// Add auxilary tags with debug information. Modifies aux in-place
func addAuxDebug(aln, primary *Alignment, aux []sam.Aux) {
	// NOTE: these statistics generally refer to the configuration of the active molecules after the
	// Arachne optimization process has finished.

	// Total number of alignments returned by BWA
	cp := AuxifyString([]byte("CP"), []byte(strconv.FormatInt(int64(aln.mapq_data.copies), 10)))
	// number of alignments in active molecules
	cm := AuxifyString([]byte("CM"), []byte(strconv.FormatInt(int64(aln.mapq_data.copies_in_active_molecules), 10)))
	// number of unique active molecules
	cu := AuxifyString([]byte("CU"), []byte(strconv.FormatInt(int64(aln.mapq_data.unique_molecules_active), 10)))
	// Alignments outside active molecules
	cs := AuxifyString([]byte("CS"), []byte(strconv.FormatInt(int64(aln.mapq_data.copies_outside_active_molecules), 10)))
	// Total number of active alignments in the molecule containing the alignment
	rd := AuxifyString([]byte("RD"), []byte(strconv.FormatInt(int64(aln.mapq_data.reads_in_molecule), 10)))
	// Alignment of the read-pair forms a 'proper' read-pair: reads have the correct relative orientation & distance.
	pp := AuxifyString([]byte("PP"), []byte(strconv.FormatBool(aln.is_proper)))
	// A string representation of the alignments for this read that fall in active molecules.
	aa := AuxifyString([]byte("AA"), []byte(aln.mapq_data.active_alignments_in_molecules))
	// Confidence score for the existence of the molecule containing this alignment
	mc := AuxifyString([]byte("MC"), []byte(strconv.FormatFloat(float64(aln.molecule_confidence), 'f', 6, 64)))
	ms := AuxifyString([]byte("MS"), []byte(strconv.FormatFloat(float64(aln.sum_move_probability_change), 'f', 6, 64)))
	// Mate alignment score
	ps := AuxifyString([]byte("PS"), []byte(strconv.FormatInt(int64(primary.mate_alignment.score), 10)))
	pl := AuxifyString([]byte("PL"), []byte(strconv.FormatFloat(float64(primary.mate_alignment.log_alignment_probability), 'f', 6, 64)))
	// Count of alignment operations in this alignment
	ac := AuxifyString([]byte("AC"), []byte("Match:"+strconv.FormatInt(int64(aln.matches), 10)+":Mismatches:"+strconv.FormatInt(int64(aln.mismatches), 10)+":Indels:"+strconv.FormatInt(int64(aln.indels), 10)+":soft_clipped:"+strconv.FormatInt(int64(aln.soft_clipped), 10)))
	// Count of alignment operations in the mate of this alignment
	pc := AuxifyString([]byte("PC"), []byte("Match:"+strconv.FormatInt(int64(primary.mate_alignment.matches), 10)+":Mismatches:"+strconv.FormatInt(int64(primary.mate_alignment.mismatches), 10)+":Indels:"+strconv.FormatInt(int64(primary.mate_alignment.indels), 10)+":soft_clipped:"+strconv.FormatInt(int64(primary.mate_alignment.soft_clipped), 10)))
	if aln.mapq_data.second_best != nil {
		second_best_log_probability := AuxifyString([]byte("XL"), []byte(strconv.FormatFloat(aln.mapq_data.second_best.log_alignment_probability, 'f', 6, 64)))
		second_best_proper_pair := AuxifyString([]byte("XP"), []byte(strconv.FormatBool(aln.mapq_data.second_best_proper_pair)))
		second_best_molecule_reads := AuxifyString([]byte("XR"), []byte(strconv.FormatInt(int64(aln.mapq_data.second_best_molecule_reads), 10)))
		second_best_molecule_confidence := AuxifyString([]byte("XC"), []byte(strconv.FormatFloat(aln.mapq_data.second_best_molecule_confidence, 'f', 6, 64)))
		if aln.mapq_data.second_best.mate_alignment != nil {
			xm := AuxifyString([]byte("XM"), []byte(strconv.FormatFloat(float64(aln.mapq_data.second_best.mate_alignment.log_alignment_probability), 'f', 6, 64)))
			aux = append(aux, sam.Aux(xm))
			xz := AuxifyString([]byte("XZ"), []byte("Match:"+strconv.FormatInt(int64(aln.mapq_data.second_best.mate_alignment.matches), 10)+":Mismatches:"+strconv.FormatInt(int64(aln.mapq_data.second_best.mate_alignment.mismatches), 10)+":Indels:"+strconv.FormatInt(int64(aln.mapq_data.second_best.mate_alignment.indels), 10)+":soft_clipped:"+strconv.FormatInt(int64(aln.mapq_data.second_best.mate_alignment.soft_clipped), 10)))
			aux = append(aux, sam.Aux(xz))
		}
		xx := AuxifyString([]byte("XX"), []byte("Match:"+strconv.FormatInt(int64(aln.mapq_data.second_best.matches), 10)+":Mismatches:"+strconv.FormatInt(int64(aln.mapq_data.second_best.mismatches), 10)+":Indels:"+strconv.FormatInt(int64(aln.mapq_data.second_best.indels), 10)+":soft_clipped:"+strconv.FormatInt(int64(aln.mapq_data.second_best.soft_clipped), 10)))
		aux = append(aux, sam.Aux(xx))
		aux = append(aux, sam.Aux(second_best_log_probability))
		aux = append(aux, sam.Aux(second_best_proper_pair))
		aux = append(aux, sam.Aux(second_best_molecule_reads))
		aux = append(aux, sam.Aux(second_best_molecule_confidence))
	}
	aux = append(aux, sam.Aux(aa))
	aux = append(aux, sam.Aux(cp))
	aux = append(aux, sam.Aux(cm))
	aux = append(aux, sam.Aux(cu))
	aux = append(aux, sam.Aux(cs))
	aux = append(aux, sam.Aux(rd))
	aux = append(aux, sam.Aux(ms))
	aux = append(aux, sam.Aux(mc))
	aux = append(aux, sam.Aux(pp))
	aux = append(aux, sam.Aux(ps))
	aux = append(aux, sam.Aux(pl))
	aux = append(aux, sam.Aux(ac))
	aux = append(aux, sam.Aux(pc))
}
