package aligner

import (
	"fmt"

	"arachne/fastqreader"
	"arachne/gominibwa"
)

// CIGAR operation codes used throughout arachne (bwa's encoding, translated
// to SAM's by fixCigar when records are written).
const (
	cigarMatch    = 0
	cigarInsert   = 1
	cigarDelete   = 2
	cigarSoftClip = 3
)

// hitAnalysis is everything arachne derives from one minibwa hit before it
// becomes an Alignment.
type hitAnalysis struct {
	cigar            []uint32 // op,length pairs in reference order, soft clips included
	matches          int
	mismatches       int
	indels           int // number of indel events
	softClipped      int // number of soft clipped ends
	softClippedLen   int // number of soft clipped bases
	mismatchLocs     []int
	mismatchReadLocs []int
}

// arachneCigar converts minibwa's CIGAR, which covers only the aligned
// region and is in reference order, into arachne's: bwa op codes as flat
// op,length pairs in reference order, with the clipped ends added as soft
// clips. Hit query coordinates are relative to the read as sequenced, so for
// a reverse-strand hit the read's start is the right-hand end of the
// alignment.
func arachneCigar(h *gominibwa.Hit, readLen int) []uint32 {
	leftClip, rightClip := h.QueryStart, readLen-h.QueryEnd
	if h.Reversed {
		leftClip, rightClip = rightClip, leftClip
	}
	cigar := make([]uint32, 0, 2*(len(h.Cigar)+2))
	if leftClip > 0 {
		cigar = append(cigar, cigarSoftClip, uint32(leftClip))
	}
	for _, op := range h.Cigar {
		switch op.Op {
		case 'M', '=', 'X':
			cigar = append(cigar, cigarMatch, uint32(op.Len))
		case 'I':
			cigar = append(cigar, cigarInsert, uint32(op.Len))
		case 'D', 'N':
			cigar = append(cigar, cigarDelete, uint32(op.Len))
		}
	}
	if rightClip > 0 {
		cigar = append(cigar, cigarSoftClip, uint32(rightClip))
	}
	return cigar
}

// analyzeHit walks the alignment of read against refSeq, which must be the
// reference slice [h.Start, h.End) reverse complemented if the hit is on the
// reverse strand (gominibwa.Index.GetSeq does this). Mismatch positions are
// reported as forward-strand reference coordinates and read offsets in the
// read as sequenced.
func analyzeHit(h *gominibwa.Hit, read, refSeq []byte) hitAnalysis {
	a := hitAnalysis{cigar: arachneCigar(h, len(read))}
	indelLength := 0
	refSeqOffset, readOffset := 0, 0

	// refSeq and the read run left to right together, so a reverse-strand
	// alignment's reference-ordered CIGAR is walked from its end.
	start, step := 0, 2
	if h.Reversed {
		start, step = len(a.cigar)-2, -2
	}
	for k := start; k >= 0 && k < len(a.cigar); k += step {
		op := a.cigar[k]
		n := int(a.cigar[k+1])
		switch op {
		case cigarMatch:
			a.matches += n
			for i := range n {
				ri, qi := refSeqOffset+i, readOffset+i
				if ri >= len(refSeq) || qi >= len(read) || refSeq[ri] == read[qi] {
					continue
				}
				if h.Reversed {
					a.mismatchLocs = append(a.mismatchLocs, reverseStrandMismatchRefPos(h.End, refSeqOffset, i))
				} else {
					a.mismatchLocs = append(a.mismatchLocs, int(h.Start)+ri)
				}
				a.mismatchReadLocs = append(a.mismatchReadLocs, qi)
			}
			refSeqOffset += n
			readOffset += n
		case cigarInsert:
			a.indels++
			indelLength += n
			readOffset += n
		case cigarDelete:
			a.indels++
			indelLength += n
			refSeqOffset += n
		case cigarSoftClip:
			a.softClipped++
			a.softClippedLen += n
			readOffset += n
		}
	}
	a.mismatches = max(h.EditDistance-indelLength, 0)
	a.matches -= a.mismatches
	return a
}

// GetAlignments aligns every read pair of one barcode with minibwa and
// returns, per read id (2i for read 1 of pair i, 2i+1 for read 2):
//   - the alignments within delta of that read's best score, which are the
//     candidates RFA and EM choose between, and
//   - all alignments, which split read detection searches.
//
// Every read has at least one entry: a read that did not align gets a
// placeholder at position -1. It also returns the barcode of the reads.
func GetAlignments(mapper *gominibwa.Mapper, reads []fastqreader.FastQRecord, delta int) (candidates, full [][]*Alignment, barcode string) {
	hits := mapReads(mapper, reads)
	matchScore := max(mapper.Options().MatchScore(), 1)
	idx := mapper.Index()

	candidates = make([][]*Alignment, 2*len(reads))
	full = make([][]*Alignment, 2*len(reads))
	hitID := 0
	for i := range reads {
		fq := &reads[i]
		barcode = string(fq.Barcode)
		for mate := range 2 {
			readID := 2*i + mate
			read1 := mate == 0
			readSeq, quals := &fq.Read1, &fq.ReadQual1
			if !read1 {
				readSeq, quals = &fq.Read2, &fq.ReadQual2
			}
			readHits := hits[readID]

			bestScore := 0
			for k := range readHits {
				bestScore = max(bestScore, readHits[k].DPScore/matchScore)
			}
			newAlignment := func() *Alignment {
				return &Alignment{
					id:                          hitID,
					comments:                    &fq.Tags,
					read_name:                   &fq.ReadInfo,
					read_seq:                    readSeq,
					read_qual:                   quals,
					mapq_data:                   &MapQData{active_alignments_in_molecules: ""},
					barcode:                     &fq.Barcode,
					read1:                       read1,
					read_id:                     readID,
					mate_id:                     readID ^ 1,
					molecule_id:                 -1,
					read_group:                  sample_id, // matches the @RG ID written by buildHeader
					sum_move_probability_change: 1.0,
					molecule_confidence:         0.00001875, //0.00075 * 0.025
				}
			}

			if len(readHits) == 0 {
				// placeholder for a read with no alignment
				aln := newAlignment()
				aln.pos = -1
				finishAlignment(aln)
				hitID++
				candidates[readID] = append(candidates[readID], aln)
				full[readID] = append(full[readID], aln)
				continue
			}
			for k := range readHits {
				h := &readHits[k]
				var refSeq []byte
				if seq, err := idx.GetSeq(h.Tid, h.Start, h.End, h.Reversed); err == nil {
					refSeq = seq
				}
				an := analyzeHit(h, *readSeq, refSeq)

				aln := newAlignment()
				aln.contig = h.Contig
				aln.pos = h.Start
				aln.aend = h.End
				aln.reversed = h.Reversed
				aln.score = h.DPScore / matchScore
				aln.cigar = an.cigar
				aln.matches = an.matches
				aln.mismatches = an.mismatches
				aln.mismatchLocs = an.mismatchLocs
				aln.mismatchReadLocs = an.mismatchReadLocs
				aln.indels = an.indels
				aln.soft_clipped = an.softClipped
				aln.soft_clipped_length = an.softClippedLen
				aln.readmap_s = h.QueryStart
				aln.readmap_e = h.QueryEnd
				finishAlignment(aln)
				hitID++

				full[readID] = append(full[readID], aln)
				if aln.score >= bestScore-delta {
					candidates[readID] = append(candidates[readID], aln)
				}
			}
		}
	}
	return candidates, full, barcode
}

// finishAlignment sets the likelihood fields, which depend on the others.
func finishAlignment(aln *Alignment) {
	aln.log_alignment_probability = scoreAlignment(aln, nil, 0.0) - *improper_pair_penalty // remove improper pair penalty
	aln.updated_log_alignment_probability = aln.log_alignment_probability + 2.0*float64(len(aln.mismatchLocs))
}

// mapReads aligns all pairs of a barcode in one batch and returns the hits
// per read id. A pair with an empty read cannot be aligned as a pair, so its
// non-empty read, if any, is aligned alone.
func mapReads(mapper *gominibwa.Mapper, reads []fastqreader.FastQRecord) [][]gominibwa.Hit {
	hits := make([][]gominibwa.Hit, 2*len(reads))
	var names []string
	var reads1, reads2 [][]byte
	var pairIdx []int
	for i := range reads {
		fq := &reads[i]
		switch {
		case len(fq.Read1) > 0 && len(fq.Read2) > 0:
			names = append(names, fq.ReadInfo)
			reads1 = append(reads1, fq.Read1)
			reads2 = append(reads2, fq.Read2)
			pairIdx = append(pairIdx, i)
		case len(fq.Read1) > 0:
			hits[2*i] = mapper.Map(fq.ReadInfo, fq.Read1)
		case len(fq.Read2) > 0:
			hits[2*i+1] = mapper.Map(fq.ReadInfo, fq.Read2)
		}
	}
	if len(pairIdx) == 0 {
		return hits
	}
	res, err := mapper.MapPairs(names, reads1, reads2)
	if err != nil {
		panic(fmt.Sprintf("minibwa failed to align a barcode: %v", err))
	}
	for k, i := range pairIdx {
		hits[2*i] = res[k][0]
		hits[2*i+1] = res[k][1]
	}
	return hits
}
