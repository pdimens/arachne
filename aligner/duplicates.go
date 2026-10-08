package aligner

// dupMinBaseQual is the lowest base quality that counts toward a pair's
// duplicate score. Picard and samtools markdup use the same cutoff.
const dupMinBaseQual = 15

// dupEnd is where one mate of a pair was placed.
type dupEnd struct {
	contig   string
	pos      int64
	reversed bool
}

// pairDupKey identifies a set of duplicate pairs. Pairs with equal keys are
// duplicates of one another.
//
// A pair with both mates mapped is keyed on the placement of read 1 and read 2.
// A pair with only one mapped mate (an orphan) is keyed on that mate alone, and
// on which mate it is, so orphans compete only with other orphans.
type pairDupKey struct {
	orphan      bool
	orphanRead1 bool // orphans only: the mapped mate is read 1
	first       dupEnd
	second      dupEnd
}

// dupPair is a read pair that is a candidate to be the representative of its
// duplicate set. Either mate may be nil if it is unmapped.
type dupPair struct {
	read1, read2 *Alignment
	score        int
}

func (p *dupPair) setDuplicate(dup bool) {
	if p.read1 != nil {
		p.read1.duplicate = dup
	}
	if p.read2 != nil {
		p.read2.duplicate = dup
	}
}

// activeAlignment returns the selected alignment of a read, or nil.
func activeAlignment(candidates []*Alignment) *Alignment {
	for _, a := range candidates {
		if a.active {
			return a
		}
	}
	return nil
}

// mappedForDup reports whether a read will be written as mapped. This mirrors
// the demotion in buildRecord, so a read written as unmapped is never treated
// as a placement. Unmapped reads all share pos=-1 and contig="", so counting
// them as duplicate-eligible would flag unrelated reads as duplicates.
func mappedForDup(a *Alignment) bool {
	return a != nil && a.pos != -1 && !a.IsUnmapped()
}

// baseQualitySum is the sum of the Phred scores of the bases at or above
// dupMinBaseQual. read_qual is Phred+33 ASCII, as read from the FASTQ.
func baseQualitySum(a *Alignment) int {
	if a == nil || a.read_qual == nil {
		return 0
	}
	sum := 0
	for _, q := range *a.read_qual {
		if q < 33 {
			continue
		}
		if phred := int(q) - 33; phred >= dupMinBaseQual {
			sum += phred
		}
	}
	return sum
}

// markDuplicates marks read pairs that are duplicates of one another and,
// within each set, keeps the pair with the highest base quality sum. Both mates
// of a pair always get the same flag. Ties go to the pair seen first.
//
// alignments holds one slice of candidates per read in read order: read 1 of
// pair i at index 2i and read 2 at 2i+1, as GetAlignments returns them. Only
// each read's active alignment is considered.
//
//   - Unmapped reads are never flagged. If one mate is unmapped, the other is
//     an orphan and is compared only with orphans on the same mate and placement.
//   - Split (supplementary) records are not scored. buildRecord gives them the
//     flag of their primary alignment.
func markDuplicates(alignments [][]*Alignment) {
	// init at 128 to mitigate performance hits by growing underlying container when too big
	best := make(map[pairDupKey]*dupPair, 128)

	for i := 0; i < len(alignments); i += 2 {
		read1 := activeAlignment(alignments[i])
		var read2 *Alignment
		if i+1 < len(alignments) {
			read2 = activeAlignment(alignments[i+1])
		}
		if !mappedForDup(read1) {
			read1 = nil
		}
		if !mappedForDup(read2) {
			read2 = nil
		}

		var key pairDupKey
		switch {
		case read1 != nil && read2 != nil:
			key.first = dupEnd{read1.contig, read1.pos, read1.reversed}
			key.second = dupEnd{read2.contig, read2.pos, read2.reversed}
		case read1 != nil:
			key.orphan, key.orphanRead1 = true, true
			key.first = dupEnd{read1.contig, read1.pos, read1.reversed}
		case read2 != nil:
			key.orphan = true
			key.first = dupEnd{read2.contig, read2.pos, read2.reversed}
		default:
			continue
		}

		pair := &dupPair{read1: read1, read2: read2, score: baseQualitySum(read1) + baseQualitySum(read2)}
		current, seen := best[key]
		switch {
		case !seen:
			pair.setDuplicate(false)
			best[key] = pair
		case pair.score > current.score:
			current.setDuplicate(true)
			pair.setDuplicate(false)
			best[key] = pair
		default:
			pair.setDuplicate(true)
		}
	}
}
