package aligner

// If two reads have the same value, then they are duplicates
type readDupTuple struct {
	read1      bool
	reversed   bool
	contig     string
	pos        int64
	mateContig string
	matePos    int64
}

// For each read, make a tuple of (bc_sequence, read.is_read1, read.is_reverse, read.tid, read.pos, read.mrnm, read.mpos)
// reads with an equal value of this tuple are defined as duplicates or one another.
// mark all but 1 read in each group as a duplicate
func markDuplicates(alignments [][]*Alignment) {

	dupSeen := make(map[readDupTuple]bool)

	//now go through every read_id and normalize all alternate alignment probabilities
	for _, alignmentArray := range alignments {
		for _, alignment := range alignmentArray {
			if alignment.active {

				mateAlignment := alignment.mate_alignment

				readTuple := readDupTuple{
					read1:      alignment.read1,
					reversed:   alignment.reversed,
					contig:     alignment.contig,
					pos:        alignment.pos,
					mateContig: mateAlignment.contig,
					matePos:    mateAlignment.pos}

				// If we have seen this tuple before, mark it as duplicate
				// Otherwise note tuple
				_, haveSeen := dupSeen[readTuple]
				if haveSeen {
					alignment.duplicate = true
				} else {
					dupSeen[readTuple] = true
				}
			}
		}
	}
}
