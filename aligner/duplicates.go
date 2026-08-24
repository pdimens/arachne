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
	// init at 128 to mitigate performance hits by growing underlying container when too big
	dupSeen := make(map[readDupTuple]struct{}, 128)

	//now go through every read_id and normalize all alternate alignment probabilities
	for _, alignmentArray := range alignments {
		for _, alignment := range alignmentArray {
			if !alignment.active {
				continue
			}
			readTuple := readDupTuple{
				read1:      alignment.read1,
				reversed:   alignment.reversed,
				contig:     alignment.contig,
				pos:        alignment.pos,
				mateContig: alignment.mate_alignment.contig,
				matePos:    alignment.mate_alignment.pos,
			}
			// If we have seen this tuple before, mark it as duplicate
			// Otherwise note tuple
			_, haveSeen := dupSeen[readTuple]
			if haveSeen {
				alignment.duplicate = true
			} else {
				dupSeen[readTuple] = struct{}{}
			}
		}
	}
}

/* alternative implementation
type dupKey struct {
	a uint64 // bit 0: read1 | bit 1: reversed | bits 2-17: contigID | bits 18-49: pos
	b uint64 // bits 0-15: mateContigID | bits 16-47: matePos
}

func markDuplicates(alignments [][]*Alignment) {
	contigIDs := make(map[string]uint16, 8) // few distinct contigs per barcode group
	var nextID uint16

	idFor := func(contig string) uint16 {
		if id, ok := contigIDs[contig]; ok {
			return id
		}
		id := nextID
		contigIDs[contig] = id
		nextID++
		return id
	}

	seen := make([]dupKey, 0, 32)

	for _, alignmentArray := range alignments {
		for _, alignment := range alignmentArray {
			if !alignment.active {
				continue
			}

			var flags uint64
			if alignment.read1 {
				flags |= 1
			}
			if alignment.reversed {
				flags |= 1 << 1
			}
			key := dupKey{
				a: flags | uint64(idFor(alignment.contig))<<2 | uint64(uint32(alignment.pos))<<18,
				b: uint64(idFor(alignment.mate_alignment.contig)) | uint64(uint32(alignment.mate_alignment.pos))<<16,
			}

			dup := false
			for _, k := range seen {
				if k == key {
					dup = true
					break
				}
			}
			if dup {
				alignment.duplicate = true
			} else {
				seen = append(seen, key)
			}
		}
	}
}
*/
