// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package aligner

import (
	"unsafe"

	sam "github.com/biogo/hts/sam"
)

const (
	SAM_CIGAR_MATCH     = 0
	SAM_CIGAR_INSERT    = 1
	SAM_CIGAR_DEL       = 2
	SAM_CIGAR_SKIP      = 3
	SAM_CIGAR_SOFT_CLIP = 4
	SAM_CIGAR_HARD_CLIP = 5
)

func AuxifyString(name []byte, data []byte) []byte {
	vec := make([]byte, len(data)+3)
	vec[0] = name[0]
	vec[1] = name[1]
	vec[2] = 'Z'
	copy(vec[3:], data)
	return vec
}

func AuxifyInt(name []byte, data int) []byte {
	vec := make([]byte, 7)
	vec[0] = name[0]
	vec[1] = name[1]
	vec[2] = byte('i')
	for i := range uint(4) {
		vec[3+i] = byte(((data) >> (8 * i)) & 0xff)
	}
	return vec
}

func AuxifyFloat(name []byte, data float32) []byte {
	vec := make([]byte, 7)
	vec[0] = name[0]
	vec[1] = name[1]
	vec[2] = byte('f')
	for i := range 4 {
		vec[3+i] = *(*byte)(unsafe.Add(unsafe.Pointer(&data), i))
	}
	return vec
}

func FixCigar(in []uint32) sam.Cigar {
	count := (len(in) / 2)
	cigar := make(sam.Cigar, count)
	for i := range count {
		cigar[i] = sam.NewCigarOp(sam.CigarOpType(in[i*2]), int(in[i*2+1]))
	}
	return cigar
}

func FixQual(in []byte) []byte {
	output := make([]byte, len(in))
	for i := range in {
		output[i] = in[i] - 33
	}
	return output
}

var cigartable = [5]uint32{
	0: 0,
	1: 1,
	2: 2,
	3: 4,
	4: 5,
}

var cigarCharacter = [4]string{
	"M",
	"I",
	"D",
	"S",
}

var complement = [256]byte{
	'A': 'T',
	'a': 'T',
	'C': 'G',
	'c': 'G',
	'G': 'C',
	'g': 'C',
	'T': 'A',
	't': 'A',
	'N': 'N',
	'n': 'N',
}

func ReverseComp(seq []byte) []byte {
	toReturn := make([]byte, len(seq))
	for i := range seq {
		toReturn[i] = complement[seq[len(seq)-i-1]]
	}
	return toReturn
}

func reverseCigar(cig []uint32) []uint32 {
	toReturn := make([]uint32, len(cig))
	for i := 0; i < len(cig); i += 2 {
		toReturn[i+1] = cig[len(cig)-i-1]
		toReturn[i] = cig[len(cig)-i-2]
	}
	return toReturn
}

func ReverseQual(qual []byte) []byte {
	toReturn := make([]byte, len(qual))
	for i := range qual {
		toReturn[i] = qual[len(qual)-i-1]
	}
	return toReturn
}

// Convert from "soft" clipping to "hard" clipping. Truncate the sequence and quality
// and convert "S" to "H" in the cigar string.
func HardClip(seq []byte, qual []byte, cigar []uint32, reversed bool) ([]byte, []byte, []uint32) {
	var start int
	end := len(seq)

	newcigar := make([]uint32, len(cigar))
	copy(newcigar, cigar)
	if len(newcigar) >= 2 {
		if newcigar[0] == SAM_CIGAR_SOFT_CLIP {
			start = int(newcigar[1])
			newcigar[0] = SAM_CIGAR_HARD_CLIP
		}
	}
	if len(newcigar) >= 4 {
		p := len(newcigar) - 2
		if newcigar[p] == SAM_CIGAR_SOFT_CLIP {
			end -= int(newcigar[p+1])
			newcigar[p] = SAM_CIGAR_HARD_CLIP
		}
	}

	newseq := seq[start:end]
	newqual := qual[start:end]
	return newseq, newqual, newcigar
}
