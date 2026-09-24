package preprocess

import (
	"github.com/biogo/hts/sam"
	"github.com/shenwei356/xopen"
)

const _mark_prefix byte = '@'
const _mark_newline byte = '\n'
const _mark_tab byte = '\t'
const _mark_plus byte = '+'

var _mark_forward = []byte{'/', '1'}
var _mark_reverse = []byte{'/', '2'}

// Write a SAM record as a FASTQ one.
func Sam2FQ(outfh *xopen.Writer, record *sam.Record, fr []byte) error {
	err := outfh.WriteByte(_mark_prefix)
	outfh.WriteString(record.Name)
	outfh.Write(fr)
	for _, aux := range record.AuxFields {
		outfh.WriteByte(_mark_tab)
		outfh.WriteString(aux.String())
	}
	outfh.WriteByte(_mark_newline)
	outfh.Write(record.Seq.Expand())
	outfh.WriteByte(_mark_newline)
	outfh.WriteByte(_mark_plus)
	outfh.WriteByte(_mark_newline)
	for _, q := range record.Qual {
		// Convert the raw Phred score (e.g., 40) to its ASCII character (e.g., 'I') by adding the Phred+33 offset.
		asciiChar := byte(q) + 33
		// clamp the value to the valid range (0-93) to avoid writing non-printable characters.
		if asciiChar < 33 {
			asciiChar = 33
		} else if asciiChar > 126 {
			asciiChar = 126
		}
		err = outfh.WriteByte(asciiChar)
	}
	outfh.WriteByte(_mark_newline)
	return err
}
