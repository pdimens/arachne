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

// Write a SAM record as a FASTQ one. Returns the first write error
// encountered, if any; once a write fails, no further writes are attempted.
func Sam2FQ(outfh *xopen.Writer, record *sam.Record, fr []byte) error {
	var err error
	write := func(b []byte) {
		if err != nil {
			return
		}
		_, err = outfh.Write(b)
	}
	writeByte := func(b byte) {
		if err != nil {
			return
		}
		err = outfh.WriteByte(b)
	}
	writeString := func(s string) {
		if err != nil {
			return
		}
		_, err = outfh.WriteString(s)
	}

	writeByte(_mark_prefix)
	writeString(record.Name)
	write(fr)
	for _, aux := range record.AuxFields {
		writeByte(_mark_tab)
		writeString(aux.String())
	}
	writeByte(_mark_newline)
	write(record.Seq.Expand())
	writeByte(_mark_newline)
	writeByte(_mark_plus)
	writeByte(_mark_newline)
	for _, q := range record.Qual {
		// Convert the raw Phred score (e.g., 40) to its ASCII character (e.g., 'I') by adding the Phred+33 offset.
		asciiChar := byte(q) + 33
		// clamp the value to the valid range (0-93) to avoid writing non-printable characters.
		if asciiChar < 33 {
			asciiChar = 33
		} else if asciiChar > 126 {
			asciiChar = 126
		}
		writeByte(asciiChar)
	}
	writeByte(_mark_newline)
	return err
}
