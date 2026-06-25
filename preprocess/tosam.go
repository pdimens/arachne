package preprocess

import (
	"bufio"

	"github.com/shenwei356/bio/seqio/fastx"
)

const fw string = "77"
const rv string = "141"
const _mark_0 byte = '0'
const _mark_ast byte = '*'

// convert the fastq record to SAM and write to buffered write
func Fq2Sam(fq *fastx.Record, fr string, w *bufio.Writer) error {
	var err error
	_, err = w.Write(fq.ID)
	w.WriteByte(_mark_tab)
	w.WriteString(fr)
	w.WriteByte(_mark_tab)
	w.WriteByte(_mark_ast)
	w.WriteByte(_mark_tab)
	w.WriteByte(_mark_0)
	w.WriteByte(_mark_tab)
	w.WriteByte(_mark_0)
	w.WriteByte(_mark_tab)
	w.WriteByte(_mark_ast)
	w.WriteByte(_mark_tab)
	w.WriteByte(_mark_ast)
	w.WriteByte(_mark_tab)
	w.WriteByte(_mark_0)
	w.WriteByte(_mark_tab)
	w.WriteByte(_mark_0)
	w.WriteByte(_mark_tab)
	w.Write(fq.Seq.Seq)
	w.WriteByte(_mark_tab)
	w.Write(fq.Seq.Qual)
	w.WriteByte(_mark_tab)
	w.Write(fq.Desc)
	w.WriteByte(_mark_newline)
	return err
}
