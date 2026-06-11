package aligner

import (
	"bufio"
	"log"
	"os"

	"github.com/biogo/hts/sam"
)

// NewSamWriterChannel writes records received on the returned channel to stdout
// Sends true on doneChan when the input channel is closed and all records have been flushed.
func NewSamWriterChannel(head *sam.Header, cp, buff, threads int) (chan *sam.Record, chan bool) {
	outChan := make(chan *sam.Record, cp)
	doneChan := make(chan bool)

	fh, err := os.Stdout, error(nil)
	sio := bufio.NewWriterSize(fh, buff)
	sw, err := sam.NewWriter(sio, head, sam.FlagDecimal)
	if err != nil {
		log.Fatalf("%v", err)
	}
	go func() {
		for rec := range outChan {
			err = sw.Write(rec)
			if err != nil {
				log.Fatalf("%v", err)
			}
		}
		err = sio.Flush() // always flush
		if err != nil {
			log.Fatalf("%v", err)
		}
		doneChan <- true
	}()

	return outChan, doneChan
}

//func flushToChannel(alignments [][]*Alignment, out chan *sam.Record, contigs map[string]*sam.Reference, debugTags bool, attach_bx bool) {
//	for _, alignmentArray := range alignments {
//		if len(alignmentArray) == 0 {
//			panic("not all read_ids are spoken for")
//		}
//		read_output := false
//		for _, alignment := range alignmentArray {
//			if alignment.active {
//				out <- buildRecord(alignment, alignment, debugTags, attach_bx, contigs)
//				if alignment.secondary != nil {
//					out <- buildRecord(alignment.secondary, alignment, debugTags, attach_bx, contigs)
//				}
//				read_output = true
//			}
//		}
//		if !read_output {
//			panic("read_id has no active alignment but more than one alignment")
//		}
//	}
//}
