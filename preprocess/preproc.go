package preprocess

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

func Preprocess(threads int, prefix, r1Path, r2Path string) error {
	if err := os.MkdirAll(filepath.Dir(prefix), 0o755); err != nil {
		return fmt.Errorf("creating output directory: %w", err)
	}
	// determine what kind of linked-read tech it is
	processBC, err := CheckFastqFormat(r1Path)
	if err != nil {
		return fmt.Errorf("%w", err)
	}

	// SAM header written as literal text; samtools sort receives SAM on stdin.
	var samHeader string = `@HD	VN:1.6	SO:unsorted	GO:query
@SQ	SN:*	LN:1
`
	var wg sync.WaitGroup
	errCh := make(chan error, 10)

	// ── samtools sort ──────────────────────────────────────────────────────────
	var tempFilePrefix string
	if filepath.Dir(prefix) == "." {
		tempFilePrefix = "." + RandString()
	} else {
		tempFilePrefix = filepath.Dir(prefix) + "." + RandString()
	}
	_threads := fmt.Sprintf("%d", threads-1)
	samSort := exec.Command("samtools", "sort", "-@", _threads, "--no-PG", "-T", tempFilePrefix, "-u", "-t", "BX", "-")
	sortIn, err := samSort.StdinPipe()
	if err != nil {
		return fmt.Errorf("samtools sort stdin pipe: %w", err)
	}
	sortOut, err := samSort.StdoutPipe()
	if err != nil {
		return fmt.Errorf("samtools sort stdout pipe: %w", err)
	}

	samSort.Stderr = os.Stderr

	if err := samSort.Start(); err != nil {
		return fmt.Errorf("starting samtools sort: %w", err)
	}

	// ---- Read input FASTQ -------------------------
	r1Reader, err := fastx.NewDefaultReader(r1Path)
	if err != nil {
		return fmt.Errorf("opening %s: %w", r1Path, err)
	}
	defer r1Reader.Close()

	r2Reader, err := fastx.NewDefaultReader(r2Path)
	if err != nil {
		return fmt.Errorf("opening %s: %w", r2Path, err)
	}
	defer r2Reader.Close()

	// ----- Write output FASTQ (Invalids) ------------------------
	r1InvWritePath := prefix + ".invalid.R1.fq.gz"
	r2InvWritePath := prefix + ".invalid.R2.fq.gz"

	r1InvWriter, err := xopen.Wopen(r1InvWritePath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", r1InvWritePath, err)
	}
	defer r1InvWriter.Close()

	r2InvWriter, err := xopen.Wopen(r2InvWritePath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", r2InvWritePath, err)
	}
	defer r2InvWriter.Close()

	// ---- process reads -------------------
	wg.Go(
		func() {
			defer sortIn.Close()
			w := bufio.NewWriterSize(sortIn, 2<<20)
			if _, err := w.WriteString(samHeader); err != nil {
				errCh <- fmt.Errorf("writing SAM header: %w", err)
				return
			}
			//errCh <- w.Flush()
			var rec1, rec2 *fastx.Record
			var err1, err2 error
			var isValid bool

			for {
				rec1, err1 = r1Reader.Read()
				if err1 != nil && err1 != io.EOF {
					errCh <- fmt.Errorf("reading R1: %w", err1)
					return
				}
				rec2, err2 = r2Reader.Read()
				if err2 != nil && err2 != io.EOF {
					errCh <- fmt.Errorf("reading R2: %w", err2)
					return
				}
				if err1 == io.EOF && err2 == io.EOF {
					break
				}
				if rec1 != nil {
					// STANDARDIZE and return if valid
					isValid = processBC(rec1)
					if !isValid {
						// send to invalid write
						if err := WriteFQ(r1InvWriter, rec1, _mark_forward); err != nil {
							errCh <- err
							return
						}
					} else {
						// CONVERT TO SAM AND WRITE TO SORT
						if err := Fq2Sam(rec1, fw, w); err != nil {
							errCh <- err
							return
						}
					}
				}
				if rec2 != nil {
					// STANDARDIZE and return if valid
					isValid = processBC(rec2)
					if !isValid {
						// send to invalid write
						if err := WriteFQ(r2InvWriter, rec2, _mark_reverse); err != nil {
							errCh <- err
							return
						}
					} else {
						// CONVERT TO SAM AND WRITE TO SORT
						if err := Fq2Sam(rec2, rv, w); err != nil {
							errCh <- err
							return
						}
					}
				}
			}
			if err := w.Flush(); err != nil {
				errCh <- fmt.Errorf("flushing SAM stream: %w", err)
				return
			}
		})

	//----- process sorted BAM records-----------------
	wg.Go(
		func() {
			samReader, err := bam.NewReader(sortOut, 1)
			if err != nil {
				errCh <- fmt.Errorf("creating BAM reader: %w", err)
				return
			}
			defer samReader.Close()
			// ----- Write output FASTQ (Valids) ------------------------
			r1WritePath := prefix + ".arachne.R1.fq.gz"
			r1Writer, err := xopen.Wopen(r1WritePath)
			if err != nil {
				errCh <- fmt.Errorf("opening %s: %w", r1WritePath, err)
			}
			defer r1Writer.Close()

			r2WritePath := prefix + "arachne.R2.fq.gz"
			r2Writer, err := xopen.Wopen(r2WritePath)
			if err != nil {
				errCh <- fmt.Errorf("opening %s: %w", r2WritePath, err)
			}
			defer r2Writer.Close()

			// ------ Process samtools sort output
			for {
				rec, err := samReader.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					errCh <- fmt.Errorf("reading sorted SAM: %w", err)
					return
				}
				if rec.Flags&sam.Read1 != 0 {
					Sam2FQ(r1Writer, rec, _mark_forward)
				} else {
					Sam2FQ(r2Writer, rec, _mark_reverse)
				}
			}
		})
	done := make(chan struct{})

	wg.Wait()
	close(errCh)
	close(done)

	if err := samSort.Wait(); err != nil {
		return fmt.Errorf("samtools sort failed: %w", err)
	}

	for err := range errCh {
		if err != nil {
			return err
		}
	}

	return nil
}
