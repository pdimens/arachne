// bxsort.go — external sort of FASTQ records by BX:Z tag
// Uses xopen (transparent compression) and fastx (zero-copy FASTQ parsing).
//
// go get github.com/brentp/xopen
// go get github.com/shenwei356/bio/io/fastx
// go get github.com/shenwei356/bio/seq
package main

import (
	"bufio"
	"bytes"
	"container/heap"
	"encoding/binary"
	"errors"
	"flag"
	"io"
	"log"
	"os"
	"runtime"
	"sort"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

// ── Record ────────────────────────────────────────────────────────────────────
//
// All fields are owned copies — safe to hold after fastx recycles its buffer.
// We drop the original '+' line; writing bare '+' is standard FASTQ.

type Record struct {
	Name []byte // full header, without leading '@'
	Seq  []byte
	Qual []byte
	BX   []byte // nil when tag is absent
}

const prefix = "BX:Z:"

// extractBX scans the comment portion of a FASTQ name line for BX:Z:<value>.
// Comment fields are TAB-separated after the first space.
// All comparisons are on []byte — no string allocation.
func extractBX(name []byte) []byte {
	//i := bytes.IndexByte(name, ' ')
	_, comment, ok := bytes.Cut(name, []byte{' '})
	if !ok {
		return nil
	}
	//comment := name[i+1:]
	tag := []byte(prefix)
	for len(comment) > 0 {
		var field []byte
		if j := bytes.IndexByte(comment, '\t'); j < 0 {
			field, comment = comment, nil
		} else {
			field, comment = comment[:j], comment[j+1:]
		}
		if bytes.HasPrefix(field, tag) {
			// Must copy: field aliases fastx's recycled buffer.
			return append([]byte(nil), field[len(tag):]...)
		}
	}
	return nil
}

// bxLess orders records by BX tag, pushing untagged records to the end.
func bxLess(a, b Record) bool {
	switch {
	case len(a.BX) == 0:
		return false
	case len(b.BX) == 0:
		return true
	default:
		return bytes.Compare(a.BX, b.BX) < 0
	}
}

// ── Binary temp-file format ───────────────────────────────────────────────────
//
// Each record: 4 length-prefixed fields (uint32 LE): BX, Name, Seq, Qual.
// Simple and unambiguous; no dependency on newline framing.

func writeField(w io.Writer, b []byte) error {
	var hdr [4]byte
	binary.LittleEndian.PutUint32(hdr[:], uint32(len(b)))
	if _, err := w.Write(hdr[:]); err != nil {
		return err
	}
	_, err := w.Write(b)
	return err
}

func readField(r io.Reader) ([]byte, error) {
	var hdr [4]byte
	if _, err := io.ReadFull(r, hdr[:]); err != nil {
		return nil, err
	}
	b := make([]byte, binary.LittleEndian.Uint32(hdr[:]))
	_, err := io.ReadFull(r, b)
	return b, err
}

func writeBinaryRecord(w io.Writer, rec Record) error {
	for _, f := range [][]byte{rec.BX, rec.Name, rec.Seq, rec.Qual} {
		if err := writeField(w, f); err != nil {
			return err
		}
	}
	return nil
}

func readBinaryRecord(r io.Reader) (Record, error) {
	var fields [4][]byte
	for i := range fields {
		b, err := readField(r)
		if err != nil {
			return Record{}, err // io.EOF on first field → clean end of stream
		}
		fields[i] = b
	}
	return Record{BX: fields[0], Name: fields[1], Seq: fields[2], Qual: fields[3]}, nil
}

// ── FASTQ output ──────────────────────────────────────────────────────────────

var plus = []byte("+")

func writeRecord(w *bufio.Writer, rec Record) error {
	if err := w.WriteByte('@'); err != nil {
		return err
	}
	if _, err := w.Write(rec.Name); err != nil {
		return err
	}
	if err := w.WriteByte('\n'); err != nil {
		return err
	}
	if _, err := w.Write(rec.Seq); err != nil {
		return err
	}
	if _, err := w.Write([]byte("\n+\n")); err != nil {
		return err
	}
	if _, err := w.Write(rec.Qual); err != nil {
		return err
	}
	return w.WriteByte('\n')
}

// ── Sort worker ───────────────────────────────────────────────────────────────

type chunkJob struct {
	chunk []Record
}

type spillResult struct {
	file *os.File
	n    int
	err  error
}

// sortAndSpill sorts chunk in-place and writes it to a new temp file.
// Called exclusively by worker goroutines — no shared state.
func sortAndSpill(chunk []Record, dir string) (*os.File, int, error) {
	sort.Slice(chunk, func(i, j int) bool {
		return bxLess(chunk[i], chunk[j])
	})

	tmp, err := os.CreateTemp(dir, "bxsort-*.tmp")
	if err != nil {
		return nil, 0, err
	}

	bw := bufio.NewWriterSize(tmp, 1<<20) // 1 MiB write buffer
	for _, rec := range chunk {
		if err := writeBinaryRecord(bw, rec); err != nil {
			tmp.Close()
			os.Remove(tmp.Name())
			return nil, 0, err
		}
	}
	if err := bw.Flush(); err != nil {
		tmp.Close()
		os.Remove(tmp.Name())
		return nil, 0, err
	}
	if _, err := tmp.Seek(0, io.SeekStart); err != nil {
		tmp.Close()
		os.Remove(tmp.Name())
		return nil, 0, err
	}
	return tmp, len(chunk), nil
}

// startWorkers launches nWorkers goroutines that consume chunkJobs, sort each
// chunk, spill it to a temp file, and send a spillResult to results.
// The results channel is closed once all workers finish.
func startWorkers(nWorkers int, dir string, jobs <-chan chunkJob, results chan<- spillResult) {
	var wg sync.WaitGroup
	for range nWorkers {
		wg.Go(func() {
			for job := range jobs {
				f, n, err := sortAndSpill(job.chunk, dir)
				results <- spillResult{file: f, n: n, err: err}
			}
		})
	}
	go func() {
		wg.Wait()
		close(results)
	}()
}

// ── K-way merge heap ──────────────────────────────────────────────────────────

type heapItem struct {
	rec Record
	r   *bufio.Reader
	f   *os.File
}

type mergeHeap []*heapItem

func (h mergeHeap) Len() int           { return len(h) }
func (h mergeHeap) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
func (h mergeHeap) Less(i, j int) bool { return bxLess(h[i].rec, h[j].rec) }
func (h *mergeHeap) Push(x any)        { *h = append(*h, x.(*heapItem)) }
func (h *mergeHeap) Pop() any {
	old := *h
	n := len(old)
	x := old[n-1]
	old[n-1] = nil
	*h = old[:n-1]
	return x
}

// mergeSpills performs a K-way merge of sorted temp files into out.
// Memory usage: O(K) — one record per stream in the heap.
func mergeSpills(temps []*os.File, out *bufio.Writer) error {
	h := make(mergeHeap, 0, len(temps))
	for _, f := range temps {
		r := bufio.NewReaderSize(f, 1<<20)
		rec, err := readBinaryRecord(r)
		if err != nil {
			if errors.Is(err, io.EOF) {
				continue // empty spill — shouldn't happen, but be safe
			}
			return err
		}
		heap.Push(&h, &heapItem{rec: rec, r: r, f: f})
	}
	heap.Init(&h)

	for h.Len() > 0 {
		item := heap.Pop(&h).(*heapItem)
		if err := writeRecord(out, item.rec); err != nil {
			return err
		}
		next, err := readBinaryRecord(item.r)
		if err != nil {
			if !errors.Is(err, io.EOF) {
				return err
			}
			continue // stream exhausted
		}
		item.rec = next
		heap.Push(&h, item)
	}
	return nil
}

// ── Main ──────────────────────────────────────────────────────────────────────

func main() {
	chunkSize := flag.Int("n", 500_000, "records per sort chunk")
	tmpDir := flag.String("T", os.TempDir(), "directory for temporary files")
	outPath := flag.String("o", "-", "output FASTQ (- = stdout); .gz/.bgz for compressed output")
	nThreads := flag.Int("@", runtime.NumCPU(), "sort worker threads (like samtools sort -@)")
	flag.Parse()

	// ── Input ─────────────────────────────────────────────────────────────────
	inputPath := "-"
	if args := flag.Args(); len(args) > 0 {
		inputPath = args[0]
	}
	seq.ValidateSeq = false // skip per-base alphabet validation (~10% faster)
	reader, err := fastx.NewReader(seq.Unlimit, inputPath, "")
	if err != nil {
		log.Fatalf("fastx reader: %v", err)
	}
	defer reader.Close()

	// ── Output ────────────────────────────────────────────────────────────────
	// xopen.Wopen writes .gz/.bgz transparently based on extension.
	wfh, err := xopen.Wopen(*outPath)
	if err != nil {
		log.Fatalf("open output: %v", err)
	}
	defer wfh.Close()

	// ── Worker pool ───────────────────────────────────────────────────────────
	// Buffer depth = nThreads so the reader is never blocked waiting for a
	// free worker. Matches samtools' model: I/O and sorting overlap.
	jobs := make(chan chunkJob, *nThreads)
	results := make(chan spillResult, *nThreads)
	startWorkers(*nThreads, *tmpDir, jobs, results)

	// Collector: gathers spill results from workers asynchronously.
	var (
		temps      []*os.File
		collectMu  sync.Mutex
		collectErr error
		collectWg  sync.WaitGroup
	)
	collectWg.Go(func() {
		spillIdx := 0
		for res := range results {
			if res.err != nil {
				collectMu.Lock()
				if collectErr == nil {
					collectErr = res.err
				}
				collectMu.Unlock()
				continue
			}
			spillIdx++
			collectMu.Lock()
			temps = append(temps, res.file)
			collectMu.Unlock()
			log.Printf("  spill %d: %d records → %s", spillIdx, res.n, res.file.Name())
		}
	})

	// ── Pass 1: read → chunk → dispatch to workers ────────────────────────────
	log.Printf("pass 1: chunking at %d records, %d sort workers", *chunkSize, *nThreads)
	chunk := make([]Record, 0, *chunkSize)

	for {
		rec, err := reader.Read()
		if err != nil {
			if errors.Is(err, io.EOF) {
				break
			}
			log.Fatalf("read: %v", err)
		}

		// rec.Name, rec.Seq.Seq, rec.Qual are all slices into fastx's internal
		// recycled buffer — copy everything before the next Read() call.
		chunk = append(chunk, Record{
			Name: append([]byte(nil), rec.Name...),
			Seq:  append([]byte(nil), rec.Seq.Seq...),
			Qual: append([]byte(nil), rec.Seq.Qual...),
			BX:   extractBX(rec.Name), // already copies the BX value
		})

		if len(chunk) == *chunkSize {
			jobs <- chunkJob{chunk: chunk}
			// New backing array — the old one is now owned by the worker.
			chunk = make([]Record, 0, *chunkSize)
		}
	}
	if len(chunk) > 0 {
		jobs <- chunkJob{chunk: chunk}
	}
	close(jobs)      // signal workers: no more chunks
	collectWg.Wait() // wait for all spill results to be collected

	if collectErr != nil {
		log.Fatalf("spill: %v", collectErr)
	}

	// Ensure temp files are removed even on merge failure.
	defer func() {
		for _, f := range temps {
			f.Close()
			os.Remove(f.Name())
		}
	}()

	// ── Pass 2: K-way merge ───────────────────────────────────────────────────
	log.Printf("pass 2: merging %d sorted chunks", len(temps))
	out := bufio.NewWriterSize(wfh, 1<<20)
	if err := mergeSpills(temps, out); err != nil {
		log.Fatalf("merge: %v", err)
	}
	if err := out.Flush(); err != nil {
		log.Fatalf("flush: %v", err)
	}
	log.Println("done.")
}
