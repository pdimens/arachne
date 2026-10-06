//go:build minibwa

package gominibwa

// #cgo CFLAGS: -I${SRCDIR}/minibwa
// #cgo LDFLAGS: -lminibwa -lz -lm -lpthread
// #include <stdlib.h>
// #include "bridge.h"
import "C"
import (
	"bytes"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"unsafe"
)

// Index is a loaded minibwa index. It is read-only once loaded and safe to
// share between goroutines; each goroutine needs its own Mapper.
type Index struct {
	idx *C.mb_idx_t
}

// LoadIndex loads the index with the given prefix, which is the path of the
// reference FASTA unless a different prefix was used when building.
func LoadIndex(prefix string) (*Index, error) {
	cprefix := C.CString(prefix)
	defer C.free(unsafe.Pointer(cprefix))
	idx := C.mb_idx_load(cprefix, 0)
	if idx == nil {
		return nil, fmt.Errorf("gominibwa: could not load index %q (expected %s.l2b and %s.mbw)", prefix, prefix, prefix)
	}
	x := &Index{idx: idx}
	runtime.SetFinalizer(x, func(x *Index) { x.Close() })
	return x, nil
}

// Close releases the index. Mappers created from it must not be used
// afterwards.
func (x *Index) Close() {
	if x.idx != nil {
		C.mb_idx_destroy(x.idx)
		x.idx = nil
	}
}

// NumContigs returns the number of reference contigs.
func (x *Index) NumContigs() int {
	return int(C.gmb_n_ctg(x.idx))
}

// ContigName returns the name of contig tid, or "" if tid is out of range.
func (x *Index) ContigName(tid int) string {
	p := C.mb_idx_ctg_name(x.idx, C.int32_t(tid))
	if p == nil {
		return ""
	}
	return C.GoString(p)
}

// ContigLength returns the length of contig tid, or -1 if out of range.
func (x *Index) ContigLength(tid int) int64 {
	return int64(C.mb_idx_ctg_len(x.idx, C.int32_t(tid)))
}

// Contigs returns all contig names and lengths in index order, like
// gobwa.GoBwaReference.GetReferenceContigsInfo.
func (x *Index) Contigs() ([]string, []int64) {
	n := x.NumContigs()
	names := make([]string, n)
	lengths := make([]int64, n)
	for i := range n {
		names[i] = x.ContigName(i)
		lengths[i] = x.ContigLength(i)
	}
	return names, lengths
}

var baseLetters = [5]byte{'A', 'C', 'G', 'T', 'N'}

// GetSeq returns reference bases [start, end) of contig tid, optionally
// reverse complemented, like gobwa.GoBwaReference.GetSeq.
func (x *Index) GetSeq(tid int, start, end int64, reversed bool) ([]byte, error) {
	if end < start {
		return nil, fmt.Errorf("gominibwa: invalid interval [%d,%d)", start, end)
	}
	n := end - start
	if n == 0 {
		return []byte{}, nil
	}
	buf := (*C.uint8_t)(C.malloc(C.size_t(n)))
	defer C.free(unsafe.Pointer(buf))
	got := int64(C.gmb_getseq(x.idx, C.int32_t(tid), C.int64_t(start), C.int64_t(end), buf))
	if got < 0 {
		return nil, fmt.Errorf("gominibwa: could not fetch contig %d [%d,%d)", tid, start, end)
	}
	codes := unsafe.Slice((*byte)(unsafe.Pointer(buf)), int(got))
	out := make([]byte, got)
	for i, c := range codes {
		if c > 4 {
			c = 4
		}
		if reversed {
			out[int(got)-1-i] = complement(baseLetters[c])
		} else {
			out[i] = baseLetters[c]
		}
	}
	return out, nil
}

func complement(b byte) byte {
	switch b {
	case 'A':
		return 'T'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'T':
		return 'A'
	}
	return b
}

// minibwaPath finds the minibwa binary next to the running executable, then
// on PATH.
func minibwaPath() (string, error) {
	if exe, err := os.Executable(); err == nil {
		p := filepath.Join(filepath.Dir(exe), "minibwa")
		if st, err := os.Stat(p); err == nil && !st.IsDir() {
			return p, nil
		}
	}
	return exec.LookPath("minibwa")
}

// Build indexes a reference FASTA by running `minibwa index`, producing
// <fasta>.l2b and <fasta>.mbw. It is the counterpart of gobwa.GoBwaIndex.
func Build(fastaPath string, threads int) error {
	bin, err := minibwaPath()
	if err != nil {
		return fmt.Errorf("gominibwa: minibwa binary not found: %w", err)
	}
	args := []string{"index"}
	if threads > 0 {
		args = append(args, fmt.Sprintf("-t%d", threads))
	}
	cmd := exec.Command(bin, append(args, fastaPath)...)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		return fmt.Errorf("minibwa index failed: %w\n%s", err, stderr.String())
	}
	return nil
}

// Flag bits for Options.SetFlag, mirroring MB_F_* in minibwa.h.
const (
	FlagPairedEnd uint64 = 0x8
	FlagEQX       uint64 = 0x20
	FlagNoAlign   uint64 = 0x80
	FlagNoPairing uint64 = 0x10000
)

// Options wraps minibwa's mb_opt_t. Create with NewOptions; the zero value
// is not usable. Options are read-only while mapping and may be shared.
type Options struct {
	opt *C.mb_opt_t
}

// NewOptions creates options from a preset: "sr" (short reads), "lr" (long
// reads), or "adap" (adaptive, minibwa's default). An empty preset gives
// "adap".
func NewOptions(preset string) (*Options, error) {
	var cpreset *C.char
	if preset != "" {
		cpreset = C.CString(preset)
		defer C.free(unsafe.Pointer(cpreset))
	}
	opt := C.gmb_opt_new(cpreset)
	if opt == nil {
		return nil, fmt.Errorf("gominibwa: unknown preset %q", preset)
	}
	o := &Options{opt: opt}
	runtime.SetFinalizer(o, func(o *Options) { C.free(unsafe.Pointer(o.opt)) })
	return o, nil
}

// SetFlag turns the given MB_F_* flag bits on or off.
func (o *Options) SetFlag(flag uint64, on bool) {
	if on {
		o.opt.flag |= C.uint64_t(flag)
	} else {
		o.opt.flag &^= C.uint64_t(flag)
	}
}

// Flags returns the current MB_F_* flag bits.
func (o *Options) Flags() uint64 { return uint64(o.opt.flag) }

// CigarOp is one CIGAR operation.
type CigarOp struct {
	Op  byte // one of "MIDNSHP=XB"
	Len int
}

const cigarChars = "MIDNSHP=XB"

// Hit is one alignment of a read.
type Hit struct {
	Tid      int    // contig index
	Contig   string // contig name
	Start    int64  // 0-based start on the forward strand
	End      int64  // exclusive end
	Reversed bool   // the read aligns to the reverse strand

	QueryStart int // aligned query interval [QueryStart, QueryEnd)
	QueryEnd   int

	Score        int // chaining score
	DPScore      int // base-level alignment score
	Matches      int // matching bases in the aligned block
	BlockLen     int // alignment block length
	EditDistance int // BlockLen - Matches, as minimap2's NM
	MapQ         int

	// Cigar covers the aligned region only: no soft or hard clipping. Use
	// QueryStart and QueryEnd, and the read length, to recover clipping.
	Cigar []CigarOp

	Secondary     bool // not the best hit of its group
	Supplementary bool // non-primary segment of a split alignment
	ProperPair    bool // paired and consistent with the insert size distribution
	Rescued       bool // found by mate rescue
}

// Primary reports whether this is the hit SAM would call primary.
func (h *Hit) Primary() bool { return !h.Secondary && !h.Supplementary }

// Mapper aligns reads against an Index. A Mapper owns thread-local buffers
// and must not be used from more than one goroutine at a time.
type Mapper struct {
	idx  *Index
	opt  *Options
	tbuf *C.mb_tbuf_t
}

// NewMapper creates a Mapper. Create one per worker goroutine.
func (x *Index) NewMapper(opt *Options) *Mapper {
	m := &Mapper{idx: x, opt: opt, tbuf: C.mb_tbuf_init(0)}
	runtime.SetFinalizer(m, func(m *Mapper) { m.Close() })
	return m
}

// Close releases the Mapper's buffers.
func (m *Mapper) Close() {
	if m.tbuf != nil {
		C.mb_tbuf_destroy(m.tbuf)
		m.tbuf = nil
	}
}

// Map aligns a single read.
func (m *Mapper) Map(name string, seq []byte) []Hit {
	if len(seq) == 0 {
		return nil
	}
	cname := C.CString(name)
	defer C.free(unsafe.Pointer(cname))
	cseq := C.CBytes(seq)
	defer C.free(cseq)

	var nHit C.int32_t
	hits := C.mb_map(m.opt.opt, m.idx.idx, C.int32_t(len(seq)), (*C.char)(cseq), 0, &nHit, m.tbuf, cname)
	return m.collect(hits, int(nHit))
}

// MapPair aligns the two reads of a pair together, so that mate rescue and
// pairing information (ProperPair) are applied. Options must have
// FlagPairedEnd set (the default for the "sr" and "adap" presets).
func (m *Mapper) MapPair(name string, read1, read2 []byte) ([]Hit, []Hit, error) {
	hits, err := m.MapPairs([]string{name}, [][]byte{read1}, [][]byte{read2})
	if err != nil {
		return nil, nil, err
	}
	return hits[0][0], hits[0][1], nil
}

// MapPairs aligns many read pairs in one batch, which is faster than
// aligning pairs one at a time. names, reads1 and reads2 must have equal
// length. The result is indexed [pair][0 for read 1, 1 for read 2].
func (m *Mapper) MapPairs(names []string, reads1, reads2 [][]byte) ([][2][]Hit, error) {
	n := len(names)
	if len(reads1) != n || len(reads2) != n {
		return nil, errors.New("gominibwa: names, reads1 and reads2 must have the same length")
	}
	if n == 0 {
		return nil, nil
	}
	for i := range n {
		if len(reads1[i]) == 0 || len(reads2[i]) == 0 {
			return nil, fmt.Errorf("gominibwa: pair %d (%s) has an empty read", i, names[i])
		}
	}
	// Pairs are consecutive entries (2k, 2k+1) of the batch.
	total := 2 * n
	ptrSize := C.size_t(unsafe.Sizeof(uintptr(0)))
	cseqs := (**C.char)(C.malloc(C.size_t(total) * ptrSize))
	cnames := (**C.char)(C.malloc(C.size_t(total) * ptrSize))
	qlens := (*C.int32_t)(C.malloc(C.size_t(total) * C.size_t(unsafe.Sizeof(C.int32_t(0)))))
	nHits := (*C.int32_t)(C.calloc(C.size_t(total), C.size_t(unsafe.Sizeof(C.int32_t(0)))))
	seqSlice := unsafe.Slice(cseqs, total)
	nameSlice := unsafe.Slice(cnames, total)
	qlenSlice := unsafe.Slice(qlens, total)
	defer func() {
		for i := range total {
			C.free(unsafe.Pointer(seqSlice[i]))
			C.free(unsafe.Pointer(nameSlice[i]))
		}
		C.free(unsafe.Pointer(cseqs))
		C.free(unsafe.Pointer(cnames))
		C.free(unsafe.Pointer(qlens))
		C.free(unsafe.Pointer(nHits))
	}()
	for i := range n {
		for j, read := range [2][]byte{reads1[i], reads2[i]} {
			k := 2*i + j
			seqSlice[k] = (*C.char)(C.CBytes(read))
			nameSlice[k] = C.CString(names[i])
			qlenSlice[k] = C.int32_t(len(read))
		}
	}

	hitArr := C.mb_map_batch(m.opt.opt, m.idx.idx, C.int32_t(total), qlens,
		(**C.char)(unsafe.Pointer(cseqs)), nHits, m.tbuf, (**C.char)(unsafe.Pointer(cnames)))
	if hitArr == nil {
		return nil, errors.New("gominibwa: mb_map_batch failed")
	}
	defer C.free(unsafe.Pointer(hitArr))
	perSeq := unsafe.Slice(hitArr, total)
	counts := unsafe.Slice(nHits, total)

	out := make([][2][]Hit, n)
	for i := range n {
		out[i][0] = m.collect(perSeq[2*i], int(counts[2*i]))
		out[i][1] = m.collect(perSeq[2*i+1], int(counts[2*i+1]))
	}
	return out, nil
}

// collect converts a C hit array into Go values and frees the C memory.
func (m *Mapper) collect(hits *C.mb_hit_t, n int) []Hit {
	if hits == nil {
		return nil
	}
	defer C.free(unsafe.Pointer(hits))
	if n == 0 {
		return nil
	}
	chits := unsafe.Slice(hits, n)
	out := make([]Hit, n)
	for i := range chits {
		h := &chits[i]
		bits := uint32(C.gmb_hit_bits(h))
		hit := Hit{
			Tid:           int(h.tid),
			Contig:        m.idx.ContigName(int(h.tid)),
			Start:         int64(h.ts),
			End:           int64(h.te),
			Reversed:      bits&C.GMB_REV != 0,
			QueryStart:    int(h.qs),
			QueryEnd:      int(h.qe),
			Score:         int(h.score),
			DPScore:       int(C.gmb_hit_dp_score(h)),
			Matches:       int(h.mlen),
			BlockLen:      int(h.blen),
			EditDistance:  int(h.blen) - int(h.mlen),
			MapQ:          int(h.mapq),
			Secondary:     h.parent != h.id,
			Supplementary: h.parent == h.id && bits&C.GMB_SAMPRI == 0,
			ProperPair:    bits&C.GMB_PROPER != 0,
			Rescued:       bits&C.GMB_RESCUE != 0,
		}
		var nCigar C.int32_t
		if cig := C.gmb_hit_cigar(h, &nCigar); cig != nil {
			ops := unsafe.Slice(cig, int(nCigar))
			hit.Cigar = make([]CigarOp, len(ops))
			for k, c := range ops {
				op := int(c & 0xf)
				ch := byte('?')
				if op < len(cigarChars) {
					ch = cigarChars[op]
				}
				hit.Cigar[k] = CigarOp{Op: ch, Len: int(c >> 4)}
			}
			C.free(unsafe.Pointer(h.p))
		}
		out[i] = hit
	}
	return out
}
