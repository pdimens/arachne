package aligner

import (
	"fmt"
	"strings"

	"arachne/gobwa"

	"github.com/biogo/hts/sam"
)

// Method selects the algorithm used to resolve multi-mapping reads within
// a barcode.
type Method string

const (
	// MethodRFA is the original read-cloud optimizer (simulated annealing
	// over candidate molecules).
	MethodRFA Method = "rfa"
	// MethodEM is the EMA-style expectation-maximization over candidate
	// clouds.
	MethodEM Method = "em"
)

// ParseMethod converts a user-supplied string to a Method. The empty
// string selects the default (RFA).
func ParseMethod(s string) (Method, error) {
	switch Method(strings.ToLower(strings.TrimSpace(s))) {
	case "", MethodRFA:
		return MethodRFA, nil
	case MethodEM:
		return MethodEM, nil
	}
	return "", fmt.Errorf("unknown method %q (valid: %s, %s)", s, MethodRFA, MethodEM)
}

// BarcodeFunc processes every read pair of one barcode: it aligns them,
// resolves ambiguity, and writes the resulting records to out. Method
// specific configuration is captured by the closure.
type BarcodeFunc func(work *WorkUnit,
	out chan *sam.Record,
	ref *gobwa.GoBwaReference,
	settings *gobwa.GoBwaSettings,
	contigs map[string]*sam.Reference,
	debugtags *bool)

// barcodeFuncFor returns the per-barcode processor for the given method.
func barcodeFuncFor(m Method, improperPenalty float64) BarcodeFunc {
	switch m {
	case MethodEM:
		config := DefaultEMConfig(improperPenalty)
		return func(work *WorkUnit, out chan *sam.Record, ref *gobwa.GoBwaReference,
			settings *gobwa.GoBwaSettings, contigs map[string]*sam.Reference, debugtags *bool) {
			DoEMForOneBarcode(work, out, ref, settings, config, contigs, debugtags, work.reads)
		}
	default:
		config := &RFAConfig{improperPenalty}
		stats := &RFAStats{}
		return func(work *WorkUnit, out chan *sam.Record, ref *gobwa.GoBwaReference,
			settings *gobwa.GoBwaSettings, contigs map[string]*sam.Reference, debugtags *bool) {
			DoRFAForOneBarcode(work, out, ref, settings, config, stats, contigs, debugtags, work.reads)
		}
	}
}
