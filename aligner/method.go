package aligner

import (
	"arachne/gominibwa"

	"github.com/biogo/hts/sam"
)

// BarcodeFunc processes every read pair of one barcode: it aligns them,
// resolves ambiguity, and writes the resulting records to out. Algorithm
// specific configuration is captured by the closure.
type BarcodeFunc func(work *WorkUnit,
	out chan *sam.Record,
	mapper *gominibwa.Mapper,
	contigs map[string]*sam.Reference,
	debugtags *bool)

// barcodeFuncFor returns the per-barcode processor: the EMA-style
// expectation-maximization (the default), or RFA if rfa is true.
func barcodeFuncFor(rfa bool, improperPenalty float64, emConfig *EMConfig) BarcodeFunc {
	switch {
	case !rfa:
		config := emConfig
		if config == nil {
			config = DefaultEMConfig(improperPenalty)
		}
		return func(work *WorkUnit, out chan *sam.Record, mapper *gominibwa.Mapper,
			contigs map[string]*sam.Reference, debugtags *bool) {
			DoEMForOneBarcode(work, out, mapper, config, contigs, debugtags, work.reads)
		}
	default:
		config := &RFAConfig{improperPenalty}
		stats := &RFAStats{}
		return func(work *WorkUnit, out chan *sam.Record, mapper *gominibwa.Mapper,
			contigs map[string]*sam.Reference, debugtags *bool) {
			DoRFAForOneBarcode(work, out, mapper, config, stats, contigs, debugtags, work.reads)
		}
	}
}
