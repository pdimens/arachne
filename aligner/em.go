package aligner

import (
	"fmt"
	"math"
	"os"

	"arachne/fastqreader"
	"arachne/gobwa"

	"github.com/biogo/hts/sam"
)

// The EM resolution method is modelled on EMA (Shajii et al., Cell Systems
// 2019; https://github.com/arshajii/ema). Within one barcode, every
// candidate alignment of every read is assigned to a "cloud" (a run of
// candidate positions separated by less than --infer-distance). Each
// candidate carries a responsibility (gamma): the posterior probability
// that the read truly originates there. The EM alternates between
//
//	E: gamma_i ∝ P(alignment_i) * w(cloud_i) * P(best mate placement in the same cloud)
//	M: w(cloud) = Σ gamma_i of the candidates in that cloud
//
// for a fixed number of iterations. Each read is then assigned its
// highest-gamma candidate and gamma is converted to a MAPQ.
//
// Differences from EMA, kept deliberately small for this first
// implementation:
//   - cloud weights are always normalised per read across that read's
//     candidate clouds (EMA's "many_clouds" branch). EMA's default branch
//     instead links the sub-clouds produced by splitting a cloud in which
//     one read collides with itself (split.c); that disjoint-set machinery
//     is not ported.
//   - alignment likelihoods are arachne's log_alignment_probability
//     (mismatch/indel/soft-clip penalties) rather than EMA's own model.
//   - MAPQ is capped using arachne's pseudo-count alignment, not BWA's MAPQ.

// EMConfig holds parameters of the EM method. Constant after creation.
type EMConfig struct {
	// log10 penalty for an improperly paired / unpaired placement. <= 0.
	ImproperPenalty float64
	// Number of EM iterations (EMA: 5).
	Iterations int
	// Minimum number of read pairs in a barcode to run the EM. Smaller
	// barcodes keep the initial best-pair placement. Defaults to 3, the same
	// threshold RFA uses; EMA itself uses 30, but on small simulated barcodes
	// that left most multi-mapping reads unresolved.
	MinPairs int
	// Stop iterating early once no gamma changes by more than this.
	Tolerance float64
}

// DefaultEMConfig returns EMA's defaults with the given improper pair penalty.
func DefaultEMConfig(improperPenalty float64) *EMConfig {
	return &EMConfig{
		ImproperPenalty: improperPenalty,
		Iterations:      5,
		MinPairs:        3,
		Tolerance:       1e-6,
	}
}

// Minimum cloud weight, so that log() of an empty cloud stays finite.
const emWeightFloor = 1e-30

// A cluster of nearby candidate alignments on one contig.
type emCloud struct {
	id     int
	contig string
	start  int64
	stop   int64
	expCov float64 // expected number of reads, Σ gamma of members
	weight float64
}

// One candidate placement of one read.
type emCand struct {
	aln      *Alignment
	cloud    int     // index into the clouds slice, -1 if unplaced
	logScore float64 // natural log alignment likelihood
	gamma    float64 // responsibility
}

// Determine if a barcode is worth running the EM on.
func worthRunningEM(barcode_fragments []fastqreader.FastQRecord, uniqueBarcode bool, config *EMConfig) bool {
	return uniqueBarcode && len(barcode_fragments) > 0 && len(barcode_fragments) >= config.MinPairs
}

func DoEMForOneBarcode(work *WorkUnit,
	out chan *sam.Record,
	ref *gobwa.GoBwaReference,
	settings *gobwa.GoBwaSettings,
	config *EMConfig,
	contigs map[string]*sam.Reference,
	debugtags *bool,
	reads []fastqreader.FastQRecord) {
	barcode_reads := work.reads
	arena := gobwa.NewArena()
	worthEM := worthRunningEM(barcode_reads, work.unique_barcode, config)
	barcode_chains, _ := GetChains(ref, settings, barcode_reads, arena, 25)
	alignments, stashed_alignments := GetAlignments(ref, settings, barcode_chains, 17, arena)

	// also sets the initial best-pair placement (active) for every read
	positions := tagBestAlignments(alignments)

	if len(barcode_reads) > 2 && *verbose {
		fmt.Fprintf(os.Stderr, "working on barcode %s  num reads: %d  doing EM: %v  unique_barcode %v\n",
			string(barcode_reads[0].Barcode),
			len(barcode_reads),
			worthEM,
			work.unique_barcode)
	}

	if worthEM {
		runEM(alignments, positions, config)
	} else {
		estimateMapQualities(alignments, nil, config.ImproperPenalty)
	}
	markDuplicates(alignments)
	CheckSplitReads(stashed_alignments, centromeres)
	flushToChannel(alignments, out, contigs, debugtags)
	ReturnBuffer(reads)
	arena.Free()
}

// runEM resolves the candidate alignments of one barcode in place: on
// return exactly one alignment per read is active, mates are linked, and
// MAPQ / mapq_data are filled in. positions are the candidate alignments
// sorted by position per contig, as returned by tagBestAlignments.
func runEM(alignments [][]*Alignment, positions [][]*Alignment, config *EMConfig) {
	clouds, cloudOf := buildEMClouds(positions, *inferDistance)
	cands := make([][]emCand, len(alignments))
	for read_id, alns := range alignments {
		cands[read_id] = make([]emCand, len(alns))
		for i, aln := range alns {
			cloud, ok := cloudOf[aln]
			if !ok {
				cloud = -1
			}
			cands[read_id][i] = emCand{
				aln:      aln,
				cloud:    cloud,
				logScore: aln.log_alignment_probability * math.Ln10,
			}
		}
	}

	emInitialize(cands, clouds)
	for range config.Iterations {
		if emIterate(cands, clouds, config) < config.Tolerance {
			break
		}
	}
	emSelect(cands, clouds)
}

// Group candidate alignments into clouds. positions must be sorted by
// position within each contig. Alignments without a position never reach
// here (tagBestAlignments excludes them).
func buildEMClouds(positions [][]*Alignment, dist int64) ([]*emCloud, map[*Alignment]int) {
	clouds := []*emCloud{}
	cloudOf := map[*Alignment]int{}
	for _, list := range positions {
		var current *emCloud
		for i, aln := range list {
			if i == 0 || aln.pos-list[i-1].pos > dist {
				current = &emCloud{id: len(clouds), contig: aln.contig, start: aln.pos}
				clouds = append(clouds, current)
			}
			current.stop = aln.pos
			cloudOf[aln] = current.id
		}
	}
	return clouds, cloudOf
}

// Initial responsibilities from alignment likelihoods alone, and the
// cloud weights they imply.
func emInitialize(cands [][]emCand, clouds []*emCloud) {
	probs := []float64{}
	for r := range cands {
		probs = probs[:0]
		for _, c := range cands[r] {
			probs = append(probs, c.logScore)
		}
		normalizeLogProbs(probs)
		for i := range cands[r] {
			cands[r][i].gamma = probs[i]
		}
	}
	emUpdateWeights(cands, clouds)
}

// M-step: cloud weight is the expected number of reads in it.
func emUpdateWeights(cands [][]emCand, clouds []*emCloud) {
	for _, c := range clouds {
		c.expCov = 0
	}
	for r := range cands {
		for _, c := range cands[r] {
			if c.cloud >= 0 {
				clouds[c.cloud].expCov += c.gamma
			}
		}
	}
	for _, c := range clouds {
		c.weight = c.expCov
	}
}

// One E-step followed by the M-step. Returns the largest change in any
// gamma. Responsibilities are updated synchronously (all reads see the
// previous iteration's gammas), so the result is independent of read order.
func emIterate(cands [][]emCand, clouds []*emCloud, config *EMConfig) float64 {
	improperNat := config.ImproperPenalty * math.Ln10
	next := make([][]float64, len(cands))
	for r := range cands {
		n := len(cands[r])
		next[r] = make([]float64, n)

		// each read normalises the weights of the clouds it could be in
		total := 0.0
		for _, c := range cands[r] {
			total += cloudWeight(clouds, c.cloud)
		}
		var mates []emCand
		if n > 0 {
			mates = cands[cands[r][0].aln.mate_id]
		}
		for i, c := range cands[r] {
			prior := 1.0
			if total > 0 {
				prior = cloudWeight(clouds, c.cloud) / total
			}
			next[r][i] = c.logScore + math.Log(math.Max(prior, emWeightFloor)) + bestMateScore(c, mates, improperNat)
		}
		normalizeLogProbs(next[r])
	}

	maxChange := 0.0
	for r := range cands {
		for i := range cands[r] {
			maxChange = math.Max(maxChange, math.Abs(cands[r][i].gamma-next[r][i]))
			cands[r][i].gamma = next[r][i]
		}
	}
	emUpdateWeights(cands, clouds)
	return maxChange
}

func cloudWeight(clouds []*emCloud, id int) float64 {
	if id < 0 {
		return emWeightFloor
	}
	return math.Max(clouds[id].weight, emWeightFloor)
}

// Best log score of the mate being placed consistently with candidate c:
// in the same cloud, on the opposite strand, weighted by the mate's own
// responsibility. An unpaired placement scores the improper pair penalty.
func bestMateScore(c emCand, mates []emCand, improperNat float64) float64 {
	best := improperNat
	if c.cloud < 0 {
		return best
	}
	for _, m := range mates {
		if m.cloud != c.cloud || m.aln.reversed == c.aln.reversed || m.gamma == 0 {
			continue
		}
		penalty := improperNat
		if isPair(c.aln, m.aln) {
			penalty = 0
		}
		if s := penalty + math.Log(m.gamma); s > best {
			best = s
		}
	}
	return best
}

// Pick the highest-responsibility candidate for each read, link mates,
// and convert responsibilities to MAPQ.
func emSelect(cands [][]emCand, clouds []*emCloud) {
	chosen := make([]int, len(cands))
	for r := range cands {
		best := 0
		for i, c := range cands[r] {
			if c.gamma > cands[r][best].gamma {
				best = i
			}
			c.aln.active = false
			c.aln.is_proper = false
		}
		chosen[r] = best
	}
	readsInCloud := map[int]int{}
	for r := range cands {
		if len(cands[r]) == 0 {
			continue
		}
		c := cands[r][chosen[r]]
		c.aln.active = true
		c.aln.bwa_pick = false
		if c.cloud >= 0 {
			readsInCloud[c.cloud]++
		}
	}
	for r := range cands {
		if len(cands[r]) == 0 {
			continue
		}
		a := cands[r][chosen[r]].aln
		if a.mate_id >= 0 && a.mate_id < len(cands) && len(cands[a.mate_id]) > 0 {
			refreshActiveMatePairing(a, cands[a.mate_id][chosen[a.mate_id]].aln)
		}
	}

	for r := range cands {
		if len(cands[r]) == 0 {
			continue
		}
		emAssignMapq(cands[r], chosen[r], clouds, readsInCloud)
	}
}

// Fill in MAPQ and the mapq_data used for SAM tags for the chosen
// candidate of one read.
func emAssignMapq(read []emCand, chosen int, clouds []*emCloud, readsInCloud map[int]int) {
	best := read[chosen]
	aln := best.aln

	// Confidence that the read belongs to any of its candidates at all, as
	// opposed to somewhere not found: weigh the candidates against arachne's
	// pseudo-count alignment (log10 domain, same as RFA's estimate).
	pseudo := psuedoCountAlignmentScore(aln, 0.0)
	shift := pseudo
	for _, c := range read {
		shift = math.Max(shift, c.aln.log_alignment_probability)
	}
	pseudoMass := math.Exp((pseudo - shift) * math.Ln10)
	candMass := 0.0
	for _, c := range read {
		candMass += math.Exp((c.aln.log_alignment_probability - shift) * math.Ln10)
	}
	confidence := 1.0 - pseudoMass/(pseudoMass+candMass)

	p := best.gamma * confidence
	mapq := 60.0
	if p < 0.999999 {
		mapq = math.Min(60.0, -10.0*math.Log10(1.0-p))
	}
	if aln.pos == -1 {
		mapq = 0
	} else if inCentromere(centromeres, aln.contig, aln.pos) {
		mapq = 0
	}
	aln.mapq = int(mapq)

	// runner-up candidate by responsibility
	second := -1
	for i, c := range read {
		if i != chosen && (second < 0 || c.gamma > read[second].gamma) {
			second = i
		}
	}
	data := aln.mapq_data
	data.copies = len(read)
	data.unique_molecules_active = 1
	data.score = scoreAlignment(aln, aln.mate_alignment, 0.0)
	data.second_best_score = pseudo
	data.second_best = nil
	data.second_best_proper_pair = false
	if second >= 0 {
		alt := read[second].aln
		data.second_best = alt
		data.second_best_score = scoreAlignment(alt, aln.mate_alignment, 0.0)
		data.second_best_proper_pair = aln.mate_alignment != nil && isPair(alt, aln.mate_alignment)
	}
	if best.cloud >= 0 {
		aln.molecule_id = best.cloud
		aln.active_molecule = true
		data.reads_in_molecule = readsInCloud[best.cloud]
	}
	aln.molecule_confidence = best.gamma
}

// normalizeLogProbs converts natural-log probabilities in place into a
// normalised distribution (softmax). Entries more than ~115 nats below the
// maximum are clamped to exactly zero, as in EMA.
func normalizeLogProbs(p []float64) {
	n := len(p)
	if n == 0 {
		return
	}
	if n == 1 {
		p[0] = 1.0
		return
	}
	thresh := math.Log(1e-50) - math.Log(float64(n))
	pmax := p[0]
	for _, v := range p[1:] {
		pmax = math.Max(pmax, v)
	}
	total := 0.0
	for i := range p {
		p[i] -= pmax
		if p[i] < thresh {
			p[i] = 0
		} else {
			p[i] = math.Exp(p[i])
		}
		total += p[i]
	}
	for i := range p {
		p[i] /= total
	}
}
