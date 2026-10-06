// Small helpers over minibwa's internal structures that its public header
// does not expose.

#ifndef GOMINIBWA_BRIDGE_H
#define GOMINIBWA_BRIDGE_H

#include <stdint.h>
#include "minibwa.h"

// number of contigs in the index
int64_t gmb_n_ctg(const mb_idx_t *idx);

// Fetch reference bases [st,en) of contig tid as 0..3 codes (A,C,G,T; 4 for
// ambiguous). Returns the number of bases written, or <0 on error.
int64_t gmb_getseq(const mb_idx_t *idx, int32_t tid, int64_t st, int64_t en, uint8_t *out);

// Allocate and initialise options; preset may be NULL. Returns NULL if the
// preset is unknown. Free with free().
mb_opt_t *gmb_opt_new(const char *preset);

// cgo cannot read C bitfields or flexible array members, so expose them.
#define GMB_REV    0x1
#define GMB_PROPER 0x2
#define GMB_SAMPRI 0x4
#define GMB_RESCUE 0x8
uint32_t gmb_hit_bits(const mb_hit_t *h);
// CIGAR of the aligned region (no clipping); NULL if the hit has none.
const uint32_t *gmb_hit_cigar(const mb_hit_t *h, int32_t *n_cigar);
// DP score of the base-level alignment; 0 if the hit has none.
int32_t gmb_hit_dp_score(const mb_hit_t *h);

#endif
