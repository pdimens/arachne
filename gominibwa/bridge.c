//go:build minibwa

#include <stdlib.h>
#include "mbpriv.h"
#include "l2bit.h"
#include "bridge.h"

int64_t gmb_n_ctg(const mb_idx_t *idx)
{
	return (int64_t)idx->l2b->n_ctg;
}

int64_t gmb_getseq(const mb_idx_t *idx, int32_t tid, int64_t st, int64_t en, uint8_t *out)
{
	return l2b_getseq(idx->l2b, tid, st, en, out);
}

mb_opt_t *gmb_opt_new(const char *preset)
{
	mb_opt_t *opt = (mb_opt_t*)malloc(sizeof(mb_opt_t));
	if (opt == NULL) return NULL;
	mb_opt_init(opt);
	if (preset != NULL && mb_opt_preset(opt, preset) < 0) {
		free(opt);
		return NULL;
	}
	return opt;
}

uint32_t gmb_hit_bits(const mb_hit_t *h)
{
	return (h->rev? GMB_REV : 0) | (h->proper_pair? GMB_PROPER : 0) |
	       (h->sam_pri? GMB_SAMPRI : 0) | (h->rescued? GMB_RESCUE : 0);
}

const uint32_t *gmb_hit_cigar(const mb_hit_t *h, int32_t *n_cigar)
{
	if (h->p == NULL) { *n_cigar = 0; return NULL; }
	*n_cigar = h->p->n_cigar;
	return h->p->cigar;
}

int32_t gmb_hit_dp_score(const mb_hit_t *h)
{
	return h->p? h->p->dp_score : 0;
}
