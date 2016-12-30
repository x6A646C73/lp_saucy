#ifndef LP_AMORPH_H
#define LP_AMORPH_H

#include "saucy.h"

struct lp_amorph_graph {
	struct saucy_graph sg;
	int *colors;
	void *data;
	void (*consumer)(int, const int *, int, const int *,
		struct lp_amorph_graph *g, char *);
	void (*free)(struct lp_amorph_graph *);
	void (*stats)(struct lp_amorph_graph *, FILE *f);
};

struct lp_amorph_graph* lp_amorph_read_build( char *filename );
int lp_warmup_theta( char *filename, struct lp_saucy *s );

#endif
