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

struct lp_saucy
{ //{{{
    /* Refinement: workspace */
    char *stuff;     /* Bit vector, but one char per bit */
    int *wstuff;     /* weight data */
    int *gamma;      /* Working permutation */

    /* Search: orbit partition */
    int *theta;      /* Running approximation of orbit partition */
    int *thsize;     /* Size of cells in theta, defined for mcrs */
    int *thnext;     /* Next rep in list (circular list) */
    int *thprev;     /* Previous rep in list */
    int *threp;      /* First rep for a given cell front */
    int *thfront;    /* The cell front associated with this rep */

    /* Polymorphic functions */
    int (*lp_is_automorphism)(struct saucy *);

     /* Statistics structure */
    struct saucy_stats *stats;
}; //}}} END saucy

struct lp_saucy *lp_saucy_alloc(int n, int w);

void lp_saucy_free(struct lp_saucy *s);

int lp_is_undirected_automorphism( struct lp_saucy *s );

int lp_is_directed_automorphism( struct lp_saucy *s );

#endif
