#ifndef WARMUP_THETA
#define WARMUP_THETA

#include "saucy.h"

struct lp_saucy
{ //{{{
    /* Graph data */
    uint64_t n;       /* Size of domain */
    uint64_t wcount;  /* number of edge colors */
    const int *adj;   /* Neighbors of k: edg[adj[k]]..edg[adj[k+1]] */
    const int *edg;   /* Actual neighbor data */
    const int *wght;  /* Actual edge colors */
    const int *dadj;  /* Fanin neighbor indices, for digraphs */
    const int *dedg;  /* Fanin neighbor data, for digraphs */
    const int *dwght; /* Fanin weight data, for digraphs */
    
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
    int ndiffs;      /* Current number of diffs */
    int *unsupp;     /* Inverted diff array */
    
    /* Statistics structure */
    struct saucy_stats *stats;
}; //}}} END lp_saucy

int lp_warmup_theta( char *filename, struct lp_saucy *s, int directed, 
                     const struct saucy_graph *g )
struct lp_saucy *lp_saucy_alloc( int n, int w, struct saucy_stats *stats );
void lp_saucy_free( struct lp_saucy *s );

#endif
