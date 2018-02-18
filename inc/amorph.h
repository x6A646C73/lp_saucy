#ifndef SAUCY_AMORPH_H
#define SAUCY_AMORPH_H

#include "saucy.h"

struct amorph_graph 
{
    struct saucy_graph sg;
    int *colors;
    void *data;
    char **var_names;
    void (*consumer)( int, const int *, int, const int *, struct amorph_graph *g, char * );
    void (*free)( struct amorph_graph * );
    void (*stats)( struct amorph_graph *, FILE *f );
};

struct amorph_graph* read_lp( const char *filename );
struct amorph_graph* read_graph( const char *filename );

//TODO: put in functionality for copying graphs, remove edges, etc.
//      basically, put in what the heuristic needs

#endif
