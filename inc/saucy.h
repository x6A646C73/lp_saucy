#ifndef SAUCY_H
#define SAUCY_H

#define SAUCY_VERSION "2.0"

typedef int saucy_consumer(int, const int *, int, int *, void *);

struct coloring
{ //{{{
    int *lab;        /* Labelling of objects */
    int *unlab;      /* Inverse of lab */
    int *cfront;     /* Pointer to front of cells */
    int *clen;       /* Length of cells (defined for cfront's) */
}; //}}} END coloring

struct saucy
{ //{{{
    /* Graph data */
    uint64_t n;       /* Size of domain */
    uint64_t wcount;  /* number of edge colors */
    const int *adj;   /* Neighbors of k: edg[adj[k]]..edg[adj[k+1]] */
    //NOTE: this is the range from adj[k] to adj[k+1], that is edges from adj[k]
    //      to adj[k+1]-1 adj[k] holds the index for the start of neighbors in
    //      *edg
    const int *edg;   /* Actual neighbor data */
    const int *wght;  /* Actual edge colors */
    const int *dadj;  /* Fanin neighbor indices, for digraphs */
    const int *dedg;  /* Fanin neighbor data, for digraphs */
    const int *dwght; /* Fanin weight data, for digraphs */
    void *arg;        /* Opaque client data */

    /* Coloring data */
    struct coloring left, right;
    int *nextnon;    /* Forward next-nonsingleton pointers */ 
    int *prevnon;    /* Backward next-nonsingleton pointers */

    /* Refinement: inducers */
    char *indmark;   /* Induce marks */
    int *ninduce;    /* Nonsingletons that might induce refinement */
    int *sinduce;    /* Singletons that might induce refinement */
    int nninduce;    /* Size of ninduce stack */
    int nsinduce;    /* Size of sinduce stack */

    /* Refinement: marked cells */
    int *clist;      /* List of cells marked for refining */
    int csize;       /* Number of cells in clist */

    /* Refinement: workspace */
    char *stuff;     /* Bit vector, but one char per bit */
    int *wstuff;     /* weight data */
    int *ccount;     /* Number of connections to refining cell */
    int *dccount;     /* Number of connections to refining cell */
    int *bucket;     /* Workspace */
    int *count;      /* Num vertices with same adj count to ref cell */
    int *junk;       /* More workspace */
    int *diffL;   /* Weight workspace */
    int *gamma;      /* Working permutation */
    int *conncnts;   /* Connection counts for cell fronts */

    /* Search data */
    int lev;         /* Current search tree level */
    int anc;         /* Level of greatest common ancestor with zeta */
    int *anctar;     /* Copy of target cell at anc */
    int kanctar;     /* Location within anctar to iterate from */
    int *start;      /* Location of target at each level */
    int indmin;      /* Used for group size computation */
    int match;       /* Have we not diverged from previous left? */

    /* Search: orbit partition */
    int *theta;      /* Running approximation of orbit partition */
    int *thsize;     /* Size of cells in theta, defined for mcrs */
    int *thnext;     /* Next rep in list (circular list) */
    int *thprev;     /* Previous rep in list */
    int *threp;      /* First rep for a given cell front */
    int *thfront;    /* The cell front associated with this rep */

    /* Search: split record */
    int *splitwho;   /* List of where splits occurred */
    int *splitfrom;  /* List of cells which were split */
    int *splitlev;   /* Where splitwho/from begins for each level */
    int nsplits;     /* Number of splits at this point */

    /* Search: differences from leftmost */
    char *diffmark;  /* Marked for diff labels */
    int *diffs;      /* List of diff labels */
    int *difflev;    /* How many labels diffed at each level */
    int ndiffs;      /* Current number of diffs */
    int *undifflev;  /* How many diff labels fixed at each level */
    int nundiffs;    /* Current number of diffs in singletons (fixed) */
    int *unsupp;     /* Inverted diff array */
    int *specmin;    /* Speculated mappings */
    int *pairs;      /* Not-undiffed diffs that can make two-cycles */
    int *unpairs;    /* Indices into pairs */
    int npairs;      /* Number of such pairs */
    int *diffnons;   /* Diffs that haven't been undiffed */
    int *undiffnons; /* Inverse of that */
    int ndiffnons;   /* Number of such diffs */

    /* Polymorphic functions */
    saucy_consumer *consumer;
    int (*split)(struct saucy *, struct coloring *, int, int);
    int (*is_automorphism)(struct saucy *);
    int (*ref_singleton)(struct saucy *, struct coloring *, int);
    int (*ref_nonsingle)(struct saucy *, struct coloring *, int);

     /* Statistics structure */
    struct saucy_stats *stats;
}; //}}} END saucy


struct saucy_stats
{ //{{{
    double grpsize_base;
    int grpsize_exp;
    int levels;
    int nodes;
    int bads;
    int gens;
    int support;
}; //}}} END saucy_stats

struct saucy_graph
{ //{{{
    int n;
    int e;
    int w;
    int *adj;
    int *edg;
    int *wght;
}; //}}} END saucy_graph

struct saucy* saucy_alloc( int n, int w, int directed,
                           const struct saucy_graph *g,
                           saucy_consumer *consumer,
                           void *arg,
                           struct saucy_stats *stats );

void saucy_search( struct saucy *s, const int *colors );

void saucy_free(struct saucy *s);

// Functions for warming up theta
void update_theta( struct saucy *s );
int find_representative( int k, int *theta );
void multiply_index( struct saucy *s, int k );
#endif
