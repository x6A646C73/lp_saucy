#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "lp_amorph.h"
#include "warmup_theta.h"

static int check_mapping( struct lp_saucy *s, const int *adj, const int *edg,
                          const int *wght, int k )
{ //{{{
    int i, gk, ret = 1;
    
    /* Mark gamma of neighbors */
    for( i = adj[k]; i != adj[k+1]; ++i )
    { //{{{
        s->stuff[s->gamma[edg[i]]] = 1;
        /* wstuff should only need n spots since a vertex */
        /* can connect to at most n vertices */
        s->wstuff[s->gamma[edg[i]]] = wght[i];
    } //}}}
    
    /* Check neighbors of gamma */
    gk = s->gamma[k];
    for( i = adj[gk]; ret && i != adj[gk+1]; ++i )
    { //{{{
        ret = s->stuff[edg[i]];
        /* TODO: verify that this is a valid test for weight data */
        ret = ret && (wght[i] == s->wstuff[edg[i]]);
    } //}}}
    
    /* Clear out bit vector before we leave */
    for( i = adj[k]; i != adj[k+1]; ++i )
    { //{{{
        s->stuff[s->gamma[edg[i]]] = 0;
        s->wstuff[s->gamma[edg[i]]] = 0;
    } //}}}
    
    return ret;
} //}}} END check_mapping

static int is_undirected_automorphism( struct lp_saucy *s )
{ //{{{
    int i, j;
    
    for( i = 0; i < s->ndiffs; ++i )
    { //{{{
        j = s->unsupp[i];
        if( !check_mapping( s, s->adj, s->edg, s->wght, j ) ) return 0;
    } //}}}
    
    return 1;
} //}}} END is_undirected_automorphism

static int is_directed_automorphism( struct lp_saucy *s )
{ //{{{
    int i, j;
    
    for( i = 0; i < s->ndiffs; ++i )
    { //{{{
        j = s->unsupp[i];
        if( !check_mapping( s, s->adj, s->edg, s->wght, j ) ) return 0;
        if( !check_mapping( s, s->dadj, s->dedg, s->dwght, j ) ) return 0;
    } //}}}
    
    return 1;
} //}}} END is_directed_automorphism

static int find_representative( int k, int *theta )
{ //{{{
    int rep, tmp;
    
    /* Find the minimum cell representative */
    for( rep = k; rep != theta[rep]; rep = theta[rep] );
    
    /* Update all intermediaries */
    while( theta[k] != rep )
    { //{{{
        tmp = theta[k]; theta[k] = rep; k = tmp;
    } //}}}
    
    return rep;
} //}}} END find_representative

static void update_theta( struct lp_saucy *s )
{ //{{{
    int i, k, x, y, tmp;
    
    for( i = 0; i < s->ndiffs; ++i )
    { //{{{
        k = s->unsupp[i];
        x = find_representative( k, s->theta );
        y = find_representative( s->gamma[k], s->theta );
        
        if( x != y )
        { //{{{
            if( x > y )
            { //{{{
                tmp = x;
                x = y;
                y = tmp;
            } //}}}
            s->theta[y] = x;
            s->thsize[x] += s->thsize[y];
            
            s->thnext[s->thprev[y]] = s->thnext[y];
            s->thprev[s->thnext[y]] = s->thprev[y];
            s->threp[s->thfront[y]] = s->thnext[y];
        } //}}}
    } //}}}
} //}}} END update_theta

static void multiply_index( struct saucy *s, int k )
{ //{{{
    if( ( s->stats->grpsize_base *= k ) > 1e10 )
    { //{{{
        s->stats->grpsize_base /= 1e10;
        s->stats->grpsize_exp += 10;
    } //}}}
} //}}} END multiply_index

//TODO: consider something other than text file?
int lp_warmup_theta( char *filename, struct lp_saucy *s, int directed, 
                     const struct saucy_graph *g )
{ //{{{
    FILE *file;
    int num_gens=0, i=0, j=0, temp;
    int rep, repsize;
    
    /* Save graph information */
    s->n = g->n;
    s->wcount = g->w;
    s->adj = g->adj;
    s->edg = g->edg;
    s->wght = g->wght;
    s->dadj = g->adj + g->n + 1;
    s->dedg = g->edg + g->e;
    s->dwght = g->wght + g->e;
    
    // generators are the whole perm: (0 1 2)(3 4) as 1 2 0 4 3
    file = fopen( filename, "r" );
    fscanf( file, "%d", &num_gens );
    while( i < num_gens && !feof( file ) )
    { //{{{
        s->ndiffs = 0;
        for( j = 0; j < s->n; ++j )
        {   
            fscanf( file, "%d", &temp );
            if( temp != j ) 
            {
                s->unsupp[s->ndiffs] = j;
                ++(s->ndiffs);
            }
            s->gamma[j] = temp;
        }
        
        if( (!directed && is_undirected_automorphism( s )) || 
            (directed && is_directed_automorphism( s )) ) 
        { //{{{
            ++s->stats->gens;
            update_theta( s ); 
        } //}}}
        ++i;
    } //}}}
    fclose( file );
    
    i = 0;
    while( i < s->n )
    { //{{{
        rep = find_representative( i, s->theta ); 
        repsize = s->thsize[rep];
        multiply_index( s, repsize );
        i += repsize;
    } //}}}
    
    /* Normalize group size */
    while( s->stats->grpsize_base >= 10.0 )
    { //{{{
        s->stats->grpsize_base /= 10;
        ++s->stats->grpsize_exp;
    } //}}}
    
    //for( i = 0; i < s->n; ++i ) s->gamma[i] = i;
    
    return num_gens;
} //}}} END warmup_theta

static int *ints(int n) { return malloc(n * sizeof(int)); }
static char *bits(int n) { return calloc(n, sizeof(char)); }

struct lp_saucy* lp_saucy_alloc( int n, int w, struct saucy_stats *stats )
{ //{{{
    int i;
    struct lp_saucy *s = (struct lp_saucy *)malloc( sizeof(struct lp_saucy) );
    if( s == NULL ) return NULL;
    
    // Allocate data {{{
    s->stuff = bits( n+1 );
    s->wstuff = ints( n+1 );
    s->gamma = ints( n );
    s->theta = ints( n );
    s->thsize = ints( n );
    s->unsupp = ints( n );
    s->thnext = ints( n );
    s->thprev = ints( n );
    s->threp = ints( n );
    s->thfront = ints( n );
    //}}}
    
    if( s->stuff && s->wstuff && s->gamma && s->theta && s->thsize
        && s->unsupp && s->thnext && s->thprev && s->threp && s->thfront )
    {
        /* The initial orbit partition is discrete */
        for( i = 0; i < n; ++i ) s->theta[i] = i;
        
        /* The initial permutation is the identity */
        for( i = 0; i < n; ++i ) s->gamma[i] = i;
        
        /* Initially every cell of theta has one element */
        for( i = 0; i < n; ++i ) s->thsize[i] = 1;
        
        /* Every theta rep list is singleton */
        for( i = 0; i < n; ++i ) s->thprev[i] = s->thnext[i] = i;
        
        /* Initialize stats */
        s->stats = stats;
        s->stats->grpsize_base = 1.0;
        s->stats->grpsize_exp = 0;
        s->stats->nodes = 1;
        s->stats->bads = s->stats->gens = s->stats->support = 0;
        
        return s;
    else
    { //{{{
        lp_saucy_free( s );
        return NULL;
    } //}}}
} //}}} END lp_saucy_alloc

void lp_saucy_free( struct lp_saucy *s )
{ //{{{
    free( s->thfront );
    free( s->threp );
    free( s->thnext );
    free( s->thprev );
    free( s->thsize );
    free( s->unsupp );
    free( s->theta );
    free( s->gamma );
    free( s->stuff );
    free( s->wstuff );
    free( s );
} //}}} END saucy_free
