#include <cstdlib>
#include <map>
#include <cstring>
#include <cstdio>
#include "CoinLpIO.hpp"
#include "CoinMpsIO.hpp"
//#include "util.h" //gives die, warn, etc.
#include "lp_amorph.h"

typedef std::map<double,int> colorMap;

static int parse_filename( char *filename )
{ //{{{
    int res = -1; 
    const char *dot = strrchr( filename, '.' );
    if( !dot || dot == filename ) res = -1; 
    else if( strcmp( dot+1, "mps" ) == 0 ) res = 0;
    else if( strcmp( dot+1, "lp" ) == 0 ) res = 1;
    return res;
} //}}} END parse_filename

static void readLP( char *file, const double *obj, const double *rhs, int &vars,
             int &cons, const CoinPackedMatrix *mat )
{ //{{{
    CoinLpIO prob;
    prob.readLp( file );
    obj = prob.getObjCoefficients(); //gives vars their color
    rhs = prob.getRightHandSide();   //holds color for other side of bipartite
    vars = prob.getNumCols();        //number of variables
    cons = prob.getNumRows();        //number of constraints
    mat = prob.getMatrixByCol();     //contains weights for edges
} //}}} END readLP

static void readMPS( char *file, const double *obj, const double *rhs, int &vars,
              int &cons, const CoinPackedMatrix *mat )
{ //{{{
    CoinMpsIO prob;
    prob.readMps( file, "mps" );
    obj = prob.getObjCoefficients(); //gives vars their color
    rhs = prob.getRightHandSide();   //holds color for other side of bipartite
    vars = prob.getNumCols();        //number of variables
    cons = prob.getNumRows();        //number of constraints
    mat = prob.getMatrixByCol();     //contains weights for edges
} //}}} END readMPS

//From: saucyio.c
static int lp_init_fixadj1( int n, int *adj )
{ //{{{
    int val, sum, i;
    
    /* Translate adj values to real locations */
    val = adj[0]; sum = 0; adj[0] = 0;
    for (i = 1; i < n; ++i) {
        sum += val;
        val = adj[i];
        adj[i] = sum;
    }
    
    return sum + val;
} //}}} END lp_init_fixadj1

//From: saucyio.c
static void lp_init_fixadj2( int n, int e, int *adj )
{ //{{{
    int i;
    
    /* Translate again-broken sizes to adj values */
    for (i = n-1; i > 0; --i) {
        adj[i] = adj[i-1];
    }
    
    adj[0] = 0;
    adj[n] = e;
} //}}} END lp_init_fixadj2

//From: saucyio.c
static void lp_amorph_print_automorphism( int n, const int *gamma, int nsupp,
                                          const int *support,
                                          struct lp_amorph_graph *g, char *marks )
{ //{{{
    int i, j, k;
    
    /* We presume support is already sorted */
    for( i = 0; i < nsupp; ++i )
    { //{{{
        k = support[i];
        
        /* Skip elements already seen */
        if( marks[k] ) continue;
        
        /* Start an orbit */
        marks[k] = 1;
        printf( "(%d", k );
        
        /* Mark and notify elements in this orbit */
        for( j = gamma[k]; j != k; j = gamma[j] )
        {
            marks[j] = 1;
            printf( " %d", j );
        }
        
        /* Finish off the orbit */
        putchar( ')' );
    } //}}}
    putchar( '\n' );
    
    /* Clean up after ourselves */
    for( i = 0; i < nsupp; ++i )
        marks[support[i]] = 0;
} //}}} END lp_amorph_print_automorphism

//From: saucyio.c
void lp_amorph_graph_free( struct lp_amorph_graph *g )
{ //{{{
    free( g->sg.adj );
    free( g->sg.edg );
    free( g->colors );
    free( g );
} //}}} END lp_amorph_graph_free

struct lp_amorph_graph* lp_amorph_read_build( char *filename )
{ //{{{
    // Variable declarations {{{
    int n, e;
    int i, j, k, *aout, *eout, *colors;
    int tempk, tempj, tempw;
    int w, *wout;
    int ndx, col;
    const double *elm;
    const int *ind, *str;
    struct lp_amorph_graph *g = NULL;
    colorMap varC, cosC, conC;
    std::pair<std::map<double,int>::iterator,bool> ret;
    int c1=0, c2=0, c3=0;
    // ILP variables
    int vars, cons, lp=1;
    const double *obj, *rhs;
    const CoinPackedMatrix *mat;
    //}}}
    
    lp = parse_filename( filename );
    if( lp == 1 ) readLP( filename, obj, rhs, vars, cons, mat );
    else if( lp == 0 ) readMPS( filename, obj, rhs, vars, cons, mat );
    else
    { //{{{
        fprintf( stderr, "ERROR: must provide lp or mps file. Terminating..." );
        return NULL;
    } //}}}
    
    n = vars+cons;
    e = mat->getNumElements();
    
    // Allocate everything {{{
    g = ( struct lp_amorph_graph* ) malloc( sizeof( struct lp_amorph_graph ) );
    aout = ( int* ) calloc( (n+1), sizeof(int) );
    eout = ( int* ) malloc( 2 * e * sizeof(int) );
    wout = ( int* ) malloc( 2 * e * sizeof(int) );
    wlist = ( int* ) malloc( 2 * e * sizeof(int) );
    colors = ( int* ) malloc( n * sizeof(int) );
    if (!g || !aout || !eout || !colors)
    { //{{{
        free(g); free(aout); free(eout);
        free(wout); free(wlist);    free(colors);
        return NULL;
    } //}}}
    //}}}
    
    // Assign everything {{{
    g->sg.n    =      n;
    g->sg.e    =      e;
    g->sg.w    =      0;
    g->sg.adj  =   aout;
    g->sg.edg  =   eout;
    g->sg.wght =   wout;
    g->colors  = colors;
    
    elm  =     mat->getElements();
    ind  =      mat->getIndices();
    str  = mat->getVectorStarts();
    ///}}}
    
    for( i = 0; i < vars; i++ )
    { //{{{
        ret = varC.insert( std::pair<double,int>( obj[i], c1 ) );
        if( ret.second == true ) c1++;
    } //}}}
    for( i = 0; i < cons; i++ )
    { //{{{
        ret = cosC.insert( std::pair<double,int>( rhs[i], c1+c2 ) );
        if( ret.second == true ) c2++;
    } //}}}
    
    col = 0;
    j = 1;
    for( i = 0; i < e; i++ )
    { //{{{
        while( i == str[j] )
        { //{{{
            col++;
            j++;
        } //}}}
        
        ++aout[col]; ++aout[vars+ind[i]];
        
        ret = conC.insert( std::pair<double,int>( elm[i], c1+c2+c3 ) );
        if( ret.second == true ) c3++;
    } //}}}
    g->sg.w = c3;
    lp_init_fixadj1(n, aout);

    for( i = 0; i < v; ++i )
        colors[i] = varC.find( obj[i] )->second;
    for( i = 0; i < c; ++i )
        colors[i] = cosC.find( rhs[i] )->second;
    
    col = 0;
    j = 1;
    for( i = 0; i < e; i++ )
    { //{{{
        while( i == str[j] )
        { //{{{
            col++;
            j++;
        } //}}}
        
        tempj       =                 aout[col]++; 
        tempk       =         aout[vars+ind[i]]++;
        eout[tempj] =                 vars+ind[i]; 
        eout[tempk] =                         col;
        tempw       = conC.find( elm[i] )->second;
        wout[tempj] =                           w; 
        wout[tempk] =                           w;
    } //}}}
    lp_init_fixadj2(n, 2 * e, aout);
    
    // Assign the functions {{{
    g->consumer = lp_amorph_print_automorphism;
    g->free     =         lp_amorph_graph_free;
    g->stats    =                         NULL; 
    // }}}
    
    return g;
} //}}} END amorph_build

//TODO: fix this and automorph functions to work here
static int check_mapping( struct lp_saucy *s, const int *adj, const int *edg, const int *wght, int k )
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

int is_undirected_automorphism( struct lp_saucy *s )
{ //{{{
    int i, j;
    
    for( i = 0; i < s->ndiffs; ++i )
    { //{{{
        j = s->unsupp[i];
        if( !check_mapping( s, s->adj, s->edg, s->wght, j ) ) return 0;
    } //}}}
    
    return 1;
} //}}} END is_undirected_automorphism

int is_directed_automorphism( struct lp_saucy *s )
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

static int lp_find_representative( int k, int *theta )
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
} //}}} END lp_find_representative

static void lp_multiply_index( struct lp_saucy *s, int k )
{ //{{{
    if( ( s->stats->grpsize_base *= k ) > 1e10 )
    { //{{{
        s->stats->grpsize_base /= 1e10;
        s->stats->grpsize_exp += 10;
    } //}}}
} //}}} END lp_multiply_index

// can optimize later...
static void lp_update_theta( struct lp_saucy *s )
{ //{{{
    int i, x, y, tmp;
    
    for( i = 0; i < s->n; ++i )
    { //{{{
        if( s->gamma[i] == i ) continue;
        
        x = lp_find_representative( i, s->theta );
        y = lp_find_representative( s->gamma[i], s->theta );
        
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

//TODO: consider something other than text file?
//TODO: use original saucy functions for this?
int lp_warmup_theta( char *filename, struct lp_saucy *s )
{ //{{{
    std::ifstream file( filename );
    int num_gens = i = 0;
    
    // generators are the whole perm: (0 1 2)(3 4) as 1 2 0 4 3
    file >> num_gens;
    while( i < num_gens && !file.eof )
    { //{{{
        for( j = 0; j < s->n; ++j ) file >> s->gamma[j];
        
        if( s->lp_is_automorphism( s ) )
        { //{{{
            ++s->stats->gens;
            //TODO: double check that this works
            //TODO: write a small test program for this with 1 2 0 4 3 as input
            //      then can determine what thnext, etc. hold as well...
            lp_update_theta( s );
        } //}}}
        ++i;
    } //}}}
    
    i = 0;
    while( i < s->n )
    { //{{{
        rep = lp_find_representative( i, s->theta );
        repsize = s->thsize[rep];
        lp_multiply_index( s, repsize );
        i += repsize;
    } //}}}
    
    /* Normalize group size */
    while( s->stats->grpsize_base >= 10.0 )
    { //{{{
        s->stats->grpsize_base /= 10;
        ++s->stats->grpsize_exp;
    } //}}}

    for( i = 0; i < s->n; ++i ) s->gamma[i] = i;
    
    return num_gens
} //}}} END warmup_theta

static int *ints( int n ) { return malloc( n * sizeof(int) ); }
static char *bits( int n ) { return calloc( n, sizeof(char) ); }

struct lp_saucy *lp_saucy_alloc( int n )
{ //{{{
    struct lp_saucy *s = (struct lp_saucy *)malloc( sizeof(struct lp_saucy) );
    if( s == NULL ) return NULL;
    
    s->stuff = bits( n+1 );
    s->wstuff = ints( n+1 );
    s->gamma = ints( n );
    s->theta = ints( n );
    s->thsize = ints( n );
    s->thnext = ints( n );
    s->thprev = ints( n );
    s->threp = ints( n );
    s->thfront = ints( n );
    
    if( s->stuff && s->wstuff && s->gamma && s->theta && s->thsize
        && s->thnext && s->thprev && s->threp && s->thfront )
    { //{{{
        return s;
    } //}}}
    else
    { //{{{
        saucy_free( s );
        return NULL;
    } //}}}
} //}}} END lp_saucy_alloc

void lp_saucy_free( struct saucy *s )
{ //{{{
    free( s->thfront );
    free( s->threp );
    free( s->thnext );
    free( s->thprev );
    free( s->thsize );
    free( s->theta );
    free( s->gamma );
    free( s->stuff );
    free( s->wstuff );
    free( s );
} //}}} END lp_saucy_free
