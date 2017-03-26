#include <cstdlib>
#include <map>
#include <cstring>
#include <cstdio>
#include "CoinLpIO.hpp"
#include "CoinMpsIO.hpp"
#include "util.h"
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

static void readMPS( char *file, const double *&obj, const double *&rhs, int &vars,
                     int &cons, const CoinPackedMatrix *&mat )
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
    int i, j, *aout, *eout, *colors;
    int tempk, tempj, tempw;
    int *wout;
    int col;
    const double *elm;
    const int *ind, *str;
    struct lp_amorph_graph *g = NULL;
    colorMap varC, cosC, conC;
    std::pair<std::map<double,int>::iterator,bool> ret;
    int c1=0, c2=0, c3=0;
    // ILP variables
    int vars, cons, lp=1;
    const double *obj=NULL, *rhs=NULL;
    const CoinPackedMatrix *mat=NULL;
    //}}}
    
    lp = parse_filename( filename );
    if( lp == 1 ) readLP( filename, obj, rhs, vars, cons, mat );
    else if( lp == 0 ) readMPS( filename, obj, rhs, vars, cons, mat );
    else die( "must provide lp or mps file" );
    if( mat == NULL || obj == NULL || rhs == NULL )
        die( "failed to read problem file" );
    
    n = vars+cons;
    e = mat->getNumElements();
    
    // Allocate everything {{{
    g = ( struct lp_amorph_graph* ) malloc( sizeof( struct lp_amorph_graph ) );
    aout = ( int* ) calloc( (n+1), sizeof(int) );
    eout = ( int* ) malloc( 2 * e * sizeof(int) );
    wout = ( int* ) malloc( 2 * e * sizeof(int) );
    //wlist = ( int* ) malloc( 2 * e * sizeof(int) );
    colors = ( int* ) malloc( n * sizeof(int) );
    if (!g || !aout || !eout || !colors)
    { //{{{
        free(g); free(aout); free(eout);
        free(wout); free(colors); //free(wlist);
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
    
    for( i = 0; i < vars; ++i )
        colors[i] = varC.find( obj[i] )->second;
    for( i = 0; i < cons; ++i )
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
        wout[tempj] =                       tempw; 
        wout[tempk] =                       tempw;
    } //}}}
    lp_init_fixadj2(n, 2 * e, aout);
    
    // Assign the functions {{{
    g->consumer = lp_amorph_print_automorphism;
    g->free     =         lp_amorph_graph_free;
    g->stats    =                         NULL; 
    // }}}
    
    return g;
} //}}} END amorph_build
