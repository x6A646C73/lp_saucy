#include <cstdlib>
#include <iostream>
#include <cmath>
#include <climits>
#include <ctime>
#include <cstring>
#include <cctype>
#include <csignal>

// Saucy headers
#include "saucy.h"
#include "amorph.h"
#include "util.h"
#include "platform.h"

//NOTE: should be possible to do with original saucy structs afterall
//      think more on how to do this, it would be preferable I suppose
//TODO: use original saucy, this seems too slow, figure out how to fix...
//      just make a check at the DCCOUNT part that wght > 0
//      also, make edge_prob emulate edg and wght in structure?
using namespace std;

static sig_atomic_t timeout_flag = 0; /* Has the alarm gone off yet? */
static int quiet_mode = 1;  /* Don't output automorphisms */
static int K = 0;
static int LP = 0;

struct edge_t
{
    int u;
    int v;
    int idxu;
    int idxv;
    int wght;
};

struct e_t
{
    int idx; //into 2D edge_prob array
    int wght;
};

struct orbit_t {
    int orbit;
    int *edges; //of length K
    int len;
};

int numEdges;
struct edge_t *edge_data;

// Saucy structs
struct saucy *s;
struct saucy_stats stats;
struct amorph_graph *ng;
static char *marks; // "Bit" vector for printing

//****************************************************************************80
// BEGIN UTILITY METHODS
//****************************************************************************80
static void arg_lp( char *arg ) { LP = 1; }

static void arg_budget( char *arg )
{
    K = atoi( arg );
    if( K < 0 ) die( "budget must be non-negative" );
}

static void arg_help( char *arg );

static struct option options[] = {
    { "budget", 'b', "N", arg_budget,
    "edge budget for the almost symmetries (default 0)" },
    { "lp", 'l', 0, arg_lp,
    "flag indicating whether the input file is an lp (default 0)" },
    { "help", 0, 0, arg_help,
    "output this help message" },
    { 0, 0, 0, 0, 0 }
};

static void arg_help( char *arg )
{
    printf( "usage: saucy [OPTION]... FILE\n" );
    print_options( options );
    exit( 0 );
}

static int on_automorphism( int n, const int *gamma, int k, int *support, void *arg )
{
    struct amorph_graph *g = (struct amorph_graph *)arg;
    if( !quiet_mode )
    {
        qsort_integers( support, k );
        g->consumer( n, gamma, k, support, g, marks );
    }
    return !timeout_flag;
}

//TODO: make this a macro
inline int min( int a, int b )
{
    return (a<b) ? a : b;
}
//****************************************************************************80
// END UTILITY METHODS
//****************************************************************************80


//****************************************************************************80
// BEGIN DYNAMIC PROGRAMMING METHODS
//****************************************************************************80

// EVALUATE implements the user-defined valuation function
int evaluate( int *edges, int len )
{
    int i, idxu, idxv;
    
    // Remove weights from population member
    for( i = 0; i < len; i++ )
    {
        idxu = edge_data[edges[i]].idxu;
        idxv = edge_data[edges[i]].idxv;
        ng->sg.wght[idxu] -= edge_data[edges[i]].wght;
        ng->sg.wght[idxv] -= edge_data[edges[i]].wght;
    }
    
    // Call saucy
    saucy_search( s, &ng->sg, 0, ng->colors, on_automorphism,
                  ng, &stats );
    
    // Return the weights to the graph
    for( i = 0; i < len; i++ )
    {
        idxu = edge_data[edges[i]].idxu;
        idxv = edge_data[edges[i]].idxv;
        ng->sg.wght[idxu] += edge_data[edges[i]].wght;
        ng->sg.wght[idxv] += edge_data[edges[i]].wght;
    }
    
    return stats.orbits;
}

// INITIALIZE initializes the genes within the variables bounds. 
//void initialize( string filename, int &seed )
void initialize( string filename )
{
    int i, j, k, u, v, temp, idx, start;
    double rnum;
    
    ng = (LP) ? read_lp( filename.c_str() ) : read_graph( filename.c_str() );
    if( !ng ) die( "unable to read input file" );
    
    numEdges = ng->sg.e;
    edge_data = (edge_t *)malloc( numEdges*sizeof(edge_t) );
    
    // Allocate some memory to facilitate printing
    marks = (char *)calloc( ng->sg.n, sizeof(char) );
    if( !marks ) die( "out of memory" );
    
    // Allocate saucy space
    s = (struct saucy*)saucy_alloc( ng->sg.n, ng->sg.w );
    if( s == NULL ) die( "saucy initialization failed" );
    
    for( u = 0, idx = 0; u < ng->sg.n; u++ )
    {
        for( j = ng->sg.adj[u]; j != ng->sg.adj[u+1]; j++ )
        {
            v = ng->sg.edg[j];
            if( v > u )
            {
                for( k = ng->sg.adj[v]; k != ng->sg.adj[v+1]; k++ )
                {
                    if( ng->sg.edg[k] == u )
                    {
                        edge_data[idx].idxu = k;
                        break;
                    }
                }
                edge_data[idx].idxv = j;
                edge_data[idx].u = u;
                edge_data[idx].v = v;
                edge_data[idx].wght = ng->sg.wght[j];
                idx++;
            }
        }
    }
    
    return;
}

//TODO: associate wghts with edges below
void findAlmostDYN( struct orbit_t *ret_orb )
{
    int e, k, i, new_orb;
    struct orbit_t A[K+1][numEdges+1];
    
    // Initialize the table A[][]
    for( k = 0; k <= K; k++ )
    {
        for( e = 0; e <= numEdges; e++ )
        {
            A[k][e].edges = (int*)malloc( K*sizeof(int) );
            A[k][e].len = 0;
        }
    }
    
    // Build table A[][] in bottom up manner
    for( k = 0; k <= K; k++ )
    {
        for( e = 0; e <= numEdges; e++ )
        {
            if( e == 0 || k == 0 )
            {
                A[k][e].orbit = evaluate( A[k][e].edges, 0 );
            }
            else
            {
                A[k-1][e-1].edges[A[k-1][e-1].len] = e-1; //since edge ids start at 0
                A[k-1][e-1].len++;
                new_orb = evaluate( A[k-1][e-1].edges, A[k-1][e-1].len );
                
                if( new_orb >= A[k][e-1].orbit )
                {
                    //copy over A[k][e-1] stuff
                    A[k][e].orbit = A[k][e-1].orbit;
                    A[k][e].len = A[k][e-1].len;
                    for( i = 0; i < A[k][e].len; i++ )
                    {
                        A[k][e].edges[i] = A[k][e-1].edges[i];
                    }
                }
                else
                {
                    //copy over A[k-1][e-1] stuff
                    A[k][e].orbit = new_orb;
                    A[k][e].len = A[k-1][e-1].len;
                    for( i = 0; i < A[k][e].len; i++ )
                    {
                        A[k][e].edges[i] = A[k-1][e-1].edges[i];
                    }
                }
                A[k-1][e-1].len--;
            }
        }
    }
    
    // Store the best soln
    ret_orb->len = A[K][numEdges].len;
    ret_orb->orbit = A[K][numEdges].orbit;
    for( i = 0; i < A[K][numEdges].len; i++ )
    {
        ret_orb->edges[i] = A[K][numEdges].edges[i];
    }
    
    // Clean up
    for( k = 0; k <= K; k++ )
    {
        for( e = 0; e <= numEdges; e++ )
        {
            free( A[k][e].edges );
        }
    }
}

//****************************************************************************80
// END DYNAMIC PROGRAMMING METHODS
//****************************************************************************80


//****************************************************************************80
//    MAIN supervises the genetic algorithm.
//****************************************************************************80
int main( int argc, char *argv[] )
{
    int i, c, idx;
    char *filename = NULL;
    struct orbit_t ret_orb;
    
    // Option handling 
    parse_arguments( &argc, &argv, options );
    if( argc < 1 ) die( "missing filename" );
    if( argc > 1 ) die( "trailing arguments" );
    filename = *argv;
    
    // Initialize data
    initialize( filename );
    ret_orb.edges = (int*)malloc( K*sizeof(int) );
    ret_orb.len = 0;
    ret_orb.orbit = 0;
    
    // Run dynamic programming algorithm
    findAlmostDYN( &ret_orb );

    // Print the best solution found
    cout << "\n";
    cout << "  Best member\n";
    cout << "\n";
    for( i = 0; i < ret_orb.len; i++ )
    {
        idx = ret_orb.edges[i];
        cout << "  edges[" << i << "] = " << idx << "( "
             << edge_data[idx].u << ", " << edge_data[idx].v << " )" << endl;
    }
    cout << "\n";
    cout << "  Best fitness = " << ret_orb.orbit << "\n";
    
    // Clean up and terminate
    saucy_free( s );
    ng->free( ng );
    free( marks );
    free( edge_data );
    free( ret_orb.edges );
    
    return 0;
}
