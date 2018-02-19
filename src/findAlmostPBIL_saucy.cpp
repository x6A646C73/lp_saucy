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
static int seed = 123456789;
static int POPSIZE = 100;
static int MAXGENS = 1000;
static double PMUTATION = 0.02;
static double MUTSHIFT = 0.05;
static double LEARNRATE = 0.1;
static int K = 0;
static int LP = 0;

struct edge_t
{
    //int u;
    //int v;
    int idxu;
    int idxv;
    double prob;
};

struct e_t
{
    int idx; //into 2D edge_prob array
    int wght;
};

struct gene_t
{
    struct e_t *gene; //hold edge data
    int fitness;
    int len;
};

int numEdges;
struct gene_t *pop;
double *edge_prob;
struct edge_t *edge_prob;

// Saucy structs
struct saucy *s;
struct saucy_stats stats;
struct amorph_graph *ng;
static char *marks; // "Bit" vector for printing

//****************************************************************************80
// BEGIN UTILITY METHODS
//****************************************************************************80
static void arg_lp( char *arg ) { LP = 1; }

static void arg_popsize( char *arg )
{
    POPSIZE = atoi( arg );
    if( POPSIZE <= 0 ) die( "population size must be positive" );
}

static void arg_maxgens( char *arg )
{
    MAXGENS = atoi( arg );
    if( MAXGENS <= 0 ) die( "number of generations must be positive" );
}

static void arg_mutation( char *arg )
{
    PMUTATION = atof( arg );
    if( PMUTATION <= 0.0 || PMUTATION >=1.0 ) die( "mutation rate  must be in (0,1)" );
}

static void arg_mutshift( char *arg )
{
    MUTSHIFT = atof( arg );
    if( MUTSHIFT <= 0.0 || MUTSHIFT >=1.0 ) die( "mutation shift must be in (0,1)" );
}

static void arg_learnrate( char *arg )
{
    LEARNRATE = atof( arg );
    if( LEARNRATE <= 0.0 || LEARNRATE >=1.0 ) die( "learning rate must be in (0,1)" );
}

static void arg_budget( char *arg )
{
    K = atoi( arg );
    if( K < 0 ) die( "budget must be non-negative" );
}

static void arg_seed( char *arg )
{
    seed = atoi( arg );
    if( seed <= 0 ) die( "seed must be positive" );
}

static void arg_help( char *arg );

static struct option options[] = {
    { "popsize", 'p', "N", arg_popsize,
    "population size for the PBIL (default 100)" },
    { "generations", 'g', "N", arg_maxgens,
    "number of generations to run (default 1000)" },
    { "mutation", 'm', "N", arg_mutation,
    "mutation rate for the PBIL (default 0.02)" },
    { "mutshift", 't', "N", arg_mutshift,
    "mutation shift amount for the PBIL (default 0.05)" },
    { "rate", 'r', "N", arg_learnrate,
    "learning rate for the PBIL (default 0.1)" },
    { "budget", 'b', "N", arg_budget,
    "edge budget for the almost symmetries (default 0)" },
    { "seed", 's', "N", arg_seed,
    "seed for random number generator (default 123456789)" },
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

// I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//int i4_uniform_ab( int a, int b, int &seed )
int i4_uniform_ab( int a, int b )
{
    int c;
    const int i4_huge = 2147483647;
    int k;
    float r;
    int value;
    
    if( seed == 0 )
    {
        cerr << "\n";
        cerr << "I4_UNIFORM_AB - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit ( 1 );
    }
    
    //  Guarantee A <= B.
    if( b < a )
    {
        c = a;
        a = b;
        b = c;
    }
    
    k = seed / 127773;
    
    seed = 16807 * ( seed - k * 127773 ) - k * 2836;
    
    if( seed < 0 )
        seed = seed + i4_huge;
    
    r = ( float ) ( seed ) * 4.656612875E-10;
    
    //  Scale R to lie between A-0.5 and B+0.5.
    r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
        +         r   * ( ( float ) b + 0.5 );
    
    //  Use rounding to convert R to an integer between A and B.
    value = round ( r );
    
    //  Guarantee A <= VALUE <= B.
    if( value < a )
        value = a;
    if ( b < value )
        value = b;
    
    return value;
}

// R8_UNIFORM_AB returns a scaled pseudorandom R8.
//double r8_uniform_ab( double a, double b, int &seed )
double r8_uniform_ab( double a, double b )
{
    int i4_huge = 2147483647;
    int k;
    double value;
    
    if( seed == 0 )
    {
        cerr << "\n";
        cerr << "R8_UNIFORM_AB - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit ( 1 );
    }
    
    k = seed / 127773;
    
    seed = 16807 * ( seed - k * 127773 ) - k * 2836;
    
    if ( seed < 0 )
        seed = seed + i4_huge;
    
    value = ( double ) ( seed ) * 4.656612875E-10;
    value = a + ( b - a ) * value;
    
    return value;
}

// ISIN determines if an edge is in a list.
bool isIn( int e, struct e_t *arr, int size )
{
    for( int i = 0; i < size; i++ )
    {
        if( e == arr[i].idx ) return true;
    }
    
    return false;
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

//****************************************************************************80
// BEGIN UTILITY METHODS
//****************************************************************************80


//****************************************************************************80
// BEGIN PBIL METHODS
//****************************************************************************80

// EVALUATE implements the user-defined valuation function
//TODO: make this work with more generic weight budget, i.e. not removing entire
//      edge, just some of the weight
void evaluate( int &idxB, int &idxW )
{
    int mem, i, idxu, idxv;
    int fitB = INT_MAX, fitW = -1;
    
    for( mem = 0; mem < POPSIZE; mem++ )
    {
        // Remove weights from population member
        for( i = 0; i < pop[mem].len; i++ )
        {
            idxu = edge_prob[pop[mem].gene[i].idx].idxu;
            idxv = edge_prob[pop[mem].gene[i].idx].idxv;
            ng->sg.wght[idxu] -= pop[mem].gene[i].wght;
            ng->sg.wght[idxv] -= pop[mem].gene[i].wght;
        }
        
        //TODO: get this right
        // Call saucy
        saucy_search( s, &ng->sg, 0, ng->colors, on_automorphism,
                      ng, &stats );
        pop[mem].fitness = stats.orbits;
        if( pop[mem].fitness < fitB )
        {
            fitB = pop[mem].fitness;
            idxB = mem;
        }
        if( pop[mem].fitness > fitW )
        {
            fitW = pop[mem].fitness;
            idxW = mem;
        }
        
        // Return the weights to the graph
        for( i = 0; i < pop[mem].len; i++ )
        {
            idxu = edge_prob[pop[mem].gene[i].idx].idxu;
            idxv = edge_prob[pop[mem].gene[i].idx].idxv;
            ng->sg.wght[idxu] -= pop[mem].gene[i].wght;
            ng->sg.wght[idxv] -= pop[mem].gene[i].wght;
        }
    }
    
    return;
}

// UPDATE modifies the probabilities for the edges and number of edges
void update( int idxB, int idxW )
{
    int i, j;
    bool inB;
    
    //pop[idx].gene = { int u; int v; int wght; }
    //len(edge_prob) == ng->sg.e
    //for( i = 0; i < numEdges; i++ )
    //TODO: for this, no longer assume edge[i] is 0 or 1 but instead can be any
    //      value between wght[i]-K and wght[i]
    //      how is the learning done then?
    for( i = 0; i < numEdges; i++ )
    {
        inB = isIn( i, pop[idxB].gene, pop[idxB].len );
        edge_prob[i].prob = edge_prob[i].prob*( 1.0 - LEARNRATE ) +
                            ( inB ? LEARNRATE : 0 );
    }
    
    return;
}

// INITIALIZE initializes the genes within the variables bounds. 
//void initialize( string filename, int &seed )
void initialize( string filename )
{
    int i, j, k, v, temp, idx, start;
    double rnum;
    
    ng = (LP) ? read_lp( filename.c_str() ) : read_graph( filename.c_str() );
    if( !ng ) die( "unable to read input file" );
    
    pop = (struct gene_t*)malloc( (POPSIZE+1)*sizeof(struct gene_t) );
    numEdges = ng->sg.e;
    edge_prob = (double *)malloc( numEdges*sizeof(double) );
    
    // Allocate some memory to facilitate printing
    marks = (char *)calloc( ng->sg.n, sizeof(char) );
    if( !marks ) die( "out of memory" );
    
    // Allocate saucy space
    s = (struct saucy*)saucy_alloc( ng->sg.n, ng->sg.w );
    if( s == NULL ) die( "saucy initialization failed" );
    
    for( i = 0, idx = 0; i < ng->sg.n; i++ )
    {
        for( j = ng->sg.adj[i]; j != ng->sg.adj[i+1]; j++ )
        {
            v = ng->sg.edg[j];
            if( v > i )
            {
                for( k = ng->sg.adj[v]; k != ng->sg.adj[v+1]; k++ )
                {
                    if( ng->sg.edg[k] == i )
                    {
                        edge_prob[idx].idxu = k;
                        break;
                    }
                }
                edge_prob[idx].idxv = j;
                edge_prob[idx].prob = 0.5;
                idx++;
            }
        }
    }
    
    //  Initialize population members 
    for( i = 0; i < POPSIZE; i++ )
    {
        pop[i].fitness = 0;
        //TODO: may end up needing different size for gene
        pop[i].gene = (struct e_t*)malloc( K*sizeof(struct e_t) );
    }
    
    // Initialize the location that holds the best seen individual
    pop[POPSIZE].fitness = INT_MAX;
    pop[POPSIZE].gene = (struct e_t*)malloc( K*sizeof(struct e_t) );
    pop[POPSIZE].len = 0;
    
    return;
}

// GENERATE creates the next generation of genes
//void generate( int &seed )
//TODO: make this work with more generic weight budget, i.e. not removing entire
//      edge, just some of the weight
void generate()
{
    int i, j, edges, temp, idx, idxw;
    int start;
    double rnum;
    
    //  Initialize variables within the bounds 
    for( i = 0; i < POPSIZE; i++ )
    {
        pop[i].fitness = 0;
        
        //TODO: for general weighted graphs, choose a budget B <= K
        //      while... choose edge, select R=rand(0,B)
        //      remove R from edge weight, repeat until full budget used up
        edges = i4_uniform_ab( 0, K );
        start = i4_uniform_ab( 0, N-1 );
        
        for( idx = start, temp = 0; idx < numEdges && temp < edges; idx++ )
        {
            rnum = r8_uniform_ab( 0.0, 1.0 );
            if( rnum < edge_prob[idx].prob )
            {
                pop[i].gene[temp].idx = idx;
                idxw = edge_prob[idx].idxu;
                pop[i].gene[temp].wght = ng->sg.wght[idxu];
                temp++;
            }
        }
        if( temp < edges && start != 0 )
        {
            for( idx = 0; idx < start && temp < edges; idx++ )
            {
                rnum = r8_uniform_ab( 0.0, 1.0 );
                if( rnum < edge_prob[idx].prob )
                {
                    pop[i].gene[temp].idx = idx;
                    idxw = edge_prob[idx].idxu;
                    pop[i].gene[temp].wght = ng->sg.wght[idxu];
                    temp++;
                }
            }
        }
        pop[i].len = temp;
    }
    
    return;
}

// MUTATE performs a random uniform mutation. 
//void mutate( int &seed )
void mutate()
{
    int i, inum;
    double rnum;
    
    for( i = 0; i < numEdges; i++ )
    {
        rnum = r8_uniform_ab( 0.0, 1.0 );
        if( rnum < PMUTATION )
        {
            inum = i4_uniform_ab( 0, 1 );
            edge_prob[i].prob = edge_prob[i].prob*(1.0-MUTSHIFT) +
                                ( inum ? MUTSHIFT : 0);
        }
    }
    
    return;
}
//****************************************************************************80
// END PBIL METHODS
//****************************************************************************80


//****************************************************************************80
//    MAIN supervises the genetic algorithm.
//****************************************************************************80
int main( int argc, char *argv[] )
{
    int generation;
    int i, c, idxB, idxW;
    char *filename = NULL;
    
    // Option handling 
    parse_arguments( &argc, &argv, options );
    if( argc < 1 ) die( "missing filename" );
    if( argc > 1 ) die( "trailing arguments" );
    filename = *argv;
    
    // Run the PBIL
    initialize( filename );
    for( generation = 0; generation < MAXGENS; generation++ )
    {
        if( generation%100 == 0 )
        {
            cout << " Generation " << generation << endl;
        }
        
        generate();
        evaluate( idxB, idxW ); //set idx's here
        if( pop[idxB].fitness < pop[POPSIZE].fitness )
        {
            pop[POPSIZE].fitness = pop[idxB].fitness;
            pop[POPSIZE].len = pop[idxB].len;
            for( i = 0; i < pop[idxB].len; i++ )
            {
                pop[POPSIZE].gene[i].idx = pop[idxB].gene[i].idx;
                pop[POPSIZE].gene[i].wght = pop[idxB].gene[i].wght;
            }
        }
        update( idxB, idxW ); //use idx's to update probs
        mutate();
    }
    
    // Print edge data
    cout << "\n";
    for( i = 0; i < numEdges; i++ )
    {
        if( edge_prob[i].prob > 0.05 )
        {
            cout << i << " " << edge_prob[i].prob << "\n";
        }
    }
    cout << "\n";
    
    // Print the best solution found
    cout << "\n";
    cout << "  Best member after " << MAXGENS << " generations:\n";
    cout << "\n";
    for( i = 0; i < pop[POPSIZE].len; i++ )
    {
        //cout << "  gene[" << i << "] = " << pop[POPSIZE].gene[i].idx << " ("
        //     << pop[POPSIZE].gene[i].idx << ", "
        //     << pop[POPSIZE].gene[i].idx << ")\n";
        cout << "  gene[" << i << "] = " << pop[POPSIZE].gene[i].idx << endl;
    }
    cout << "\n";
    cout << "  Best fitness = " << pop[POPSIZE].fitness << "\n";
    
    // Clean up and terminate
    saucy_free( s );
    ng->free( ng );
    free( marks );
    for( i = 0; i <= POPSIZE; i++ )
    {
        free( pop[i].gene );
    }
    free( pop );
    free( edge_prob );
    
    return 0;
}
