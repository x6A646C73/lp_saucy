#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <fstream>
#include <map>
#include <vector>
#include "saucy.h"
#include "util.h"
#include "amorph.h"
#include "platform.h"
#include "CoinMpsIO.hpp"


//TODO: move reading functions out of here...
//      create lp_saucyio.{h,cpp}
void amorph_graph_free( struct amorph_graph *g )
{
    free(g->sg.adj);
    free(g->sg.edg);
    free(g->colors);
    free(g);
}

static int init_fixadj1( int n, int *adj )
{
    int val, sum, i;

    /* Translate adj values to real locations */
    val = adj[0]; sum = 0; adj[0] = 0;
    for (i = 1; i < n; ++i) {
        sum += val;
        val = adj[i];
        adj[i] = sum;
    }
    return sum + val;
}

static void init_fixadj2( int n, int e, int *adj )
{
    int i;

    /* Translate again-broken sizes to adj values */
    for (i = n-1; i > 0; --i) {
        adj[i] = adj[i-1];
    }
    adj[0] = 0;
    adj[n] = e;
}

static void amorph_print_automorphism(
    int n, const int *gamma, int nsupp, const int *support,
    struct amorph_graph *g, char *marks )
{
    int i, j, k;

    /* We presume support is already sorted */
    for (i = 0; i < nsupp; ++i) {
        k = support[i];

        /* Skip elements already seen */
        if (marks[k]) continue;

        /* Start an orbit */
        marks[k] = 1;
        printf( "(%s", g->var_names[k] );
        //printf("(%d", k);

        /* Mark and notify elements in this orbit */
        for (j = gamma[k]; j != k; j = gamma[j]) {
            marks[j] = 1;
            printf( " %s", g->var_names[j] );
            //printf(" %d", j);
        }

        /* Finish off the orbit */
        putchar(')');
    }
    putchar('\n');

    /* Clean up after ourselves */
    for (i = 0; i < nsupp; ++i) {
        marks[support[i]] = 0;
    }
}

/* return value >1 indicates error */
static int dupe_check( int n, int *adj, int *edg )
{
    int i, j, self_loop_ctr;
    int *dupe_tmp = (int *)calloc(n, sizeof(int));
    if (!dupe_tmp) {
        warn("can't allocate memory");
        free(dupe_tmp);
        return 2;
    }

    /* Check outgoing edges of each vertex for duplicate endpoints */
    for (i = 0; i < n; ++i) {
        self_loop_ctr = 0;
        for (j = adj[i] ; j < adj[i+1] ; j++) {
            /* Self-loops lead to two entries of edg[j]==i,
             * so we check for those and only worry if we see
             * 3 hits of edg[j]==i (which means 2 self-loops).
             */
            if (edg[j] == i) {
                ++self_loop_ctr;
                if (self_loop_ctr > 2) {
                    warn("duplicate edge in input");
                    free(dupe_tmp);
                    return 1;
                }
            }
            /* If we have recorded this vertex as connected to i,
             * we have a dupe.
             * Using i+1 because we used calloc above, and 0 is
             * a valid vertex index.
             */
            else if (dupe_tmp[edg[j]] == i+1) {
                warn("duplicate edge in input");
                free(dupe_tmp);
                return 1;
            }
            dupe_tmp[edg[j]] = i+1;
        }
    }

    free(dupe_tmp);
    return 0;
}

//TODO: duplicate for lp files or regular graph files
struct amorph_graph* amorph_read( const char *filename )
{
    int i, j, vars, cons;
    const double *obj, *rhs, *val;
    const int *row_ind, *col_ind;
    const CoinPackedMatrix *mat;
    CoinMpsIO prob;
    int n, e, *aout, *ain, *eout, *ein, *colors;
    int tempk, tempj;
    int w, *wout, *win; /* weight data */
    int ndx;
    struct amorph_graph *g = NULL;
    char **var_names;
    
    //TODO: use separate functions for LP and MPS, but use memcpy to set vals,
    //      including names?
    //      Test this in LP2Graph_full.cpp(?)
    //TODO: need to add back in ain, ein, etc
    prob.readMps( filename, "mps" );
    obj = prob.getObjCoefficients(); //gives vars their color
    rhs = prob.getRightHandSide();   //holds color for other side of bipartite
    vars = prob.getNumCols();        //number of variables
    cons = prob.getNumRows();        //number of constraints
    mat = prob.getMatrixByCol();     //contains weights for edges
    
    e = mat->getNumElements();
    val = mat->getElements();
    row_ind = mat->getIndices();
    col_ind = mat->getVectorStarts();
    
    std::map<double,int> varColors;
    std::map<double,int> conColors;
    std::pair<std::map<double,int>::iterator,bool> ret;
    int c1 = 0;
    for( i = 0; i < vars; i++ )
    {
        ret = varColors.insert( std::pair<double,int>( obj[i], c1 ) );
        if( ret.second == true ) c1++;
    }
    
    int c2 = 0;
    for( i = 0; i < cons; i++ )
    {
        //ret = conColors.insert( std::pair<double,int>( rhs[i], c1 + c2 + 1 ) );
        ret = conColors.insert( std::pair<double,int>( rhs[i], c1 + c2 ) );
        if( ret.second == true ) c2++;
    }
    
    std::map<double,int> coColors;
    int c3 = 0;
    for( i = 0; i < e; i++ )
    {
        //ret = coColors.insert( std::pair<double,int>( elm[i], c1+c2+c3+1 ) );
        ret = coColors.insert( std::pair<double,int>( val[i], c1+c2+c3 ) );
        if( ret.second == true ) c3++;
    }
    
    n = vars+cons; 
    var_names = (char **)malloc( n*sizeof(char*) );
    for( i = 0; i < vars; i++ )
    {
        tempk = strlen( prob.columnName(i) );
        var_names[i] = (char *)malloc( tempk*sizeof(char) );
        strcpy( var_names[i], prob.columnName(i) );
    }
    for( i = 0; i < cons; i++ )
    {
        tempk = strlen( prob.rowName(i) );
        var_names[i+vars] = (char *)malloc( tempk*sizeof(char) );
        strcpy( var_names[i+vars], prob.rowName(i) );
    }
    
    /* Allocate everything */
    //TODO: eventually allocate smarter, only allocate when needed, delete
    //      objects as they are no longer needed
    g = (struct amorph_graph*)malloc( sizeof(struct amorph_graph) );
    //aout = (int *)calloc( digraph ? (2*n+2) : (n+1), sizeof(int) );
    aout = (int *)calloc( (n+1), sizeof(int) );
    eout = (int *)malloc( 2 * e * sizeof(int) );
    wout = (int *)malloc( 2 * e * sizeof(int) ); /* weight data */
    colors = (int *)malloc( n * sizeof(int) );
    if( !g || !aout || !eout || !colors ) goto out_free;
    
    g->sg.n = n;
    g->sg.e = e;
    g->sg.w = coColors.size();
    g->sg.adj = aout;
    g->sg.edg = eout;
    g->sg.wght = wout; /* weight data */
    g->colors = colors;
    g->var_names = var_names;
    
    ain = aout;
    ein = eout;
    win = wout;
    
    for( i = 0; i < vars; i++ )
        colors[i] = varColors.find(obj[i])->second;
    for( i = 0; i < cons; i++ )
        colors[vars+i] = conColors.find(rhs[i])->second;
    
    /* Count the size of each adjacency list */
    //j gives column (variable), i should give row (constraint)
    for( j = 0; j < vars; j++ )
    {
        for( i = col_ind[j]; i < col_ind[j+1]; i++ )
        {
            ++ain[row_ind[i]+vars]; ++aout[j];
            //std::cout << j << " " << row_ind[i]+vars << "\n";
        }
    }
    
    /* Fix that */
    init_fixadj1( n, aout );
    //if( digraph ) init_fixadj1( n, ain );
    
    /* Insert adjacencies */
    for( j = 0; j < vars; j++ )
    {
        for( i = col_ind[j]; i < col_ind[j+1]; i++ )
        {
            tempk = ain[row_ind[i]+vars]++; tempj = aout[j]++;
            eout[tempj] = row_ind[i]+vars; ein[tempk] = j;
            w = coColors.find(val[i])->second;
            wout[tempj] = w; win[tempk] = w;
            //std::cout << j << " " << row_ind[i]+vars << " " << w << "\n";
        }
    }

    /* Fix that too */
    //if( digraph )
    //{
    //    init_fixadj2( n, e, aout );
    //    init_fixadj2( n, e, ain );
    //}
    //else
        init_fixadj2( n, 2 * e, aout );
    
    /* Check for duplicate edges */
    if( dupe_check( n, aout, eout ) ) goto out_free;
    
    /* Assign the functions */
    /* TODO: determine if these functions need to be altered at all */
    g->consumer = amorph_print_automorphism;
    g->free = amorph_graph_free;
    g->stats = NULL;
    goto out;
    
out_free:
    free( g );
    free( aout );
    free( eout );
    free( wout );
    free( colors );
    g = NULL;
out:
    return g;
}

static int timeout = 0;      /* Seconds before quitting after refinement */
static sig_atomic_t timeout_flag = 0; /* Has the alarm gone off yet? */
static int stats_mode;  /* Print out stats when we're done */
static int quiet_mode;  /* Don't output automorphisms */
static int repeat = 0; /* Repeat count, for benchmarking */
static int first;      /* Have we seen the first automorphism? (for gap) */
static char *marks;    /* "Bit" vector for printing */

/* Stats are global so we can print them from the signal handler */
struct saucy_stats stats;

static void arg_stats( char *arg ) { stats_mode = 1; }
static void arg_quiet( char *arg ) { quiet_mode = 1; }

static void arg_timeout( char *arg )
{
    timeout = atoi( arg );
    if( timeout <= 0 ) die( "timeout must be positive" );
}

static void arg_repeat( char *arg )
{
    repeat = atoi( arg );
    if( repeat <= 0 ) die( "repeat count must be positive" );
}

static void arg_version( char *arg )
{
    printf( "saucy %s\n", SAUCY_VERSION );
    exit( 0 );
}

static void arg_help( char *arg );

//TODO: file options, generators and graph...
static struct option options[] = {
    { "stats", 's', 0, arg_stats,
    "output various statistics after the generators" },
    { "quiet", 'q', 0, arg_quiet,
    "do not print the generators found" },
    { "timeout", 't', "N", arg_timeout,
    "after N seconds, the next generator will be the last" },
    { "repeat", 'r', "N", arg_repeat,
    "run saucy N times; used for benchmarking (implies -sq)" },
    { "help", 0, 0, arg_help,
    "output this help message" },
    { "version", 0, 0, arg_version,
    "version information" },
    { 0, 0, 0, 0, 0 }
};

static void arg_help( char *arg )
{
    printf( "usage: saucy [OPTION]... FILE\n" );
    print_options( options );
    exit( 0 );
}

static void timeout_handler( void )
{
    /* Do nothing but set a flag to be tested during the search */
    timeout_flag = 1;
}

static void print_stats( FILE *f )
{
    fprintf( f, "group size = %fe%d\n",
             stats.grpsize_base, stats.grpsize_exp );
    fprintf( f, "levels = %d\n", stats.levels );
    fprintf( f, "nodes = %d\n", stats.nodes );
    fprintf( f, "generators = %d\n", stats.gens );
    fprintf( f, "total support = %d\n", stats.support );
    fprintf( f, "average support = %.2f\n",
             divide( stats.support, stats.gens ) );
    fprintf( f, "nodes per generator = %.2f\n",
             divide( stats.nodes, stats.gens ) );
    fprintf( f, "bad nodes = %d\n", stats.bads );
}

static void stats_handler( void )
{
    fprintf( stderr, "========= intermediate stats ===========\n" );
    print_stats( stderr );
    fprintf( stderr, "========================================\n" );
}

static int on_automorphism( int n, const int *gamma, int k, int *support, void *arg )
{
    struct amorph_graph *g = (struct amorph_graph *)arg;
    if( !quiet_mode )
    {
        qsort_integers( support, k );
        //if( gap_mode )
        //{
        //    putchar( !first ? '[' : ',' );
        //    putchar( '\n' );
        //    first = 1;
        //}
        g->consumer( n, gamma, k, support, g, marks );
    }
    return !timeout_flag;
}

int main( int argc, char **argv )
{
    struct saucy *s;
    struct amorph_graph *g = NULL;
    long cpu_time;
    int i, n, e;
    int tw;
    FILE *f;
    char *filename;   /* Graph file we're reading */
    
    /* Option handling */
    parse_arguments( &argc, &argv, options );
    if( argc < 1 ) die( "missing filename" );
    if( argc > 1 ) die( "trailing arguments" );
    filename = *argv;
    
    /* Repeating is for benchmarking */
    if( repeat > 1 ) quiet_mode = stats_mode = 1;
    
    /* Read the input file */
    g = amorph_read( filename );
    if( !g ) die( "unable to read input file" );
    n = g->sg.n;
    e = g->sg.e;
    tw = g->sg.w;
    
    /* Allocate some memory to facilitate printing */
    marks = (char *)calloc( n, sizeof(char) );
    if( !marks ) die( "out of memory" );
    
    /* Allocate saucy space */
    s = (struct saucy*)saucy_alloc( n, tw );
    if( s == NULL ) die( "saucy initialization failed" );
    
    /* Set up the alarm for timeouts */
    if( timeout > 0 ) platform_set_timer( timeout, timeout_handler );
    
    /* Print statistics when signaled */
    platform_set_user_signal( stats_handler );
    
    /* Start timing */
    /* cpu_time = platform_clock(); */
    
    /* Run the search */
    fflush( stdout );
    f = quiet_mode ? stdout : stderr;
    for( i = 0; i < repeat; ++i )
    {
        cpu_time = platform_clock();
        
        saucy_search( s, &g->sg, 0, g->colors, on_automorphism,
                      g, &stats );
        
        fprintf( f, "%e\n", divide( platform_clock() - cpu_time,
                                    PLATFORM_CLOCKS_PER_SEC ) );
    }
    
    quiet_mode = 0;
    saucy_search( s, &g->sg, 0, g->colors, on_automorphism,
                  g, &stats );
    
    /* Finish timing */
    /* cpu_time = platform_clock() - cpu_time; */
    
    /* Warn if timeout */
    if( timeout_flag ) warn( "search timed out" );
    
    /* Print out stats if requested */
    if( stats_mode )
    {
        fprintf( f, "input file = %s\n", filename);
        if( g->stats ) g->stats( g, f );
        fprintf( f, "vertices = %d\n", n );
        fprintf( f, "edges = %d\n", g->sg.e );
        fprintf( f, "different weights = %d\n", tw );
        print_stats( f );
    }
    
    /* Cleanup */
    saucy_free( s );
    g->free( g );
    free( marks );
    
    /* That's it, have a nice day */
    return EXIT_SUCCESS;
}
