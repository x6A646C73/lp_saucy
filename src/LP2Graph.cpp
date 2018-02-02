#include <cstdio>
#include <cstdlib>
#include <csignal>
#include "saucy.h"
#include "util.h"
#include "amorph.h"
#include "platform.h"
#include "CoinMpsIO.hpp"

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
    g = read_lp( filename );
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
