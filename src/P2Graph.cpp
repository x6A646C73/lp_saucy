//Included libraries
//{{{
#include <cstdlib>
#include <cstdio>
#include <csignal>
#include <cstring>
#include "saucy.h"
#include "lp_amorph.h"
#include "util.h"
#include "platform.h"
//}}}

/* Static variables */
//{{{
static char *filename=NULL;   /* Graph file we're reading */
static char *genfile=NULL;   /* Generator file we're reading */
static char *outfile=NULL;   /* Output file we're using */
static int timeout = 0;      /* Seconds before quitting after refinement */
static sig_atomic_t timeout_flag = 0; /* Has the alarm gone off yet? */
static int stats_mode;  /* Print out stats when we're done */
static int quiet_mode;  /* Don't output automorphisms */
static int repeat = 1; /* Repeat count, for benchmarking */
static char *marks;    /* "Bit" vector for printing */
//static int wght_mode; /* use weighted graph */
//}}}

/* Stats are global so we can print them from the signal handler */
struct saucy_stats stats;

static void arg_stats( char *arg ) { stats_mode = 1; }
static void arg_quiet( char *arg ) { quiet_mode = 1; }
//static void arg_wght( char *arg ) { wght_mode = 1; }

static void arg_timeout( char *arg )
{ //{{{
    timeout = atoi( arg );
    if( timeout <= 0 ) die( "timeout must be positive" );
} //}}} END arg_timeout

static void arg_repeat(char *arg)
{ //{{{
    repeat = atoi( arg );
    if( repeat <= 0 ) die( "repeat count must be positive" );
} //}}} END arg_repeat

static void arg_version( char *arg )
{ //{{{
    printf( "saucy %s\n", SAUCY_VERSION );
    exit( 0 );
} //}}} END arg_version

static void arg_infile( char *arg )
{ //{{{
    filename = arg;
    if( !filename ) die( "failed to retrieve input file" );
} //}}}

static void arg_outfile( char *arg )
{ //{{{
    outfile = arg;
    if( !outfile ) die( "failed to retrieve output file" );
} //}}}

static void arg_genfile( char *arg )
{ //{{{
    genfile = arg;
    if( !genfile ) die( "failed to retrieve generator file" );
} //}}}

static void arg_help( char *arg );

static struct option options[] = 
{ //{{{
    { "infile", 'f', 0, arg_infile,
      "file containing linear program as .lp or .mps" },
    { "genfile", 'g', 0, arg_genfile,
      "file containing a head start on orbit partition" },
    { "outfile", 'o', 0, arg_outfile,
      "file to store results of search" },
    { "quiet", 'q', 0, arg_quiet,
      "do not print the generators found" },
    { "repeat", 'r', "N", arg_repeat,
      "run saucy N times; used for benchmarking (implies -sq)" },
    { "stats", 's', 0, arg_stats,
      "output various statistics after the generators" },
    { "timeout", 't', "N", arg_timeout,
      "after N seconds, the next generator will be the last" },
    //{ "weight", 'w', 0, arg_wght,
    //  "use weighted graph version of linear program" },
    { "help", 0, 0, arg_help,
      "output this help message" },
    { "version", 0, 0, arg_version,
      "version information" },
    { 0, 0, 0, 0, 0 }
}; //}}}

static void arg_help( char *arg )
{ //{{{
    printf( "usage: saucy [OPTION]... FILE\n" );
    print_options( options );
    exit( 0 );
} //}}} END arg_help

static void timeout_handler(void)
{ //{{{
    /* Do nothing but set a flag to be tested during the search */
    timeout_flag = 1;
} //}}} END timeout_handler

static void print_stats( FILE *f )
{ //{{{
    fprintf(f, "group size = %fe%d\n",
        stats.grpsize_base, stats.grpsize_exp);
    fprintf(f, "levels = %d\n", stats.levels);
    fprintf(f, "nodes = %d\n", stats.nodes);
    fprintf(f, "generators = %d\n", stats.gens);
    fprintf(f, "total support = %d\n", stats.support);
    fprintf(f, "average support = %.2f\n",
        divide(stats.support, stats.gens));
    fprintf(f, "nodes per generator = %.2f\n",
        divide(stats.nodes, stats.gens));
    fprintf(f, "bad nodes = %d\n", stats.bads);
} //}}} END print_stats

static void stats_handler(void)
{ //{{{
    fprintf(stderr, "========= intermediate stats ===========\n");
    print_stats(stderr);
    fprintf(stderr, "========================================\n");
} //}}} END stats_handler

static int on_automorphism(int n, const int *gamma, int k, int *support, void *arg)
{ //{{{
    struct lp_amorph_graph *g = (struct lp_amorph_graph *)arg;
    
    if( !quiet_mode ) 
    { //{{{
        qsort_integers(support, k);
        //if( gap_mode ) 
        //{ //{{{
        //    putchar(!first ? '[' : ',');
        //    putchar('\n');
        //    first = 1;
        //} //}}}
        g->consumer(n, gamma, k, support, g, marks);
    } //}}}
    
    return !timeout_flag;
} //}}} END on_automorphism

int main( int argc, char **argv )
{ //{{{
    // Declare variables
    //{{{
    struct saucy *s, *s_bak;
    struct lp_amorph_graph *g = NULL;
    long cpu_time;
    int i, num_gens;
    int temp, rep, repsize;
    int n, e, tw;
    FILE *f;
    //}}}
    
    /* Option handling */
    //TODO: make sure this fixes this since files have options now
    parse_arguments( &argc, &argv, options );
    if( argc > 1 ) die( "trailing arguments" );
    
    /* Repeating is for benchmarking */
    if (repeat > 1) quiet_mode = stats_mode = 1;
    
    // Produce Weighted Graph from MPS of LP
    g = lp_amorph_read_build( filename );
    //if( wght_mode ) g = lp_amorph_read_build( filename );
    // Produce Unweighted Graph from MPS
    //else  g = lp_amorph_read_build( filename ); //not currently implemented
    if( !g ) die( "unable to build graph" );
    n = g->sg.n;
    e = g->sg.e;
    tw = g->sg.w;
    
    /* Allocate some memory to facilitate printing */
    marks = (char *)calloc( n, sizeof( char ) );
    if( !marks ) die( "out of memory" );
    
    /* Allocate saucy space */
    //NOTE: digraph_mode set to 0
    s = saucy_alloc( n, tw, 0, &g->sg, on_automorphism, g, &stats );
    if( s == NULL ) die( "saucy initialization failed" );
    
    s_bak = saucy_alloc( n, tw, 0, &g->sg, on_automorphism, g, &stats );
    if( s_bak == NULL ) die( "saucy backup initialization failed" );
    
    /* Collect provided generators */
    if( genfile != NULL )
    { //{{{
        num_gens = lp_warmup_theta( genfile, s );
        memcpy( s_bak, s, sizeof( struct saucy ) ); //TODO: versify this is correct
    } //}}}
    
    /* Set up the alarm for timeouts */
    if( timeout > 0 ) platform_set_timer( timeout, timeout_handler );
    
    /* Print statistics when signaled */
    platform_set_user_signal( stats_handler );
    
    //TODO: fix output stuff
    /* Run the search */
    if( outfile == NULL ) f = (quiet_mode) ? stderr : stdout;
    else f = fopen( outfile, "r" );
    fflush( f );
    //f = quiet_mode ? stdout : stderr;
    for( i = 0; i < repeat; ++i )
    { //{{{
        cpu_time = platform_clock();

        saucy_search( s, g->colors );
        
        fprintf( f, "%e\n", divide( platform_clock() - cpu_time,
                 PLATFORM_CLOCKS_PER_SEC ) );
        
        memcpy( s, s_bak, sizeof( struct saucy ) ); //TODO: verify this is correct
    } //}}}
    
    quiet_mode = 0;
    saucy_search( s, g->colors );
    
    /* Finish timing */
    /* cpu_time = platform_clock() - cpu_time; */
    
    //if( gap_mode && !quiet_mode ) printf( "\n]\n" );
    
    /* Warn if timeout */
    if( timeout_flag ) warn( "search timed out" );
    
    /* Print out stats if requested */
    if( stats_mode )
    { //{{{
        // f is the output file?
        fprintf( f, "input file = %s\n", filename );
        if( g->stats ) g->stats( g, f );
        fprintf( f, "vertices = %d\n", n );
        fprintf( f, "edges = %d\n", g->sg.e );
        fprintf( f, "different weights = %d\n", tw );
        print_stats( f );
        /*
        fprintf(f, "cpu time (s) = %e\n",
            divide(cpu_time, PLATFORM_CLOCKS_PER_SEC));
        */
    } //}}}
    
    /* Cleanup */
    saucy_free( s );
    saucy_free( s_bak );
    g->free( g );
    free( marks );
    
    return EXIT_SUCCESS;
} //}}} END main
