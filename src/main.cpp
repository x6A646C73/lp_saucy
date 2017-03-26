//Included libraries
//{{{
#include <cstdlib>
#include <cstdio>
#include <csignal>
#include <cstring>
#include <getopt.h>
#include "saucy.h"
#include "lp_amorph.h"
#include "warmup_theta.h"
#include "util.h"
#include "platform.h"
//}}}

// Static variables //
//{{{
static char *filename=NULL;   /* Graph file we're reading */
static char *genfile=NULL;   /* Generator file we're reading */
static int timeout = 0;      /* Seconds before quitting after refinement */
static sig_atomic_t timeout_flag = 0; /* Has the alarm gone off yet? */
static int stats_mode;  /* Print out stats when we're done */
static int quiet_mode;  /* Don't output automorphisms */
static int repeat = 1; /* Repeat count, for benchmarking */
static char *marks;    /* "Bit" vector for printing */
//static int wght_mode; /* use weighted graph */
//}}}

// Stats are global so we can print them from the signal handler //
struct saucy_stats stats;

struct option_desc {
    char *name;
    char letter;
    char *argname;
    char *description;
};

void print_options( const struct option_desc *opt )
{ //{{{
    int len, max;
    const struct option_desc *p;
    const char *s;

    max = 0;
    for( p = opt; p->name; ++p )
    { //{{{
        if( *p->description == '*' ) continue;
        len = strlen( p->name );
        if( p->argname ) len += strlen( p->argname ) + 1;
        if( max < len ) max = len;
    } //}}}

    for( p = opt; p->name; ++p )
    { //{{{
        if( p->letter ) printf( " -%c,", p->letter );
        else space( 4 );
                
        if( p->argname ) len = printf( " --%s=%s", p->name, p->argname );
        else len = printf( " --%s", p->name );
        
        s = p->description;
        if( *s == '*' )
        { //{{{
            space( 3 );
            puts( s+1 );
        } //}}}
        else
        { //{{{
            space( max - len + 6 );
            puts( s );
        } //}}}
    } //}}}
} //}}} END print_options

void parse_arguments( int argc, char **argv )
{ //{{{
    int c;
    static struct option_desc opts_desc[] = 
    { //{{{
        { "infile", 'f', 0,
            "file containing linear program as .lp or .mps" },
        { "genfile", 'g', 0,
            "file containing a head start on orbit partition" },
        //{ "outfile", 'o', 0,
        //    "file to store results of search" },
        { "quiet", 'q', 0,
            "do not print the generators found" },
        { "repeat", 'r', "N",
            "run saucy N times; used for benchmarking (implies -sq)" },
        { "stats", 's', 0,
            "output various statistics after the generators" },
        { "timeout", 't', "N",
            "after N seconds, the next generator will be the last" },
        { "help", 'h', 0,
            "output this help message" },
        { "version", 'v', 0,
            "version information" },
        { 0, 0, 0, 0 }
    }; //}}}
    static struct option long_options[] =
    { //{{{
            /* These options set a flag. */
            {"infile",  required_argument, 0, 'f'},
            {"genfile", required_argument, 0, 'g'},
            {"quiet",         no_argument, 0, 'q'},
            {"repeat",  required_argument, 0, 'r'},
            {"stats",         no_argument, 0, 's'},
            {"timeout", required_argument, 0, 't'},
            {"help",          no_argument, 0, 'h'},
            {"version",       no_argument, 0, 'v'},
            {0, 0, 0, 0}
    }; //}}}

    
    while( 1 )
    { //{{{
        // getopt_long stores the option index here. //
        int option_index = 0;
        
        c = getopt_long( argc, argv, "f:g:qr:st:hv",
                         long_options, &option_index );
        
        // Detect the end of the options. //
        if( c == -1 )
            break;
        
        switch( c )
        { //{{{
            case 0:
                // If this option set a flag, do nothing else now. //
                if( long_options[option_index].flag != 0 )
                    break;
                printf( "option %s", long_options[option_index].name );
                if( optarg )
                    printf( " with arg %s", optarg );
                printf( "\n" );
                break;
            case 'f':
                filename = (char *)malloc( strlen(optarg)*sizeof(char) );
                strcpy( filename, optarg );
                break;
            case 'g':
                genfile = (char *)malloc( strlen(optarg)*sizeof(char) );
                strcpy( genfile, optarg );
                break;
            case 'q':
                quiet_mode = 1;
                break;
            case 'r':
                repeat = atoi( optarg );
                if( repeat <= 0 ) die( "repeat count must be positive" );
                break;
            case 's':
                stats_mode = 1;
                break;
            case 't':
                timeout = atoi( optarg );
                if( timeout <= 0 ) die( "timeout must be positive" );
                break;
            case 'h':
                printf( "usage: saucy [OPTION]... FILE\n" );
                print_options( opts_desc );
                exit( 0 );
            case 'v':
                printf( "saucy %s\n", SAUCY_VERSION );
                exit( 0 );
            case '?':
                // getopt_long already printed an error message. //
                break;
            
            default:
                abort ();
        } //}}}
    } //}}}
    
    // Print any remaining command line arguments (not options). //
    if( optind < argc )
    { //{{{
        printf( "non-option ARGV-elements: " );
        while( optind < argc )
            printf( "%s ", argv[optind++] );
        putchar( '\n' );
    } //}}}
} //}}}

static void timeout_handler(void)
{ //{{{
    // Do nothing but set a flag to be tested during the search //
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
    struct saucy *s;
    struct lp_saucy *s_bak;
    struct saucy_stats stats_bak;
    struct lp_amorph_graph *g = NULL;
    long cpu_time;
    int i, num_gens;
    int temp, rep, repsize;
    int n, e, tw;
    FILE *f;
    //}}}
    
    // Option handling //
    parse_arguments( argc, argv );
    if( filename == NULL ) die( "no problem file provided" );
    
    // Repeating is for benchmarking //
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
    
    // Allocate some memory to facilitate printing //
    marks = (char *)calloc( n, sizeof( char ) );
    if( !marks ) die( "out of memory" );
    
    // Allocate saucy space //
    //NOTE: digraph_mode set to 0
    //s = saucy_alloc( n, tw, 0, &g->sg, on_automorphism, g, &stats );
    s = saucy_alloc( n, tw )
    if( s == NULL ) die( "saucy initialization failed" );
    
    s_bak = lp_saucy_alloc( n, tw, &stats_bak );
    if( s_bak == NULL ) die( "saucy backup initialization failed" );
    
    // Collect provided generators //
    if( genfile != NULL )
        num_gens = lp_warmup_theta( genfile, &s_bak, 0, &g->sg );
    
    // Set up the alarm for timeouts //
    if( timeout > 0 ) platform_set_timer( timeout, timeout_handler );
    
    // Print statistics when signaled //
    platform_set_user_signal( stats_handler );
    
    // Run the search //
    fflush( stdout );
    f = quiet_mode ? stderr : stdout;
    for( i = 0; i < repeat; ++i )
    { //{{{
        // Copy stats data //
        stats.grpsize_base = stats_bak.grpsize_base;
        stats.grpsize_exp = stats_bak.grpsize_exp;
        stats.levels = stats_bak.levels;
        stats.nodes = stats_bak.nodes;
        stats.bads = stats_bak.bads;
        stats.gens = stats_bak.gens;
        stats.support = stats_bak.support;
        
        cpu_time = platform_clock();
        
        saucy_search( &s, &g->sg, 0, g->colors, on_automorphism, g, &stats,
                      s_bak.thprev, s_bak.thnext, s_bak.thfront, s_bak.threp,
                      s_bak.thsize, s_bak.theta );
        
        fprintf( f, "%e\n", divide( platform_clock() - cpu_time,
                 PLATFORM_CLOCKS_PER_SEC ) );
    } //}}}
    
    // Copy stats data //
    stats.grpsize_base = stats_bak.grpsize_base;
    stats.grpsize_exp = stats_bak.grpsize_exp;
    stats.levels = stats_bak.levels;
    stats.nodes = stats_bak.nodes;
    stats.bads = stats_bak.bads;
    stats.gens = stats_bak.gens;
    stats.support = stats_bak.support;
    quiet_mode = 0;
    saucy_search( &s, &g->sg, 0, g->colors, on_automorphism, g, &stats,
                  s_bak.thprev, s_bak.thnext, s_bak.thfront, s_bak.threp,
                  s_bak.thsize, s_bak.theta );
    
    // Warn if timeout //
    if( timeout_flag ) warn( "search timed out" );
    
    // Print out stats if requested //
    if( stats_mode )
    { //{{{
        fprintf( f, "input file = %s\n", filename );
        if( g->stats ) g->stats( g, f );
        fprintf( f, "vertices = %d\n", n );
        fprintf( f, "edges = %d\n", g->sg.e );
        fprintf( f, "different weights = %d\n", tw );
        print_stats( f );
    } //}}}
    
    // Cleanup //
    saucy_free( s );
    lp_saucy_free( s_bak );
    g->free( g );
    free( marks );
    
    return EXIT_SUCCESS;
} //}}} END main
