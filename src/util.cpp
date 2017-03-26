#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>

#include "util.h"

static const char *progname;

static void vwarn( const char *fmt, va_list args )
{ //{{{
    if (progname) fprintf( stderr, "%s: ", progname );
    vfprintf( stderr, fmt, args );
    putc( '\n', stderr );
} //}}} END vwarn

void warn( const char *fmt, ... )
{ //{{{
    va_list args;
    va_start( args, fmt );
    vwarn( fmt, args );
    va_end( args );
} //}}} END warn

void die( const char *fmt, ... )
{ //{{{
    va_list args;
    va_start( args, fmt );
    vwarn( fmt, args );
    va_end( args );
    exit( EXIT_FAILURE );
} //}}} END die

void bang( const char *fmt, ... )
{ //{{{
    const char *err = strerror( errno );
    va_list args;
    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    va_end( args );
    fprintf( stderr, ": %s\n", err );
    exit( EXIT_FAILURE );
} //}}} END bang

void space( int k )
{ //{{{
    while( k-- ) putchar( ' ' );
} //}}} END space

/*
void print_options( const struct option_desc *opt )
{ //{{{
    int len, max;
    const struct option *p;
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

void parse_options( int argc, char **argv )
{ //{{{
    int c;
    
    while( 1 )
    { //{{{
        //TODO: move this outside the loop?
        static struct option long_options[] =
        { //{{{
            // These options set a flag. 
            {"infile",  required_argument, 0, 'f'},
            {"genfile", required_argument, 0, 'g'},
            {"outfile", required_argument, 0, 'o'},
            {"quiet",         no_argument, 0, 'q'},
            {"repeat",  required_argument, 0, 'r'},
            {"stats",         no_argument, 0, 's'},
            {"timeout", required_argument, 0, 't'},
            {"help",          no_argument, 0, 'h'},
            {"version",       no_argument, 0, 'v'},
            {0, 0, 0, 0}
        }; //}}}
        // getopt_long stores the option index here. 
        int option_index = 0;
        
        c = getopt_long( argc, argv, "f:g:o:qr:st:hv",
                         long_options, &option_index );
        
        // Detect the end of the options. 
        if( c == -1 )
            break;
        
        switch( c )
        { //{{{
            case 0:
                // If this option set a flag, do nothing else now. 
                if( long_options[option_index].flag != 0 )
                    break;
                printf( "option %s", long_options[option_index].name );
                if( optarg )
                    printf( " with arg %s", optarg );
                printf( "\n" );
                break;
            case 'f':
                strcpy( filename, optarg );
                break;
            case 'g':
                strcpy( genfile, optarg );
                break;
            case 'o':
                strcpy( outfile, optarg );
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
                timeout = atoi( arg );
                if( timeout <= 0 ) die( "timeout must be positive" );
                break;
            case 'h':
                printf( "usage: saucy [OPTION]... FILE\n" );
                print_options( options );
                exit( 0 );
            case 'v':
                printf( "saucy %s\n", SAUCY_VERSION );
                exit( 0 );
            case '?':
                // getopt_long already printed an error message. 
                break;
            
            default:
                abort ();
        } //}}}
    } //}}}
    
    // Print any remaining command line arguments (not options). 
    if( optind < argc )
    { //{{{
        printf( "non-option ARGV-elements: " );
        while( optind < argc )
            printf( "%s ", argv[optind++] );
        putchar( '\n' );
    } //}}}
} //}}}
*/
