#include <cstdlib>
#include <cstdio>
#include <map>
#include <zlib.h>
#include "saucy.h"
#include "util.h"
#include "amorph.h"
#include "CoinMpsIO.hpp"


void amorph_graph_free( struct amorph_graph *g )
{
    free(g->sg.adj);
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

static void amorph_print_automorphism(
    int n, const int *gamma, int nsupp, const int *support,
    struct amorph_graph *g, char *marks )
{
    int i, j, k;
    
    /* We presume support is already sorted */
    for( i = 0; i < nsupp; ++i )
    {
        k = support[i];
        
        /* Skip elements already seen */
        if( marks[k] ) continue;
        
        /* Start an orbit */
        marks[k] = 1;
        if( g->var_names )
            printf( "(%s", g->var_names[k] );
        else
            printf( "(%d", k );
        
        /* Mark and notify elements in this orbit */
        for( j = gamma[k]; j != k; j = gamma[j] )
        {
            marks[j] = 1;
            if( g->var_names )
                printf( " %s", g->var_names[j] );
            else
                printf( " %d", j );
        }
        
        /* Finish off the orbit */
        putchar( ')' );
    }
    putchar( '\n' );
    
    /* Clean up after ourselves */
    for( i = 0; i < nsupp; ++i )
    {
        marks[support[i]] = 0;
    }
}

struct amorph_graph* read_lp( const char *filename )
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
        }
    }
    
    /* Fix that */
    init_fixadj1( n, aout );
    
    /* Insert adjacencies */
    for( j = 0; j < vars; j++ )
    {
        for( i = col_ind[j]; i < col_ind[j+1]; i++ )
        {
            tempk = ain[row_ind[i]+vars]++; tempj = aout[j]++;
            eout[tempj] = row_ind[i]+vars; ein[tempk] = j;
            w = coColors.find(val[i])->second;
            wout[tempj] = w; win[tempk] = w;
        }
    }

    /* Fix that too */
    init_fixadj2( n, 2 * e, aout );
    
    /* Check for duplicate edges */
    if( dupe_check( n, aout, eout ) ) goto out_free;
    
    /* Assign the functions */
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

static int read_int( gzFile f, int *k )
{
    int c, r = 0, neg = 0;
    for (c = gzgetc(f); c != '-' && !isdigit(c); ) {
        for (; isspace(c); c = gzgetc(f));
        for (; c == 'c'; c = gzgetc(f)) {
            while ((c = gzgetc(f)) != '\n') {
                if (c == EOF) return 0;
            }
        }
    }
    if (c == '-') {
        neg = 1;
        c = gzgetc(f);
    }
    if (!isdigit(c)) return 0;
    for (; isdigit(c); c = gzgetc(f)) {
        r = r * 10 + c - '0';
    }
    *k = neg ? -r : r;
    return isspace(c);
}

struct amorph_graph* read_graph( const char *filename )
{
    int i, j, k, n, e, *aout, *eout, *ain, *ein, *colors;
    int tempk, tempj;
    int w, *wout, *win; /* weight data */
    int ndx;
    struct amorph_graph *g = NULL;
    gzFile f;
    
    /* Open the file */
    f = gzopen(filename, "r");
    if( f == NULL ) goto out;
    
    /* Read the sizes */
    if( !read_int(f, &n) || !read_int(f, &e) || !read_int(f, &w) )
    {
        goto out_close;
    }
    
    /* Allocate everything */
    g = (struct amorph_graph *)malloc( sizeof(struct amorph_graph) );
    aout = (int *)calloc( n+1, sizeof(int));
    eout = (int *)malloc(2 * e * sizeof(int));
    wout = (int *)malloc(2 * e * sizeof(int)); /* weight data */
    //wlist = (int *)calloc( 2 * e, sizeof(int) );
    colors = (int *)malloc(n * sizeof(int));
    if( !g || !aout || !eout || !colors ) goto out_free;
    
    g->sg.n = n;
    g->sg.e = e;
    g->sg.w = w;
    g->sg.adj = aout;
    g->sg.edg = eout;
    g->sg.wght = wout; /* weight data */
    g->colors = colors;
    
    ain = aout;
    ein = eout;
    win = wout; /* weight data */
    
    /* Initial coloring with provided splits */
    for( i = 0; i < n; i++)
    {
        if( !read_int(f, &j) || !read_int(f, &k) ) goto out_free;
        colors[j] = k;
    }
    
    /* Count the size of each adjacency list */
    for( i = 0; i < e; ++i )
    {
        if( !read_int(f, &j) || !read_int(f, &k) || !read_int(f, &w) ) goto out_free;
        ++aout[j]; ++ain[k];
    }
    
    /* Fix that */
    init_fixadj1(n, aout);
    
    /* Go back to the front of the edges */
    gzrewind(f);
    for( i = 0; i < 3 + 2*n; ++i )
    {
        if( !read_int(f, &k) ) goto out_free;
    }
    
    /* Insert adjacencies */
    for( i = 0; i < e; ++i )
    {
        if( !read_int(f, &j) || !read_int(f, &k) || !read_int(f, &w) ) goto out_free;
        
        /* Simple input validation: check vertex values */
        if( j >= n || j < 0)
        {
            warn( "invalid vertex in input: %d", j );
            goto out_free;
        }
        if( k >= n || k < 0 )
        {
            warn( "invalid vertex in input: %d", k );
            goto out_free;
        }
        
        tempj = aout[j]++; tempk = ain[k]++;
        eout[tempj] = k; ein[tempk] = j;
        
        /* Collect weight data */
        wout[tempj] = w; win[tempk] = w;
    }
    
    /* Fix that too */
    init_fixadj2(n, 2 * e, aout);
    
    /* Check for duplicate edges */
    if( dupe_check( n, aout, eout ) ) goto out_free;
    
    /* Assign the functions */
    g->consumer = amorph_print_automorphism;
    g->free = amorph_graph_free;
    g->stats = NULL;
    goto out_close;
    
out_free:
    free(g);
    free(aout);
    free(eout);
    free(wout);
    free(colors);
    g = NULL;
out_close:
    gzclose(f);
out:
    return g;
}
