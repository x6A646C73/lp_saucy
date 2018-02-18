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

//TODO: make work with .lp and .mps
struct amorph_graph* read_lp( const char *filename )
{
    int i, j, vars, cons;
    const double *obj, *rhs, *val;
    const int *row_ind, *col_ind;
    const CoinPackedMatrix *mat;
    CoinMpsIO prob;
    int n, e, w, *aout, *colors;
    int tempk;
    struct amorph_graph *g = NULL;
    char **var_names;
    
    //TODO: use separate functions for LP and MPS, but use memcpy to set vals,
    //      including names?
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
    g = (struct amorph_graph*)malloc( sizeof(struct amorph_graph) );
    aout = (int *)calloc( n*n, sizeof(int) );
    colors = (int *)malloc( n * sizeof(int) );
    if( !g || !aout || !colors ) goto out_free;
    
    g->sg.n = n;
    g->sg.e = e;
    g->sg.w = coColors.size();
    g->sg.adj = aout;
    g->colors = colors;
    g->var_names = var_names;
    
    for( i = 0; i < vars; i++ )
        colors[i] = varColors.find(obj[i])->second;
    for( i = 0; i < cons; i++ )
        colors[vars+i] = conColors.find(rhs[i])->second;
    
    /* Insert adjacencies */
    for( j = 0; j < vars; j++ )
    {
        for( i = col_ind[j]; i < col_ind[j+1]; i++ )
        {
            w = coColors.find(val[i])->second;
            aout[ j*n + row_ind[i]+vars ] = w;
            aout[ j + (row_ind[i]+vars)*n ] = w;
        }
    }
    
    /* Assign the functions */
    g->consumer = amorph_print_automorphism;
    g->free = amorph_graph_free;
    g->stats = NULL;
    goto out;
    
out_free:
    free( g );
    free( aout );
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
    int i, j, k, n, e, *aout, *colors;
    int w;
    struct amorph_graph *g = NULL;
    gzFile f;
    
    /* Open the file */
    f = gzopen( filename, "r" );
    if (f == NULL) goto out;
    
    /* Read the sizes */
    //NOTE: this graph format starts with #vertices #edges #edge weights 
    if( !read_int(f, &n) || !read_int(f, &e) || !read_int(f, &w) )
    {
        goto out_close;
    }
    
    /* Allocate everything */
    g = (struct amorph_graph*)malloc( sizeof(struct amorph_graph) );
    aout = (int *)calloc( n*n, sizeof(int) );
    colors = (int *)malloc( n * sizeof(int) );
    if (!g || !aout || !colors ) goto out_free;
    
    g->sg.n = n;
    g->sg.e = e;
    g->sg.w = w;
    g->sg.adj = aout;
    g->colors = colors;
    g->var_names = NULL;
    
    /* Initial coloring with provided splits */
    for( i = 0; i < n; i++ )
    {
        if( !read_int(f, &j) || !read_int(f, &k) ) goto out_free;
        colors[j] = k;
    }
    
    /* Count the size of each adjacency list */
    for( i = 0; i < e; ++i )
    {
        if( !read_int(f, &j) || !read_int(f, &k) || !read_int(f, &w) ) goto out_free;
        aout[i*n + j] = w;
        aout[j*n + i] = w;
    }
    
    /* Assign the functions */
    g->consumer = amorph_print_automorphism;
    g->free = amorph_graph_free;
    g->stats = NULL;
    goto out_close;

out_free:
    free(g);
    free(aout);
    free(colors);
    g = NULL;
out_close:
    gzclose(f);
out:
    return g;
}
