#include "ph.h"

void init_bsgs( bsgs_t *b, const mpz_t base, const mpz_t modulus, const unsigned long order )
{
    float s = ceil( sqrt( (float) order ) );
    unsigned long sqrt_order = (unsigned long) s;

    init_vec( &(b->bs_vals), sqrt_order );

    for ( unsigned long i = 0; i < sqrt_order; i++ ) {
        mpz_powm_ui( b->bs_vals.el[i], base, i, modulus );
    }

    mpz_t base_inv;
    mpz_init( base_inv );
    mpz_invert( base_inv, base, modulus );

    mpz_t base_minus_m;
    mpz_init( base_minus_m );
    mpz_powm_ui( base_minus_m, base_inv, sqrt_order, modulus );

    mpz_init_set( b->modulus, modulus );
    mpz_init_set( b->base, base );
    b->order = order;
    b->sqrt_order = sqrt_order;
    mpz_init_set( b->a_minus_m, base_minus_m );
    
    mpz_clear( base_inv );
    mpz_clear( base_minus_m );
}

void free_bsgs( bsgs_t *b )
{
    mpz_clear( b->modulus );
    mpz_clear( b->base );
    mpz_clear( b->a_minus_m );
    free_vec( &(b->bs_vals) );
}

int compute_dlog( bsgs_t *b, unsigned long *rop, const mpz_t value )
{
    mpz_t gamma;
    mpz_init_set( gamma, value );

    unsigned long i = 0;
    unsigned long j = 0;
    while( mpz_cmp( b->bs_vals.el[i], gamma ) ) {
        i++;

        if( i == b->sqrt_order ) {
            i = 0;
            j++;
            mpz_mul( gamma, gamma, b->a_minus_m );
            mpz_mod( gamma, gamma, b->modulus );
        }

        if( j*(b->sqrt_order)+i > b->order ) {
            mpz_clear( gamma );
            return 1;
        }
    }

    *rop = j*(b->sqrt_order)+i;
    return 0;
}

/*******************************************************************/

void print_bsgs_t( bsgs_t *b ) {
    fprint_bsgs_t( stdout, b );
}

void fprint_bsgs_t( FILE *f, bsgs_t *b )
{
    fprintf(f, "bsgs_t\n\n");

    gmp_fprintf( f, " modulus: %Zd\n base: %Zd\n order: %lu\n sqrt_order:%lu\n a_minus_m: %Zd\n",
                 b->modulus, b->base, b->order, b->sqrt_order, b->a_minus_m );

    fprintf( f, "\n base^i\n" );
    for ( unsigned i = 0; i < b->sqrt_order; i++ ) {
        gmp_fprintf( f, "  %u: %Zd\n", i, b->bs_vals.el[i] );
    }
    fprintf( f, "\n\n" );
}  

int test_bsgs_basic( verbose ) 
{
    int ret_val = 0;

    bsgs_t b;

    mpz_t base;
    mpz_init_set_ui( base, 2 );

    mpz_t modulus;
    mpz_init_set_ui( modulus, 11 );

    init_bsgs( &b, base, modulus, 11 );

    if( b.sqrt_order != 4 || 
        mpz_cmp_ui( b.a_minus_m, 9 ) || 
        mpz_cmp_ui( b.bs_vals.el[3], 8 ) ) 
    {
        if( verbose )
            print_bsgs_t( &b );

        ret_val = 1;
    }     

    unsigned long res;
    mpz_t value;
    mpz_init_set_ui( value, 6 );
    if( compute_dlog( &b, &res, value ) ) {
        if( verbose )
            printf( "Did not find log\n" );
        ret_val = 1;
    }

    if( res != 9 ) {
        if( verbose ) 
            printf( "Returned %lu, expected 9\n", res );
        ret_val = 1;
    }

    free_bsgs( &b );

    return ret_val;
}

int test_bsgs_small( verbose )
{
    int ret_val = 0;

    bsgs_t b;

    mpz_t base;
    mpz_init_set_ui( base, 3 );

    mpz_t modulus;
    mpz_init_set_ui( modulus, 77 );

    init_bsgs( &b, base, modulus, 30 );

    unsigned long res;
    mpz_t value;
    mpz_init_set_ui( value, 69 );   
    compute_dlog( &b, &res, value );

    if( res != 21 ) {
        if( verbose ) { 
            printf( "Returned %lu, expected 21\n", res );
            print_bsgs_t( &b );
        }
        ret_val = 1;
    }

    mpz_init_set_ui( value, 58 );   
    compute_dlog( &b, &res, value );

    if( res != 26 ) {
        if( verbose ) 
            printf( "Returned %lu, expected 26\n", res );
            print_bsgs_t( &b );
        ret_val = 1;
    }

    free_bsgs( &b );

    mpz_init_set_ui( modulus, 82 );
    mpz_init_set_ui( value, 55 );
    init_bsgs( &b, base, modulus, 8 );
    compute_dlog( &b, &res, value );

    if( res != 7 ) {
        if( verbose ) {
            printf( "Returned %lu, expected 7\n", res );
            print_bsgs_t( &b );
        }
        ret_val = 1;
    }

    free_bsgs( &b );

    return ret_val;
}
