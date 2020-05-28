#include "ph.h"

void init_ph( ph_t *rop,
              const mpz_t generator,
              const unsigned primorial_order )
{
    rop->primes = malloc( primorial_order * sizeof( unsigned long ) );
    assert( lazy_sieve( primorial_order, rop->primes ) == primorial_order );
    rop->n_primes = primorial_order;

    //Modulus = \prod p_i
    mpz_init_set_ui( rop->modulus, rop->primes[0] );
    for( unsigned i = 1; i < primorial_order; i++ ) {
        mpz_mul_ui( rop->modulus, rop->modulus, rop->primes[i] );
    }

    mpz_init_set( rop->generator, generator );

    init_vec( &( rop->subgroup_generators ), rop->n_primes );
    init_vec( &( rop->mpi ), rop->n_primes );

    //mpi = M / p_i
    for( unsigned i = 0; i < rop->n_primes; i++ ) {
        mpz_tdiv_q_ui( rop->mpi.el[i], rop->modulus, rop->primes[i] );
    }

    //precompute mpi^{-1} mod p_i for CRT
    init_vec( &( rop->mpi_inv ), rop->n_primes );
    mpz_t p; 
    mpz_init( p ); 
    for( unsigned i = 0; i < rop->n_primes; i++ ) {
        mpz_set_ui( p, rop->primes[i] );
        mpz_invert( rop->mpi_inv.el[i], rop->mpi.el[i], p );
    }

    //compute g_i = g^(M/p_i)
    for ( unsigned i = 0; i < rop->n_primes; i++ ) {
        mpz_powm( rop->subgroup_generators.el[i], generator, rop->mpi.el[i], rop->modulus );
    }

    //setup bsgs algorithm for each subgroup
    rop->subgroup_bsgs = malloc(  rop->n_primes * sizeof( bsgs_t ) );
    for( unsigned i = 0; i < rop->n_primes; i++ ) {
        init_bsgs( &( rop->subgroup_bsgs[i] ), 
                   rop->subgroup_generators.el[i],
                   rop->modulus,
                   rop->primes[i] );

    }
    
}

void free_ph( ph_t *op ) 
{
    mpz_clear( op->modulus );
    mpz_clear( op->generator );
    free_vec( &( op->subgroup_generators ) );
    free_vec( &( op->mpi ) );
    free_vec( &( op->mpi_inv ) );
    
    for( unsigned i = 0; i < op->n_primes; i++ ) 
        free_bsgs( &( op->subgroup_bsgs[i] ) );

    free( op->primes );
    memset( op, 0, sizeof( ph_t ) );
}

//Assumes rop and op are properly initialized
void compute_subgroup_dlogs( unsigned long *rop, 
                             ph_t *op,
                             mpz_t value )
{
    vec_t vals;
    init_vec( &vals, op->n_primes );

    //h_i = h^{M/p_i} mod M
    for ( unsigned i = 0; i < op->n_primes; i++ ) {
        mpz_powm( vals.el[i], value, op->mpi.el[i], op->modulus );

    }

    //solve g_i^{x_i} = h_i mod M
    for( unsigned i = 0; i < op->n_primes; i++ ) {
        compute_dlog( &( op->subgroup_bsgs[i] ),
                      &( rop[i] ),
                      vals.el[i] );


    }

    free_vec( &vals );
}

//Assumes rop is initialized
void solve_congruences( mpz_t rop,
                        ph_t *op,
                        unsigned long *x_i ) 
{
     mpz_t tmp;
     mpz_init( tmp );
     mpz_set_ui( rop, 0 );

     //sum x_i * mpi * mpi_inv 
     for( unsigned i = 0; i < op->n_primes; i++ ) {
        mpz_mul( tmp, op->mpi.el[i], op->mpi_inv.el[i] );
        mpz_mul_ui( tmp, tmp, x_i[i] );
        mpz_mod( tmp, tmp, op->modulus );
        mpz_add( rop, rop, tmp ); 
     }

    mpz_mod( rop, rop, op->modulus );
}

/****************************************************************************/

void print_ph( ph_t *op )
{
    printf("-----------ph_t----------\n");
    gmp_printf("Modulus:\t\t%Zd\n", op->modulus);
    gmp_printf("Generator:\t\t%Zd\n", op->generator);
    printf("Subgroup generators:\n");
    pp_vec( &(op->subgroup_generators) );
    printf("Mpi:\n");
    pp_vec( &(op->mpi) );
    printf("mpi_inv:\n");
    pp_vec( &(op->mpi_inv) );

    printf("\nbsgs_t's:\n");
    for( unsigned i = 0; i < op->n_primes; i++ ) {
        printf("(%u, %lu):\n", i, op->primes[i] );
        print_bsgs_t( &( op->subgroup_bsgs[i] ) );
        printf("\n");
    }
}

/****************************************************************************/

int test_ph_basic( int verbose ) 
{
    int ret_val = 0;

    unsigned long order = 4;

    unsigned long exp_sg[] = {153, 39, 99, 99};
    unsigned long exp_mpi[] = {105, 70, 42, 30};
    unsigned long exp_mpi_inv[] = {1, 1, 3, 4};

    mpz_t generator;
    mpz_init_set_ui( generator, 3 );

    ph_t ph;

    init_ph( &ph, generator, order );

    if( verbose )
        print_ph( &ph );

    if( mpz_cmp_ui( ph.modulus, 210 ) )
        ret_val = 1;

    if( cmp_vec_ui( &ph.subgroup_generators, exp_sg, 4 ) )
        ret_val = 2;

    if( cmp_vec_ui( &ph.mpi, exp_mpi, 4 ) )
        ret_val = 3;

    if( cmp_vec_ui( &ph.mpi_inv, exp_mpi_inv, 4 ) )
        ret_val = 4;

    free_ph( &ph );

    return ret_val;
}

int test_ph_crt( int verbose ) 
{
    int ret_val = 0;

    unsigned long order = 4;

    unsigned long xi[] = {1, 0, 2, 1};

    mpz_t generator;
    mpz_init_set_ui( generator, 3 );

    ph_t ph;

    init_ph( &ph, generator, order );

    mpz_t res;
    mpz_init( res );
    
    solve_congruences( res, &ph, xi );

    if( mpz_cmp_ui( res, 57 ) ) {
        if( verbose )
            gmp_printf("Expected 57, got %Zd\n", res);
        ret_val = 1;
    }

    free_ph( &ph );

    return ret_val;
}
