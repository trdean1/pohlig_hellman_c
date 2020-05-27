#pragma once

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <assert.h>
#include <gmp.h>

typedef struct vec_ {
    mpz_t *el;
    int count;
} vec_t;

typedef struct bsgs_ {
    mpz_t modulus;
    mpz_t base;
    unsigned long order;
    unsigned long sqrt_order;
    mpz_t a_minus_m;
    vec_t bs_vals;
} bsgs_t;

/***********************************************/

void init_vec( vec_t *v, 
               int count );

void free_vec( vec_t *v );

/***********************************************/

void init_bsgs( bsgs_t *b, 
                const mpz_t base, 
                const mpz_t modulus, 
                const unsigned long order );

void free_bsgs( bsgs_t *b );

int compute_dlog( bsgs_t *b,
                  unsigned long *rop,
                  const mpz_t value );

void print_bsgs_t( bsgs_t *b );

void fprint_bsgs_t( FILE *f,
                    bsgs_t *b );

/***********************************************/

unsigned lazy_sieve( unsigned n, 
                     unsigned *primes ); 

/***********************************************/
////Test functions
int test_sieve_basic( int verbose );
int test_sieve_return_length( int verbose );
int test_bsgs_basic( int verbose );
