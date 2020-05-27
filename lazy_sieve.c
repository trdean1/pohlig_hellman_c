#include "ph.h"

//Return the first n primes
//Assumes primes is allocated and is of length n*sizeof(unsigned)
//returns number of primes actually found which should never be different
//unless my padding of prime density estimation is too small...
unsigned lazy_sieve(unsigned n, unsigned *primes) 
{ 
    float n_fl = (float) n;
    unsigned max_sieve = (unsigned) ceil( 1.33 * n_fl * log( n_fl ) );

    char *sieve = malloc( max_sieve * sizeof( char ) ); 
    memset(sieve, 1, max_sieve*sizeof(char)); 
  
    for (unsigned p=2; p*p<=max_sieve; p++) 
    { 
        if (sieve[p] == 1) 
        { 
            for (unsigned i=p*p; i<=max_sieve; i += p) 
                sieve[i] = 0; 
        } 
    } 
  
    unsigned count = 0;
    for (unsigned p=2; p<=max_sieve; p++) { 
        if (sieve[p]) {
            primes[count] = p;
            count++;
        }

        if( count == n ) 
            break;
    }

    free( sieve );
    
    return count;
} 

/****************************************************************************/
  
int test_sieve_basic( int verbose ) 
{ 
    int n = 1000; 
    unsigned *primes = malloc( n * sizeof( unsigned ) );
    unsigned m = lazy_sieve(n, primes); 

    if( verbose ) 
        printf( "Returned %u primes (requested %u) \n", m, n);

    if ( n != m ) {
        return 1;
    }

    if( verbose ) {
    for( unsigned i = 0; i < n; i++ ) {
        printf( "%4u ", primes[i] );
        if ( (i+1) % 10 == 0 ) {
            printf( "\n" );
        }
        }
        printf( "\n" );
    }

    if ( primes[999] == 7919 ) {
        return 0;
    } else {
        return 1;
    }
} 

int test_sieve_return_length( int verbose )
{
    for ( unsigned i = 10; i < 10000; i *= 2 ) {
        unsigned *primes = malloc( i * sizeof( unsigned ) );
        unsigned m = lazy_sieve( i, primes );

        if ( m != i ) {
            if( verbose )
                printf( "\nRequested %u, got %u\n", i, m);
            free( primes );
            return 1;
        }

        free( primes );
    }

    return 0;
}
