#include "ph.h"

int main(int argc, char *argv[]) {
    int verbose = 0;
    if( argc == 2 ) {
        if( !strcmp( argv[1], "-v" ) ) {
            verbose = 1;
        } else {
            fprintf(stderr, "usage: %s [-v]\n", argv[0] );
            exit(1);
        }
    }

    printf("test_sieve_basic...\t");
    if ( !test_sieve_basic( verbose ) ) {
        printf("OK!\n");
    } else {
        printf("Failed\n");
    }

    printf("test_sieve_return_length...\t");
    if ( !test_sieve_return_length( verbose ) ) {
        printf("OK!\n");
    } else {
        printf("Failed\n");
    }

    printf("test_bsgs_basic...\t");
    if ( !test_bsgs_basic( verbose ) ) {
        printf("OK!\n");
    } else {
        printf("Failed\n");
    }

    printf("test_bsgs_small...\t");
    if ( !test_bsgs_small( verbose ) ) {
        printf("OK!\n");
    } else {
        printf("Failed\n");
    }

    printf("test_ph_basic...\t");
    if ( !test_ph_basic( verbose ) ) {
        printf("OK!\n");
    } else {
        printf("Failed\n");
    }

    printf("test_ph_crt...\t");
    if ( !test_ph_crt( verbose ) ) {
        printf("OK!\n");
    } else {
        printf("Failed\n");
    }

}
