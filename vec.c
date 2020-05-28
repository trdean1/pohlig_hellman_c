#include "ph.h"

void init_vec(vec_t *v, int count)
{
	assert(v);
	v->count = count;
	v->el = malloc(count * sizeof(mpz_t));
	assert(v->el);
	for (int i=0; i < v->count; i++)
		mpz_init(v->el[i]);
}

void free_vec(vec_t *v)
{
	assert(v);
	for (int i=0; i < v->count; i++)
		mpz_clear(v->el[i]);
	free(v->el);
}

void pp_vec( vec_t *v ) {
    printf("vec_t: %d elements", v->count );
    int i;
    for( i = 0; i < v->count; i++ ) {
        if( i % 5 == 0 ) 
            printf("\n\t");
        gmp_printf("%Zd  ", v->el[i] );
    }
    if( i % 5 != 0 )
        printf("\n");
}

int cmp_vec( vec_t *a, vec_t *b )
{
    if( a->count != b->count ) 
        return -1;

    for( unsigned i = 0; i < a->count; i++ ) {
        if( mpz_cmp( a->el[i], b->el[i] ) )
            return 1;
    }

    return 0;
}

int cmp_vec_ui( vec_t *a, unsigned long *b, int count )
{
    if( a->count != count ) 
        return -1;

    for( unsigned i = 0; i < a->count; i++ ) {
        if( mpz_cmp_ui( a->el[i], b[i] ) )
            return 1;
    }

    return 0;
}
