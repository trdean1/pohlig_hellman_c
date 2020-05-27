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
