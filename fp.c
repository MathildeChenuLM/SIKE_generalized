/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
* Contains definition and methods to initialize, clear, 
* print, normalize and format elements of Fp², with a
* set of tests. 
*
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include "Header.h"


void fp_inv( mpz_t *res, mpz_t a, mpz_t p ) {
	/* Sets res with the invrese of a mod p.
	*/
	mpz_t g, t;
	mpz_inits( g, t, NULL );
	mpz_gcdext ( g, *res, t, a, p );
	mpz_clears( g, t, NULL );
}

void fp2_init_none( fp2 *res ) {
	/* Initializes res with value 0 for every fields.
	*/
	mpz_inits( res->s0, res->s1, NULL );
}

void fp2_init( fp2 *res, mpz_t x, mpz_t y, mpz_t p ) {
	/* Initializes res with fields s0 = x and s1 = y.
	*/
	mpz_t tmp; 
	mpz_init(tmp);
	mpz_inits( res->s0, res->s1, NULL );

	// The tmp variable is necessary to avoid 
	// changing x, y, or (worse !) p value.
	mpz_mod( tmp, x, p );
	mpz_set( res->s0, tmp );
	mpz_mod( tmp, y, p );
	mpz_set( res->s1, tmp );

	mpz_clear(tmp);
}

void fp2_clear( fp2 *res ) {
	/* Clears the fields of res.
	*/
	mpz_clears( res->s0, res->s1, NULL );
}

void fp2_print( fp2 *my_fp2 ) {
	/* Prints my_fp2 fields. 
	*/
	gmp_printf("%Zd + %Zd i \n", my_fp2->s0, my_fp2->s1);
}

void fp2_add( fp2 *res, fp2 *a, fp2 *b, mpz_t p) {
	/* Adds a and b and writes the result in res.
	*/
	/*
	res->s0 = (a.s0 + b.s0) % p;
	res->s1 = (a.s1 + b.s1) % p;
	*/
	mpz_t tmp1, tmp2;
	mpz_inits( tmp1, tmp2, NULL );
	mpz_add( tmp1, a->s0, b->s0 );
	mpz_mod( tmp1, tmp1, p);
	mpz_add( tmp2, a->s1, b->s1 );
	mpz_mod( tmp2, tmp2, p);

	mpz_set( res->s0, tmp1 );
	mpz_set( res->s1, tmp2);

	mpz_clears(tmp1, tmp2, NULL);
}

void fp2_sub( fp2 *res, fp2 *a, fp2 *b, mpz_t p) {
	/* Substractes a and b and writes the result in res.
	*/
	/*
	res->s0 = (a.s0 - b.s0) % p;
	res->s1 = (a.s1 - b.s1) % p;
	*/
	mpz_t tmp1, tmp2, minus_b;
	mpz_inits( tmp1, tmp2, minus_b, NULL );

	mpz_neg( minus_b, b->s0 );
	mpz_add( tmp1, a->s0, minus_b );
	mpz_mod( tmp1, tmp1, p);
	mpz_neg( minus_b, b->s1 );
	mpz_add( tmp2, a->s1, minus_b );
	mpz_mod( tmp2, tmp2, p);

	mpz_set( res->s0, tmp1 );
	mpz_set( res->s1, tmp2);

	mpz_clears(tmp1, tmp2, minus_b, NULL);
}

void fp2_mult( fp2 *res, fp2 *a, fp2 *b, mpz_t p ) {
	/* Multiplies a and b and writes the result in res.
	*/
	/*
	res->s0 = ( a.s0 * b.s0 - a.s1 * b.s1 ) % p;
	res->s1 = ( a.s0 * b.s1 + a.s1 * b.s0 ) % p;
	*/
	mpz_t tmp1, tmp2, tmpL, tmpR;
	mpz_inits( tmp1, tmp2, tmpL, tmpR, NULL );

	mpz_mul( tmpL, a->s0, b->s0 );
	mpz_mul( tmpR, a->s1, b->s1 );
	mpz_sub( tmp1, tmpL, tmpR );
	mpz_mod( tmp1, tmp1, p );

	mpz_mul( tmpL, a->s0, b->s1 );
	mpz_mul( tmpR, a->s1, b->s0 );
	mpz_add( tmp2, tmpL, tmpR );
	mpz_mod( tmp2, tmp2, p );

	mpz_set( res->s0, tmp1 );
	mpz_set( res->s1, tmp2 );

	mpz_clears(tmp1, tmp2, tmpL, tmpR, NULL);

}

void fp2_addinv( fp2 *res, fp2 *a, mpz_t p ) {
	/* Fills res with the additive inverse of a. 
	Note : The gmp_mod function always returns a positive value.
	*/
	/*
	res->s0 = (- a.s0) %p;
	res->s1 = (- a.s1) %p;
	*/
	mpz_t tmp1, tmp2;
	mpz_inits( tmp1, tmp2, NULL );

	mpz_neg( tmp1, a->s0 );
	mpz_mod( tmp1, tmp1, p );
	mpz_neg( tmp2, a->s1 );
	mpz_mod( tmp2, tmp2, p );

	mpz_set( res->s0, tmp1 );
	mpz_set( res->s1, tmp2 );

	mpz_clears( tmp1, tmp2, NULL );

}

void fp2_multinv( fp2 *res, fp2 *a, mpz_t p ) {
	/* Fills res with the multiplicative inverse of a, ie res*a = 1.
	*/
	/*
	res->s0 = s0 * ( s0² + s1² )^(-1)
	res->s1 = -s1 * ( s0² + s1² )^(-1)
	*/
	mpz_t tmp1, tmp2, tmpL, tmpR, g, t;
	mpz_inits( tmp1, tmp2, tmpL, tmpR, g, t, NULL );

	mpz_mul( tmpL, a->s0, a->s0 );
	mpz_mul( tmpR, a->s1, a->s1 );
	mpz_add( tmpR, tmpL, tmpR );
	mpz_mod( tmpR, tmpR, p );
	mpz_gcdext( g, tmpR, t, tmpR, p ); // Multiplicative inverse of tmpR mod p.

	mpz_mul( tmp1, a->s0, tmpR );
	mpz_mod( tmp1, tmp1, p );
	mpz_neg( tmp2, a->s1 );
	mpz_mul( tmp2, tmp2, tmpR );
	mpz_mod( tmp2, tmp2, p );

	mpz_set( res->s0, tmp1 );
	mpz_set( res->s1, tmp2 );

	mpz_clears(tmp1, tmp2, tmpL, tmpR, g, t, NULL);
}

int fp2_are_equal( fp2 *a, fp2 *b) {
	/* Returns 1 if a and b are equal, and 0 otherwise.
	*/
	if( mpz_cmp(a->s0, b->s0) == 0 && mpz_cmp(a->s1, b->s1) == 0 ) {
		return 1;
	}
	else {
		return 0;
	}
}

int fp2_is_zero( fp2 *a, mpz_t p ) {
	/* Returns 1 if a is (0 : 0), and 0 otherwise.
	*/
	mpz_mod( a->s0, a->s0, p );
	mpz_mod( a->s1, a->s1, p);
	if( mpz_cmp_ui(a->s0, 0) == 0 && mpz_cmp_ui(a->s1, 0) == 0 ) {
		return 1;
	}
	else {
		return 0;
	}
}

void fp2_set( fp2 *res, fp2 *a ) {
	/* Sets fields from res with fields from a.
	*/
	mpz_set( res->s0, a->s0 );
	mpz_set( res->s1, a->s1 );
}

/* 
HEY ! LISTEN : the prime p has to be 3 mod 4, to avoid
having -1 in Fp. Think about it while testing...
Anyway p = 2^e*3^e' - 1.
*/

// The main is only for testing
/*
int main(int argc, char const *argv[])
{	
	printf("Start ! \n");
	fp2 a, b, c;
	mpz_t p, un, deux;

	// Test fp2_init_none, fp2_clear and fp2_print.
	fp2_init_none(&a);
	printf("Null a : \n");
	fp2_print(&a);
	fp2_clear(&a);
	printf("\n");

	// Test fp2_init.
	mpz_inits(p, un, deux, NULL);
	mpz_set_ui(p, 7);
	mpz_set_ui(un, 1);
	mpz_set_ui(deux, 2);
	fp2_init(&a, un, deux, p);
	printf("a = 1 + 2i \n");
	fp2_print(&a);
	printf("\n");

	// Test fp2_add.
	fp2_init_none(&b);
	fp2_add(&b, &a, &a, p);
	printf("b = a+a = 2 + 4i \n");
	fp2_print(&b);
	printf("\n");

	// Test fp2_sub.
	fp2_sub(&a, &b, &a, p);
	printf("b-a = 2a-a = 1 + 2i\n");
	fp2_print(&a);
	printf("\n");

	// Test fp2_mult.
	fp2_mult(&b, &a, &b, p);
	printf("a*b = 1 + i\n");
	fp2_print(&b);
	printf("\n");

	// Test fp2_addinv.
	fp2_addinv(&b, &a, p);
	fp2_add(&b, &a, &b, p);
	printf("(-a) + a = 0\n");
	fp2_print(&b);
	printf("\n");

	// Test fp2_multinv.
	fp2_multinv(&b, &a, p);
	fp2_mult(&b, &a, &b, p);
	printf("(1/a) * a = 1\n");
	fp2_print(&b);
	printf("\n");

	// Test fp2_set et fp2_are_equal.
	int boolean;
	printf("Are a and 1 equal ? Should return 0.\n");
	boolean = fp2_are_equal( &a, &b );
	printf("%d\n", boolean);
	printf("Set b = a. Are a and b equal ? Should return 1.\n");
	fp2_set(&b, &a);
	boolean = fp2_are_equal( &a, &b );
	printf("%d\n", boolean);
	printf("\n");

	fp2_init_none( &c );
	printf("Is_zero a ? Should return 0.\n");
	boolean = fp2_is_zero( &a, p );
	printf("%d\n", boolean);
	printf("Is_zero c ? Should return 1.\n");
	boolean = fp2_is_zero( &c, p );
	printf("%d\n", boolean);

	mpz_clears(p, un, deux, NULL);
	fp2_clear(&a);
	fp2_clear(&b);
	fp2_clear(&c);

	printf("Finish !\n");

}
*/

