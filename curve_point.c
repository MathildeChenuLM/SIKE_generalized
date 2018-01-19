/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
* Contains definition and methods to initialize, clear, 
* print, normalize and format curves and points, with a
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


/* ---------------------------------------------------------
                       CURVE PART
-----------------------------------------------------------*/

void curve_init_none( curve *E ) {
	/* Initializes a curve E with parameters
	(A, C) = ( 0 + i*0, 0 + i*0).
	*/
	fp2_init_none( &(E->A) );
	fp2_init_none( &(E->C) );
}

void curve_init_fp( curve *E, fp2 *A, fp2 *C ) {
	/* Initializes a curve P with parameters (A, C);
	*/
	fp2_init_none( &(E->A) );
	fp2_init_none( &(E->C) );
	fp2_set( &(E->A), A );
	fp2_set( &(E->C), C );
}

void curve_init_mpz( curve *E, mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t p) {
	/* Initializes a curve E with 
	parameters ( a + ib, c + id).
	*/
	fp2 A, C;
	fp2_init( &A, a, b, p );
	fp2_init( &C, c, d, p );

	fp2_init_none( &(E->A) );
	fp2_init_none( &(E->C) );
	fp2_set( &(E->A), &A );
	fp2_set( &(E->C), &C );

	fp2_clear( &A );
	fp2_clear( &C );
}

void curve_init_ui( curve *E, int a, int b, int c, int d, mpz_t p ) {
	/* Initializes a curve E with 
	parameters ( a + ib, c + id), with a,b,c,d of type int.
	*/
	mpz_t mp_a, mp_b, mp_c, mp_d;
	mpz_inits(mp_a, mp_b, mp_c, mp_d, NULL);
	mpz_set_si( mp_a, a );
	mpz_set_si( mp_b, b );
	mpz_set_si( mp_c, c );
	mpz_set_si( mp_d, d );

	curve_init_mpz( E, mp_a, mp_b, mp_c, mp_d, p );

	mpz_clears( mp_a, mp_b, mp_c, mp_d, NULL);
}

void curve_init_normalize( curve *E, fp2 *A, mpz_t p ) {
	/* Initializes a curve E with 
	parameters ( A, 1 ).
	*/
	fp2 one;
	fp2_init_none(&one);
	mpz_set_ui( one.s0, 1);  

	curve_init_fp( E, A, &one );

	fp2_clear(&one);
}

void curve_init_fp2_int( curve *E, fp2 *A, int a, mpz_t p ) {
	/* Initializes a curve E with 
	parameters ( A, a ), with A of type fp2 and a of type int.
	*/
	fp2 tmp;
	fp2_init_none(&tmp);
	mpz_set_ui( tmp.s0, a);  

	curve_init_fp( E, A, &tmp );

	fp2_clear(&tmp);
}

void curve_set( curve *E, curve *F ) {
	/* Sets curve E as curve F.
	*/ 
 	fp2_set( &(E->A), &(F->A) );
 	fp2_set( &(E->C), &(F->C) );
}

void curve_clear( curve *E ) {
	/* Clears fields of curve E.
	*/
	fp2_clear( &(E->A) );
	fp2_clear( &(E->C) );
}

void curve_print( curve *E ) {
	/* Prints E on two lines, 
	respectively E.A and E.C with fp representation.
	*/
	fp2_print( &(E->A) );
	fp2_print( &(E->C) );
	printf("\n");
}

void curve_Aplus_C( curve *F, curve *E, mpz_t p ) {
	/* Sets the curve F such that F = (A+2C, 4C),
	where E = (A, C).
	*/
	fp2 tmp;
	fp2_init_none( &tmp );

	// C = 4C
	fp2_add( &tmp, &(E->C), &(E->C), p );
	fp2_add( &(F->C), &tmp, &tmp, p );

	// Aplus = A + 2*C;
	fp2_add( &(F->A), &tmp, &(E->A), p );

 	fp2_clear(&tmp);
}

void curve_Aplus_Aminus( curve *F, curve *E, mpz_t p ) {
	/* Sets the curve F such that F = (A+2C, A-2C),
	where E = (A, C).
	*/
	fp2 tmp;
	fp2_init_none( &tmp );

	// A = A + 2C
	fp2_add( &tmp, &(E->C), &(E->C), p );
	fp2_add( &(F->A), &(E->A), &tmp, p );

	// C = A - 2C;
	fp2_sub( &(F->C), &(E->A), &tmp, p );

 	fp2_clear(&tmp);
}

void curve_aplus( curve * F, curve *E, mpz_t p ) {
	/* Sets the curve F such that F = (A+2C / 4C, 1),
	where E = (A, C);
	*/
	curve_Aplus_C( F, E, p );
	fp2_multinv( &(F->C), &(F->C), p );
	fp2_mult( &(F->A), &(F->A), &(F->C), p );

	mpz_t un, zero;
	mpz_inits(un, zero, NULL);
	mpz_set_ui( un, 1);
	fp2 one;
	fp2_init( &one, un, zero, p );

	fp2_set( &(F->C), &one );

	fp2_clear( &one );
	mpz_clears(un, zero, NULL);

}

void curve_normalize( curve *F, curve *E, mpz_t p ) {
	/* Sets F as ( A/C, 1 ), where E = ( A, C ).
	Note : C should NEVER be zero.
	*/
	if( fp2_is_zero( &(E->C), p ) ) {
		printf("Don't do curve normalization, C is zero !!\n");
	}
	else {
		mpz_t un, zero;
		mpz_inits(un, zero, NULL);
		mpz_set_ui( un, 1);
		fp2 one;
		fp2_init( &one, un, zero, p );
		curve G;
		curve_init_none(&G);

		fp2_multinv( &(G.A), &(E->C), p );
		fp2_mult( &(G.A), &(E->A), &(G.A), p );
		fp2_set( &(G.C), &one );

		curve_set( F, &G );

		curve_clear(&G);
		fp2_clear( &one );
		mpz_clears(un, zero, NULL);
	}
}

/* ---------------------------------------------------------
                       POINT PART
-----------------------------------------------------------*/

void point_init_none( point *P ) {
	/* Initializes a point P as ( 0 + i*0, 0 + i*0).
	*/
	fp2_init_none( &(P->X) );
	fp2_init_none( &(P->Z) );
}

void point_init_fp( point *P, fp2 *X, fp2 *Z ) {
	/* Initializes a point P as (X, Z).
	*/
	fp2_init_none( &(P->X) );
	fp2_init_none( &(P->Z) );
	fp2_set( &(P->X), X );
	fp2_set( &(P->Z), Z );
}

void point_init_mpz( point *P, mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t p) {
	/* Initializes a point P as ( a + ib, c + id).
	*/
	fp2 X, Z;
	fp2_init( &(X), a, b, p );
	fp2_init( &(Z), c, d, p );

	fp2_init_none( &(P->X) );
	fp2_init_none( &(P->Z) );
	fp2_set( &(P->X), &X );
	fp2_set( &(P->Z), &Z );

	fp2_clear( &(X) );
	fp2_clear( &(Z) );
}

void point_init_ui( point *P, int a, int b, int c, int d, mpz_t p ) {
	/* Initializes a point P with 
	parameters ( a + ib, c + id), with a,b,c,d of type int.
	*/
	mpz_t mp_a, mp_b, mp_c, mp_d;
	mpz_inits(mp_a, mp_b, mp_c, mp_d, NULL);
	mpz_set_ui( mp_a, a );
	mpz_set_ui( mp_b, b );
	mpz_set_ui( mp_c, c );
	mpz_set_ui( mp_d, d );

	point_init_mpz( P, mp_a, mp_b, mp_c, mp_d, p );

	mpz_clears( mp_a, mp_b, mp_c, mp_d, NULL);
}

void point_init_normalize( point *P, fp2 * xP, mpz_t p ) {
	/* Initializes a point P as ( xP, 1 ).
	*/
	fp2 one;
	fp2_init_none(&one);
	mpz_set_ui( one.s0, 1); 

	point_init_fp( P, xP, &one );

	fp2_clear(&one);
}

void point_set( point *R, point *P ) {
	/* Sets point R as point P.
	*/ 
 	fp2_set( &(R->X), &(P->X) );
 	fp2_set( &(R->Z), &(P->Z) );
}

void point_clear( point *P ) {
	/* Clears fields of point P.
	*/
	fp2_clear( &(P->X) );
	fp2_clear( &(P->Z) );
}

void point_print( point *P ) {
	/*Prints P on two lines, 
	respectively P->X and P->y with fp representation.
	*/
	fp2_print( &(P->X) );
	fp2_print( &(P->Z) );
	printf("\n");
}

int points_are_equal( point *P, point *Q ) {
	/* Returns 1 if points are equal, and 0 otherwise.
	
	Note : Might need some normalization before.
	*/
	if( fp2_are_equal(&(P->X), &(Q->X)) == 1 && fp2_are_equal(&(P->Z), &(Q->Z)) == 1 ) {
		return 1;
	}
	else {
		return 0;
	}
}

void point_normalize( point *R, point *P, mpz_t p ) {
	/* Sets R as ( X/Z, 1 ), where P = ( X, Z ), with Z non zero.

	WATCH OUT ! If Z is 0, sets R as (0, 0).
	*/
	if( fp2_is_zero( &(P->Z), p ) ) {
		point Z;
		point_init_none( &Z );
		point_set( R, &Z );
		point_clear(&Z);
	}
	else {
		mpz_t un, zero;
		mpz_inits(un, zero, NULL);
		mpz_set_ui( un, 1);
		fp2 one, tmp;
		fp2_init( &one, un, zero, p );
		fp2_init_none(&tmp); 

		fp2_multinv( &tmp, &(P->Z), p );

		fp2_mult( &(R->X), &(P->X), &tmp, p );
		fp2_set( &(R->Z), &one );

		fp2_clear( &one );
		fp2_clear(&tmp);
		mpz_clears(un, zero, NULL);
	}
}

void point_normalize_X( fp2 *X, point *P, mpz_t p ) {
	/* Sets X as  X/Z, where P = ( X, Z ), with Z non zero.

	WATCH OUT ! If Z is 0, sets R as (0, 0).
	*/
	if( fp2_is_zero( &(P->Z), p ) ) {
		fp2 z;
		fp2_init_none( &z );
		fp2_set( X, &z );
		fp2_clear(&z);
	}
	else {
		mpz_t un, zero;
		mpz_inits(un, zero, NULL);
		mpz_set_ui( un, 1);
		fp2 one;
		fp2_init( &one, un, zero, p );

		fp2_multinv( X, &(P->Z), p );
		fp2_mult( X, &(P->X), X, p );

		fp2_clear( &one );
		mpz_clears(un, zero, NULL);
	}
}

/* 
HEY ! LISTEN : the prime p has to be 3 mod 4, to avoid
having -1 in Fp. Think about it while testing...
Anyway p = 2^e*3^e' - 1
*/

// Main is only for testing !
/*
int main(int argc, char const *argv[])
{	
	mpz_t p, quatre, trois, deux, un;
	mpz_inits(p, quatre, trois, deux, un, NULL);
	mpz_set_ui(p, 7);
	mpz_set_ui(quatre, 4);
	mpz_set_ui(trois, 3);
	mpz_set_ui(deux, 2);
	mpz_set_ui(un, 1);

	fp2 X, Z, A, C;
	fp2_init_none(&X);
	fp2_init(&Z, quatre, quatre, p);
	fp2_init(&A, trois, deux, p);
	fp2_init(&C, deux, trois, p);

	printf("TESTS CURVES :\n");
	curve E, F;
	printf("curve_init_none E = (0: 0) : \n");
	curve_init_none(&E);
	curve_print(&E);

	printf("curve_init_fp E = (3+2i: 2+3i) : \n");
	curve_init_fp(&F, &A, &C);
	curve_print(&F);
	curve_clear(&F);

	printf("curve_init_mpz E = (1+2i: 3+4i) : \n");
	curve_init_mpz(&F, un, deux, trois, quatre, p);
	curve_print(&F);
	curve_clear(&F);

	printf("curve_init ui E = (1+2i: 3+4i) : \n");
	curve_init_ui(&F, 1, 2, 3, 4, p);
	curve_print(&F);
	curve_clear(&F);

	printf("curve_init_normalize E = (3+2i: 1) : \n");
	curve_init_normalize(&F, &A, p);
	curve_print(&F);
	curve_clear(&F);

	printf("curve_init_fp2_int E (3+2i: 4) : \n");
	curve_init_fp2_int(&F, &A, 4, p);
	curve_print(&F);

	printf("curve_set F = E = (3+2i: 4) : \n");
	curve_set(&E, &F);
	curve_print(&E);

	printf("curve_AplusC E = (4+2i : 2) : \n");
	curve_Aplus_C(&F, &E, p);
	curve_print(&F);

	printf("curve_Aplus_Aminus E = (4+2i : 2+2i) : \n");
	curve_Aplus_Aminus(&F, &E, p);
	curve_print(&F);

	printf("curve_aplus E = (2+i : 1) : \n");
	curve_aplus(&F, &E, p);
	curve_print(&F);

	printf("curve_normalize E = (6+4i : 1)\n");
	curve_normalize(&F, &E, p);
	curve_print(&F);

	curve_clear(&E);
	curve_clear(&F);


	printf("\nTESTS POINTS : \n");
	point P, Q;
	printf("point_init_none P = (0: 0) : \n");
	point_init_none(&P);
	point_print(&P);

	printf("point_init_fp P = (0: 4+4i) : \n");
	point_init_fp(&Q, &X, &Z);
	point_print(&Q);
	point_clear(&Q);

	printf("point_init_normalize P = (4+4i: 1) : \n");
	point_init_normalize(&Q, &Z, p);
	point_print(&Q);
	point_clear(&Q);

	printf("point_init_mpz P = (4+2i: 2+4i) : \n");
	point_init_mpz(&Q, quatre, deux, deux, quatre, p);
	point_print(&Q);
	point_clear(&Q);

	printf("point_init_ui P = (4+2i: 2+4i) : \n");
	point_init_ui(&Q, 4, 2, 2, 4, p);
	point_print(&Q);
	
	int boolean;
	printf("P and Q are equal ? Should return 0.\n");
	boolean = points_are_equal( &P, &Q );
	printf("points are equal : %d\n\n", boolean);
	
	printf("Set Q = P. P and Q are equal ? Should return 1.\n");
	point_set(&P, &Q);
	boolean = points_are_equal( &P, &Q );
	printf("points are equal : %d\n\n", boolean);

	printf("point normalize P = (5+5i : 1) : \n");
	point_normalize( &P, &Q, p );
	point_print(&P);

	printf("point normalize_X P = 5+5i  : \n");
	point_normalize_X( &X, &Q, p );
	fp2_print(&X);

	point R;
	point_init_ui( &R, 1, 2, 0, 0, p );
	printf("\nR : \n");
	point_print(&R);

	printf("point normalize_X R = 0  : \n");
	point_normalize_X( &X, &R, p );
	fp2_print(&X);

	printf("point normalize R = (0 : 0) : \n");
	point_normalize( &R, &R, p );
	point_print(&R);

	point_clear(&P);
	point_clear(&Q);
	point_clear(&R);
	fp2_clear(&X);
	fp2_clear(&Z);
	fp2_clear(&A);
	fp2_clear(&C);
	mpz_clears(p, quatre, trois, deux, un, NULL);
}
*/
