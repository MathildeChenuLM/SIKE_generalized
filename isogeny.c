/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
* Contains methods to compute the image from curves
* and points by an l-isogeny, in regular (2,3) case
* and generalized (l_a, l_b) case, with a set of 
* tests.
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
                     REGULAR PART
-----------------------------------------------------------*/

void four_iso_curve( curve *F, fp2 *K1, fp2 *K2, fp2 *K3, 
	point *P4, mpz_t p ) {
	/* Sets F as the four-isogenous curve, ie F = E/<P4>,
	where P4 has order 4.

	WATCH OUT : F will be in the form AplusC = (A+2C : 4C)
	*/
	fp2_sub( K2, &(P4->X), &(P4->Z), p );
	fp2_add( K3, &(P4->X), &(P4->Z), p );
	fp2_mult( K1, &(P4->Z), &(P4->Z), p );

	fp2_add( K1, K1, K1, p );
	fp2_mult( &(F->C), K1, K1, p );
	fp2_add( K1, K1, K1, p );

	fp2_mult( &(F->A), &(P4->X), &(P4->X), p );
	fp2_add( &(F->A), &(F->A), &(F->A), p );
	fp2_mult( &(F->A), &(F->A), &(F->A), p );

}

void four_iso_eval( point* Q, 
	fp2 *K1, fp2 *K2, fp2 *K3, point *P, mpz_t p ) {
	/* Sets Q as phi(P), where phi is a four-isogeny,
	described by K1, K2 and K3.

	Note : the spec gives Q', but after checking with source code, it's just Q.
	*/

	fp2 t0, t1;
	fp2_init_none(&t0);
	fp2_init_none(&t1);

	point T;
	point_init_none(&T);
	point_set(&T, P);

	fp2_add( &t0, &(T.X), &(T.Z), p );
	fp2_sub( &t1, &(T.X), &(T.Z), p );
	fp2_mult( &(T.X), &t0, K2, p );
	fp2_mult( &(T.Z), &t1, K3, p );

	fp2_mult( &t0, &t0, &t1, p );
	fp2_mult( &t0, &t0, K1, p );
	fp2_add( &t1, &(T.X), &(T.Z), p );
	fp2_sub( &(T.Z), &(T.X), &(T.Z), p );

	fp2_mult( &t1, &t1, &t1, p );
	fp2_mult( &(T.Z), &(T.Z), &(T.Z), p );
	fp2_add( &(T.X), &t0, &t1, p );
	fp2_sub( &t0, &(T.Z), &t0, p );

	fp2_mult( &(T.X), &(T.X), &t1, p );
	fp2_mult( &(T.Z), &(T.Z), &t0, p );

	point_set(Q, &T);

	point_clear(&T);
	fp2_clear(&t0);
	fp2_clear(&t1);
}

void three_iso_curve( curve *F, fp2 *K1, fp2 *K2, 
	point* P, mpz_t p ) {
	/* Sets F as the three-isogenous curve, ie F = E/<P>,
	where P has order 3.
	WATCH OUT : curve will be in the form (A+ : A-) = (A+2C : A-2C)
	*/
	fp2 t0, t1, t2, t3, t4;
	fp2_init_none(&t0);
	fp2_init_none(&t1);
	fp2_init_none(&t2);
	fp2_init_none(&t3);
	fp2_init_none(&t4);

	fp2_sub( K1, &(P->X), &(P->Z), p );
	fp2_mult( &t0, K1, K1, p );
	fp2_add( K2, &(P->X), &(P->Z), p );
	fp2_mult( &t1, K2, K2, p );
	fp2_add( &t2, &t0, &t1, p );
	fp2_add( &t3, K1, K2, p );

	fp2_mult( &t3, &t3, &t3, p );
	fp2_sub( &t3, &t3, &t2, p );
	fp2_add( &t2, &t1, &t3, p );
	fp2_add( &t3, &t3, &t0, p );
	fp2_add( &t4, &t3, &t0, p );
	fp2_add( &t4, &t4, &t4, p );

	fp2_add( &t4, &t1, &t4, p );
	fp2_mult( &(F->C), &t2, &t4, p ); // pass auf, C is Aminus !
	fp2_add(&t4, &t1, &t2, p );
	fp2_add( &t4, &t4, &t4, p );
	fp2_add( &t4, &t0, &t4, p );
	fp2_mult( &t4, &t3, &t4, p );

	fp2_sub( &t0, &t4, &(F->C), p );
	fp2_add( &(F->A), &(F->C), &t0, p);

	fp2_clear(&t0);
	fp2_clear(&t1);
	fp2_clear(&t2);
	fp2_clear(&t3);
	fp2_clear(&t4);
}

void three_iso_eval( point *phiP, 
	fp2 *K1, fp2* K2, point *P, mpz_t p ) {
	/* Sets phiP as phi(P), where phi is a three-isogeny,
	described by K1 and K2.
	Note : Does not follows exactly the specification, error on line 3.
	*/
	fp2 t0, t1, t2;
	fp2_init_none(&t0);
	fp2_init_none(&t1);
	fp2_init_none(&t2);

	point T;
	point_init_none(&T);
	point_set(&T, P);

	fp2_add(&t0, &(T.X), &(T.Z), p );
	fp2_sub( &t1, &(T.X), &(T.Z), p );
	fp2_mult( &t0, K1, &t0, p ); //typo dans la spec !

	fp2_mult(&t1, K2, &t1, p );
	fp2_add( &t2, &t0, &t1, p );
	fp2_sub( &t0, &t1, &t0, p );

	fp2_mult( &t2, &t2, &t2, p );
	fp2_mult( &t0, &t0, &t0, p );
	fp2_mult( &(T.X), &(T.X), &t2, p );

	fp2_mult( &(T.Z), &(T.Z), &t0, p );

	point_set( phiP, &T );

	point_clear(&T);
	fp2_clear(&t0);
	fp2_clear(&t1);
	fp2_clear(&t2);
}

void two_e_iso( curve *F, point *phiP1, point *phiP2, point *phiP3,
	curve *E, point *S, point *P1, point *P2, point *P3, 
	int e2, mpz_t p ) {
	/* Sets F as the (2^e2)-isogenous curve, ie F = E/<S>,
	where S has order 2^e2 in E.

	WATCH OUT ! E and F are in AplusC form.
	*/
	point T;
	point_init_none(&T);

	point_set(phiP1, P1);
	point_set(phiP2, P2);
	point_set(phiP3, P3);

	curve G;
	curve_init_none(&G);
	curve_set( &G, E );

	fp2 K1, K2, K3;
	fp2_init_none(&K1);
	fp2_init_none(&K2);
	fp2_init_none(&K3);

	for( int e = e2-2; e >=0; e = e - 2 ) {
		xDBLe( &T, S, &G, p, e );
		four_iso_curve( &G, &K1, &K2, &K3, &T, p );
		four_iso_eval( S, &K1, &K2, &K3, S, p );
		// optional :
		four_iso_eval( phiP1, &K1, &K2, &K3, phiP1, p );
		four_iso_eval( phiP2, &K1, &K2, &K3, phiP2, p );
		four_iso_eval( phiP3, &K1, &K2, &K3, phiP3, p );
	}

	curve_set( F, &G );

	fp2_clear(&K1);
	fp2_clear(&K2);
	fp2_clear(&K3);
	curve_clear(&G);
	point_clear(&T);
}

void three_e_iso( curve *F, point *phiP1, point *phiP2, point *phiP3,
	curve *E, point *S, point *P1, point *P2, point *P3, 
	int e3, mpz_t p ) {
	/* Sets F as the (3^e3)-isogenous curve, ie F = E/<S>,
	where S has order 3^e3 in E.
	WATCH OUT E and F are in the form AplusAmoins !
	*/
	point T;
	point_init_none(&T);

	point_set(phiP1, P1);
	point_set(phiP2, P2);
	point_set(phiP3, P3);

	curve G;
	curve_init_none(&G);
	curve_set( &G, E );

	fp2 K1, K2;
	fp2_init_none(&K1);
	fp2_init_none(&K2);

	for( int e = e3-1; e >=0; e-- ) {
		xTPLe( &T, S, &G, p, e );
		three_iso_curve( &G, &K1, &K2, &T, p );
		three_iso_eval( S, &K1, &K2, S, p );
		// optional :
		three_iso_eval( phiP1, &K1, &K2, phiP1, p );
		three_iso_eval( phiP2, &K1, &K2, phiP2, p );
		three_iso_eval( phiP3, &K1, &K2, phiP3, p );
	}

	curve_set( F, &G );
	
	fp2_clear(&K1);
	fp2_clear(&K2);
	curve_clear(&G);
	point_clear(&T);
}

/* ---------------------------------------------------------
                     GENRALIZED PART
-----------------------------------------------------------*/

void kernel_point( int d, point ker[d], point *G, curve *F, mpz_t p ) {
	/*
	Computes the d points from the subgroup of the curve F 
	generated by G, with G a [2d+1] torsion point. Sets ker as <G>.

	Asserts d > 0 (won't work for two torsion points, who are alone anyway).
	The kernel is of size [2d+1], but stored in an array of size d since
	zero is not memorized, and P = -P on the Kummer line.

	WATCH OUT ! F has to be in the form (A+2C : 4C) !
	*/
	point_set( &(ker[0]), G );

	if( d >= 2 ){
		xDBL( &(ker[1]), G, F, p );
	}

	for( int i=3; i<=d; i++ ) {
		xADD( &(ker[i-1]), &(ker[i-2]), G, &(ker[i-3]), p );
	}
}

void kernel_reshape( int d, point ker[d], mpz_t p ) {
	/*
	Reshapes the points in the kernel as (X+Z : X-Z).
	The kernel is then ready for odd_isogeny.
	*/
	fp2 t1, t2;
	fp2_init_none( &t1 );
	fp2_init_none( &t2 );

	for(int i=0; i<d; i++) {
		fp2_add( &t1, &(ker[i].X), &(ker[i].Z), p );
		fp2_sub( &t2, &(ker[i].X), &(ker[i].Z), p );
		fp2_set( &(ker[i].X), &t1 );
		fp2_set( &(ker[i].Z), &t2 );
	}

	fp2_clear(&t1);
	fp2_clear(&t2);
}

void odd_isogeny( point *S, int d, point ker[d], point *P, mpz_t p ) {
	/* Sets S as phi(P), where phi is the isogeny of kernel ker.
	Uses the formula from Costello and Hisil, to compute the images
	from the point efficiently.

	WATCH OUT ! The kernel points have to be in reshaped form (X+Z, X-Z).
	*/
	fp2 t1, t2;
	fp2_init_none( &t1 );
	fp2_init_none( &t2 );
	point Hat, T, U;
	point_init_none( &Hat ); // (X^ : Z^).
	point_init_none( &T ); // (X' : Z').
	point_init_none( &U ); // (t1 : t2).

	fp2_add( &t1, &( P->X ), &( P->Z ), p );
	fp2_sub( &t2, &( P->X ), &( P->Z ), p );
	fp2_set( &( Hat.X ), &t1 );
	fp2_set( &( Hat.Z ), &t2 );
	// Now (X^ : Z^) = (X+Z : X-Z), with P = (X : Z).

	criss_cross( &T, &(ker[0]), &Hat, p );
	
	for( int i=2; i<=d; i++ ){
		criss_cross( &U, &(ker[i-1]), &Hat, p );
		fp2_mult( &(T.X), &(U.X), &(T.X), p);
		fp2_mult( &(T.Z), &(U.Z), &(T.Z), p);
		// Now (X' : Z') = (t1*X' : t2*Z').
	}

	fp2_mult( &(T.X), &(T.X), &(T.X), p);
	fp2_mult( &(T.Z), &(T.Z), &(T.Z), p);
	fp2_mult( &(T.X), &(P->X), &(T.X), p);
	fp2_mult( &(T.Z), &(P->Z), &(T.Z), p);
	// Now (X' : Z') = (X*(X')² : Z*(Z')²).

	point_set( S, &T );

	point_clear(&T);
	point_clear(&Hat);
	fp2_clear(&t1);
	fp2_clear(&t2);
}


void simultaneous_odd_isogeny( point *R, point *Beta, 
	point *S1, point *S2, point *S3,
	point *G, curve *F, point *P, point *Q, point *QminusP, 
	mpz_t p, int d ) {
	/* Sets R as phi(R), where R is the 'local' image from the 'global' kernel generator, 
	the two torsion point Beta as phi(Beta),
	S1 as phi(P), S2 as phi(Q), S3 as phi(QminusP),
	where phi is an odd isogeny with degree 2d+1 and kernel <G>.

	WATCH OUT ! F is in the form ( A+2C : 4C ) !
	*/
	point_set( S1, P );
	point_set( S2, Q );
	point_set( S3, QminusP );

	point ker[d];
	for( int i = 0; i<d; i++ ) {
	 	point_init_none( &(ker[i]) );
	}
	kernel_point( d, ker, G, F, p );

	kernel_reshape( d, ker, p );
	// Now (Xi : Zi) = (Xi+Zi : Xi-Zi).
	
	odd_isogeny( R, d, ker, R, p );
	odd_isogeny( Beta, d, ker, Beta, p );	
	odd_isogeny( S1, d, ker, S1, p ); 
	odd_isogeny( S2, d, ker, S2, p );
	odd_isogeny( S3, d, ker, S3, p );

	for( int i = 0; i<d; i++ ) {
		point_clear( &(ker[i]) );
	}
}

void simultaneous_odd_isogeny_without_points( point *R, point *Beta, 
	point *G, curve *F,
	mpz_t p, int d ) {
	/* Sets R as phi(R), the two torsion point Beta as phi(Beta),
	where phi is an odd isogeny with degree 2d+1.

	WATCH OUT ! F is in the form ( A+2C : 4C ) !
	*/
	point ker[d];
	for( int i = 0; i<d; i++ ) {
		point_init_none( &(ker[i]) );
	}
	kernel_point( d, ker, G, F, p );

	kernel_reshape( d, ker, p );
	// Now (Xi : Zi) = (Xi+Zi : Xi-Zi).
	
	odd_isogeny( R, d, ker, R, p );
	odd_isogeny( Beta, d, ker, Beta, p );	

	for( int i = 0; i<d; i++ ) {
		point_clear( &(ker[i]) );
	}
}



// Main is only for testing. Uncomment to run the tests on these functions.
/*
int main(int argc, char const *argv[])
{	mpz_t p;
	curve E, F;
	point G;
	point P, S;
	point R, S1, S2, S3, Beta;
	
	// --------------------- TESTS ODD PRIMES ---------------//
	printf("TESTS ODD PRIMES : \n");
	mpz_init(p);
	mpz_set_ui( p, 19 ); 
	// WATCH OUT ! p is 23 later for the second part of testing.
	int d = 2;
	
	curve_init_ui( &E, 0, 0, 1, 0, p );
	curve_init_none( &F );
	curve_Aplus_C( &F, &E, p );

	point_init_ui( &G, 5, 0, 1, 0, p ); // point (5 : 4 : 1).

	printf("Kernel test : \n");
	printf("Should return (5 : 1) and (9 : 1)\n");
	point ker[d];
	for( int i = 0; i<d; i++ ) {
	 	point_init_none( &(ker[i]) );
	}

	kernel_point( d, ker, &G, &F, p );

	for( int i = 0; i<d; i++ ) {
		point_normalize( &(ker[i]), &(ker[i]), p );
	 	point_print( &(ker[i]) );
	}

	printf("Kernel reshape test : \n");
	printf("Should return (6 : 4) and (10 : 8)\n");
	kernel_reshape( d, ker, p );
	for( int i = 0; i<d; i++ ) {
	 	point_print( &(ker[i]) );
	}

	printf("Odd isogeny test : \n ");
	printf("Should return phi( 3 : 1 ) = (18 : 1)\n");
	point_init_ui( &P, 3, 0, 1, 0, p ); // P = (3 : 7 : 1)
	point_init_none( &S );
	odd_isogeny( &S, d, ker, &P, p );
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("Simultaneous odd isogeny test : \n");
	printf("Should return (0 : 0), (16 + 7i : 1), (18 : 1)\n");
	point_init_none( &R );
	point_set( &R, &G );
	point_init_ui( &Beta, 0, 1, 1, 0, p ); // Beta = (0 : 0 : 1)
	point_init_none( &S1 );
	point_init_none( &S2 );
	point_init_none( &S3 );
	
	simultaneous_odd_isogeny(  &R, &Beta, &S1, &S2, &S3,
		&G, &F, &P, &P, &P, p, d );

	point_normalize( &R, &R, p );
	point_normalize( &Beta, &Beta, p ); 
	point_normalize( &S2, &S2, p );
	point_print( &R );
	point_print( &Beta );
	point_print( &S2 );

	printf("Simultaneous odd isogeny without point test : \n");
	printf("Should return (0 : 0), (16 + 7i : 1)\n");
	point_set( &R, &G );
	point_clear(&Beta);
	point_init_ui( &Beta, 0, 1, 1, 0, p ); // Beta = (0 : 0 : 1)
	
	simultaneous_odd_isogeny_without_points(  &R, &Beta,
		&G, &F, p, d );

	point_normalize( &R, &R, p );
	point_normalize( &Beta, &Beta, p ); 
	point_print( &R );
	point_print( &Beta );

	point_clear(&Beta);
	point_clear(&S3);
	point_clear(&S2);
	point_clear(&S1);
	point_clear(&R);
	point_clear(&S);
	point_clear(&P);
	for( int i = 0; i<d; i++ ) {
	 	point_clear( &(ker[i]) );
	}
	point_clear(&G);
	curve_clear(&E);
	curve_clear(&F);
	mpz_clear(p);


	// --------------------- TESTS 2, 3 PRIMES ---------------//
	printf("TESTS 2, 3 PRIMES : \n");
	point Q;
	mpz_init(p);
	mpz_set_ui( p, 23 ); // 2³ * 3 - 1
	curve_init_ui( &E, 0, 0, 1, 0, p );
	curve_init_none( &F );
	curve_Aplus_Aminus( &F, &E, p );
	point_init_ui( &G, 5, 0, 1, 0, p ); // (5 : 10*j : 1) 3 torsion point.
	point_init_ui( &Q, 0, 4, 1, 0, p ); // (4i : 4+4i : 1) 4 torsion point.
	point_init_none( &S );

	printf("[3]G should return (0 : 0)\n");
	xTPL( &S, &G, &F, p );
	point_normalize( &S, &S, p );
	point_print(&S);

	printf("3_iso_curve test : \n");
	printf("Should return (A+ : A-) = (21 : 1) \n");
	fp2 K1, K2, K3;
	fp2_init_none( &K1 );
	fp2_init_none( &K2 );
	fp2_init_none(&K3);
	three_iso_curve( &F, &K1, &K2, &G, p );
	curve_normalize( &F, &F, p );
	curve_print( &F );

	printf("3_iso_eval test : \n");
	printf("Should return phi(Q) = (7+3i : 1) \n");
	three_iso_eval( &S, &K1, &K2, &Q, p );
	point_normalize( &S, &S, p );
	point_print(&S);

	printf("4_iso_curve test : \n");
	printf("Should return 2Q = (22i : 1), (A+2C : 4C) = (12 : 4) \n");
	curve_Aplus_C( &F, &E, p );
	xDBL(&S, &Q, &F, p);
	//xDBL(&S, &S, &F, p);
	point_normalize(&S, &S, p);
	point_print(&S);
	four_iso_curve( &F, &K1, &K2, &K3, &Q, p );
	//curve_normalize( &F, &F, p );
	curve_print( &F );


	fp2_clear(&K1);
	fp2_clear(&K2);
	point_clear(&G);
	point_clear(&Q);
	point_clear(&S);
	curve_clear(&E);
	curve_clear(&F);
	mpz_clear(p);
	

	printf("here !\n");
}
*/
