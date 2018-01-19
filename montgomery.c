/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
* Contains methods for Montgomery efficient arithmetic,
* as well as auxiliary methods to compute j invariant,
* and recovering the curve from three points or from a
* 2 torsion point, with a set of tests.
*
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include "Header.h"


/* --------------------------------------------------------
                MONTGOMERY ARITHMETIC PART
----------------------------------------------------------*/

void xADD( point *A, point *P, point *Q, point *R, mpz_t p ) {
	/* Sets A = P+Q on group F, using Montgomery differential addition.
	R = P-Q.
	*/
	fp2 t1, t2, t3, t4;
	fp2_init_none(&t1);
	fp2_init_none(&t2);
	fp2_init_none(&t3);
	fp2_init_none(&t4);

	fp2_sub( &t1, &(P->X), &(P->Z), p );
	fp2_add( &t2, &(Q->X), &(Q->Z), p );
	fp2_mult( &t3, &t1, &t2, p );

	fp2_add( &t1, &(P->X), &(P->Z), p );
	fp2_sub( &t2, &(Q->X), &(Q->Z), p );
	fp2_mult( &t4, &t1, &t2, p );

	fp2_add( &t1, &t3, &t4, p );
	fp2_mult( &t1, &t1, &t1, p );
	fp2_mult( &(A->X), &t1, &(R->Z), p );

	fp2_sub( &t2, &t3, &t4, p );
	fp2_mult( &t2, &t2, &t2, p );
	fp2_mult( &(A->Z), &t2, &(R->X), p );

	fp2_clear(&t1);
	fp2_clear(&t2);
	fp2_clear(&t3);
	fp2_clear(&t4);

}

void xDBL( point *R, point *P, curve *F, mpz_t p ) {
	/* Sets R = 2P on group F.

	WATCH OUT ! F has to be in the form AplusC = (A+2C : 4C).
	*/
	fp2 t0, t1;
	fp2_init_none(&(t0));
	fp2_init_none(&(t1));
	point T;
	point_init_none(&T);
	
	fp2_sub( &t0, &(P->X), &(P->Z), p );
	fp2_add( &t1, &(P->X), &(P->Z), p );
	fp2_mult( &t0, &t0, &t0, p );
	fp2_mult(&t1, &t1, &t1, p );
	fp2_mult( &(T.Z), &(F->C), &t0, p );
	fp2_mult( &(T.X), &(T.Z), &t1, p );
	fp2_sub( &t1, &t1, &t0, p );
	fp2_mult( &t0, &(F->A), &t1, p );
	fp2_add( &(T.Z), &(T.Z), &t0, p );
	fp2_mult( &(T.Z), &(T.Z), &t1, p );

	point_set( R, &T );
	
	point_clear(&T);
	fp2_clear(&t0);
	fp2_clear(&t1);
}

void xDBLe( point *R, point *P, curve *F, mpz_t p, int e ) {
	/* Sets R = 2^eP on group F.

	WATCH OUT ! F has to be in the form AplusC = (A+2C : 4C).
	*/

	point T;
	point_init_none(&T);

	point_set( &T, P );
	for( int i=0; i < e; i++ ) {
		xDBL( &T, &T, F, p );
	}

	point_set( R, &T );
	point_clear(&T);
}

void xDBLADD( point * twoP, point *PplusQ,
	point *P, point *Q, point *QminusP, curve *F, mpz_t p ) {
	/* Sets twoP as P+P, and PplusQ as P+Q.

	WATCH OUT ! F has to be in the form aplus.
	*/
	fp2 t0, t1, t2;
	fp2_init_none(&t0);
	fp2_init_none(&t1);
	fp2_init_none(&t2);
	point T, S;
	point_init_none(&T);
	point_init_none(&S);

	fp2_add( &t0, &(P->X), &(P->Z), p );
	fp2_sub( &t1, &(P->X), &(P->Z), p );
	fp2_mult( &(T.X), &t0, &t0, p );
	fp2_sub( &t2, &(Q->X), &(Q->Z), p );
	fp2_add( &(S.X), &(Q->X), &(Q->Z), p );
	fp2_mult( &t0, &t0, &t2, p );
	fp2_mult( &(T.Z), &t1, &t1, p );

	fp2_mult( &t1, &t1, &(S.X), p );
	fp2_sub( &t2, &(T.X), &(T.Z), p );
	fp2_mult( &(T.X), &(T.X), &(T.Z), p );
	fp2_mult( &(S.X), &(F->A), &t2, p );
	fp2_sub( &(S.Z), &t0, &t1, p );
	fp2_add( &(T.Z), &(S.X), &(T.Z), p );
	fp2_add( &(S.X), &t0, &t1, p );

	fp2_mult( &(T.Z), &(T.Z), &t2, p );
	fp2_mult( &(S.Z), &(S.Z), &(S.Z), p );
	fp2_mult( &(S.X), &(S.X), &(S.X), p );
	fp2_mult( &(S.Z), &(QminusP->X), &(S.Z), p );
	fp2_mult( &(S.X), &(QminusP->Z), &(S.X), p );

	point_set( twoP, &T );
	point_set( PplusQ, &S );

	point_clear(&T);
	point_clear(&S);
	fp2_clear(&t0);
	fp2_clear(&t1);
	fp2_clear(&t2);
}

void xTPL( point *R, point *P, curve *F, mpz_t p ) {
	/*
	Sets R = 3P on group F.

	WATCH OUT ! F has to be in the form AplusAminus.
	*/
	fp2 t0, t1, t2, t3, t4, t5, t6;
	fp2_init_none(&t0);
	fp2_init_none(&t1);
	fp2_init_none(&t2);
	fp2_init_none(&t3);
	fp2_init_none(&t4);
	fp2_init_none(&t5);
	fp2_init_none(&t6);

	fp2_sub( &t0, &(P->X), &(P->Z), p );
	fp2_mult( &t2, &t0, &t0, p );
	fp2_add( &t1, &(P->X), &(P->Z), p );
	fp2_mult( &t3, &t1, &t1, p );
	fp2_add( &t4, &t1, &t0, p );
	fp2_sub( &t0, &t1, &t0, p );

	fp2_mult( &t1, &t4, &t4, p );
	fp2_sub( &t1, &t1, &t3, p );
	fp2_sub( &t1, &t1, &t2, p );
	fp2_mult( &t5, &t3, &(F->A), p );
	fp2_mult( &t3, &t5, &t3, p );
	fp2_mult( &t6, &t2, &(F->C), p );

	fp2_mult( &t2, &t2, &t6, p );
	fp2_sub( &t3, &t2, &t3, p );
	fp2_sub( &t2, &t5, &t6, p );
	fp2_mult( &t1, &t2, &t1, p );
	fp2_add( &t2, &t3, &t1, p );
	fp2_mult( &t2, &t2, &t2, p );

	fp2_mult( &(R->X), &t2, &t4, p );
	fp2_sub( &t1, &t3, &t1, p );
	fp2_mult( &t1, &t1, &t1, p );
	fp2_mult( &(R->Z), &t1, &t0, p );

	fp2_clear(&t0);
	fp2_clear(&t1);
	fp2_clear(&t2);
	fp2_clear(&t3);
	fp2_clear(&t4);
	fp2_clear(&t5);
	fp2_clear(&t6);
}

void xTPLe( point *R, point *P, curve *F, mpz_t p, int e ) {
	/* Sets R = 3eP on group F.

	WATCH OUT ! F has to be in the form AplusAminus.
	*/
	point T;
	point_init_none(&T);

	point_set( &T, P );
	for( int i=0; i < e; i++ ) {
		xTPL( &T, &T, F, p );
	}

	point_set( R, &T );
	point_clear(&T);

}

void Ladder( point *R, point *P, curve *F, mpz_t m, mpz_t p ) {
	/* 
	Sets R as [m]P, with addition on curve F.

	WATCH OUT ! F has to be in the form ( A+2C : 4C).
	*/
	point R0, R1;
	point_init_none(&R0);
	point_init_none(&R1);
	point_set( &R0, P );
	xDBL( &R1, P, F, p );

	size_t l;
	l = mpz_sizeinbase ( m, 2 );
	for( long i=l-2; i>= 0; i-- ) {
		// Invariant : R1 - R0 = P.
		if( mpz_tstbit( m, i ) ) { 
			xADD( &R0, &R1, &R0, P, p );
			xDBL( &R1, &R1, F, p );
		}
		else {
			xADD( &R1, &R1, &R0, P, p );
			xDBL( &R0, &R0, F, p );
		}
	}
	point_set( R, &R0 );

	point_clear(&R0);
	point_clear(&R1);
}

void Ladder3pt( point *QplusmP, mpz_t m, 
	fp2 *xP, fp2 *xQ, fp2 *xQminusP, curve *F, mpz_t p) {
	/* Sets QplusmP as Q + mP, with xP, xQ and x(P-Q) as inputs,
	using the Montgomery ladder on curve F = (A : C).

	WATCH OUT : the spec suggest that F is in the form aplus, but
	uses Ladder3pt on (A : C) form. Conversion is thus added HERE.
	*/
	// Conversion from (A : C) to (A+2C / 4C : 1)
	curve G;
	curve_init_none( &G );
	curve_aplus( &G, F, p );

	point P0, P1, P2; 
	mpz_t un, zero;
	mpz_inits(un, zero, NULL);
	mpz_set_ui( un, 1);
	fp2 one;
	fp2_init( &one, un, zero, p );

	point_init_fp( &P0, xQ, &one );
	point_init_fp( &P1, xP, &one );
	point_init_fp( &P2, xQminusP, &one );

	size_t l;
	l = mpz_sizeinbase ( m, 2 );
	for( int i=0; i< (l); i++ ) {
		if( mpz_tstbit( m, i ) ) { 
			xDBLADD( &P0, &P1, &P0, &P1, &P2, &G, p );
		}
		else{
			xDBLADD( &P0, &P2, &P0, &P2, &P1, &G, p );
		}
	}
	point_set( QplusmP, &P1 );

	point_clear(&P0);
	point_clear(&P1);
	point_clear(&P2);
	fp2_clear(&one);
	mpz_clears(un, zero, NULL);
	curve_clear(&G);
}

void Ladder3pt_without_conversion( point *QplusmP, mpz_t m, 
	fp2 *xP, fp2 *xQ, fp2 *xQminusP, curve *F, mpz_t p) {
	/* Sets QplusmP as Q + mP, with xP, xQ and x(P-Q) as inputs,
	using the Montgomery ladder on curve F = (A : C).

	WATCH OUT : assert that F is in the form (A+2C / 4C : 1) ! 
	*/
	curve G;
	curve_init_none( &G );
	//curve_aplus( &G, F, p );
	curve_set(&G, F);

	point P0, P1, P2; 
	mpz_t un, zero;
	mpz_inits(un, zero, NULL);
	mpz_set_ui( un, 1);
	fp2 one;
	fp2_init( &one, un, zero, p );

	point_init_fp( &P0, xQ, &one );
	point_init_fp( &P1, xP, &one );
	point_init_fp( &P2, xQminusP, &one );

	size_t l;
	l = mpz_sizeinbase ( m, 2 );
	for( int i=0; i< (l); i++ ) {
		if( mpz_tstbit( m, i ) ) { 
			xDBLADD( &P0, &P1, &P0, &P1, &P2, &G, p );
		}
		else{
			xDBLADD( &P0, &P2, &P0, &P2, &P1, &G, p );
		}
	}
	point_set( QplusmP, &P1 );

	point_clear(&P0);
	point_clear(&P1);
	point_clear(&P2);
	fp2_clear(&one);
	mpz_clears(un, zero, NULL);
	curve_clear(&G);
}

/* --------------------------------------------------------
                		TOOLS PART
----------------------------------------------------------*/

void jInvariant( fp2 *j, curve *E, mpz_t p ) {
	/* Returns j invariant from curve E.
	*/
	fp2 tmp, t0, t1;
	fp2_init_none(&tmp);
	fp2_init_none(&t0);
	fp2_init_none(&t1);

	fp2_mult( &tmp, &(E->A), &(E->A), p );
	fp2_mult( &t1, &(E->C), &(E->C), p );
	fp2_add( &t0, &t1, &t1, p );
	fp2_sub( &t0, &tmp, &t0, p );
	fp2_sub( &t0, &t0, &t1, p );

	fp2_sub( &tmp, &t0, &t1, p );
	fp2_mult( &t1, &t1, &t1, p );
	fp2_mult( &tmp, &tmp, &t1, p );
	fp2_add( &t0, &t0, &t0, p );
	fp2_add( &t0, &t0, &t0, p );

	fp2_mult( &t1, &t0, &t0, p );
	fp2_mult( &t0, &t0, &t1, p );
	fp2_add( &t0, &t0, &t0, p );
	fp2_add( &t0, &t0, &t0, p );
	fp2_multinv( &tmp, &tmp, p );

	fp2_mult( &tmp, &t0, &tmp, p );

	fp2_set( j, &tmp );

	fp2_clear(&tmp);
	fp2_clear(&t0);
	fp2_clear(&t1);
}

void get_A( fp2 *A, fp2 *xP, fp2 *xQ, fp2 *xQminusP, mpz_t p ) {
	/* Sets A as the first parameter of a curve passing through
	xP, xQ and x(Q-P).
	*/
	fp2 t, t0, t1;
	fp2_init_none(&t);
	fp2_init_none(&t0);
	fp2_init_none(&t1);

	mpz_t un, zero;
	mpz_inits(un, zero, NULL);
	mpz_set_ui( un, 1);
	fp2 one;
	fp2_init( &one, un, zero, p );

	fp2_add( &t1, xP, xQ, p );
	fp2_mult( &t0, xP, xQ, p );
	fp2_mult( &t, xQminusP, &t1, p );
	fp2_add( &t, &t, &t0, p );

	fp2_mult( &t0, &t0, xQminusP, p );
	fp2_sub( &t, &t, &one, p );
	fp2_add( &t0, &t0, &t0, p );
	fp2_add( &t1, &t1, xQminusP, p );

	fp2_add( &t0, &t0, &t0, p );
	fp2_mult( &t, &t, &t, p );
	fp2_multinv( &t0, &t0, p );
	fp2_mult( &t, &t, &t0, p );

	fp2_sub( &t, &t, &t1, p );

	fp2_set( A, &t );

	fp2_clear(&t);
	fp2_clear(&t0);
	fp2_clear(&t1);
	fp2_clear(&one);
	mpz_clears(un, zero, NULL);

}

void criss_cross( point *S, point *P, point *Q, mpz_t p ) {
	/* Sets S as ( XP*ZQ + ZP*XQ : XP*ZQ - ZP*XQ ).
	*/
	fp2 t1, t2;
	fp2_init_none( &t1 );
	fp2_init_none( &t2 );

	fp2_mult( &t1, &(P->X), &(Q->Z), p );
	fp2_mult( &t2, &(P->Z), &(Q->X), p );

	fp2_add( &(S->X), &t1, &t2, p );
	fp2_sub( &(S->Z), &t1, &t2, p );

	fp2_clear(&t1);
	fp2_clear(&t2);
}

void curve_from_alpha( curve *F, fp2 *alpha, mpz_t p ) {
	/* Sets F as the curve having (alpha : 1) as two torsion point.

	WATCH OUT ! F will be in the form ( A+2C, 4C ) !
	*/
	fp2 t1, t2, one;
	fp2_init_none( &t1 );
	fp2_init_none( &t2 );
	fp2_init_none( &one );
	mpz_set_ui( one.s0, 1 );

	fp2_sub( &t1, alpha, &one, p );
	fp2_mult( &t1, &t1, &t1, p );

	fp2_add( &t2, alpha, &one, p );
	fp2_mult( &t2, &t2, &t2, p );
	fp2_sub( &t2, &t1, &t2, p );

	fp2_set( &(F->A), &t1 );
	fp2_set( &(F->C), &t2 );
	//F is now in the form (A+2C : 4C).

	fp2_clear(&t1);
	fp2_clear(&t2);
	fp2_clear(&one);
}

void curve_from_Alpha( curve *F, point *Alpha, mpz_t p ) {
	/* Sets F as the curve having Alpha as two torsion point.
	Note : doesn't follows the paper which gives ( A-2C/4, C).

	WATCH OUT ! F will be in the form ( A+2C, 4C) !
	
	*/
	fp2 t1, t2;
	fp2_init_none( &t1 );
	fp2_init_none( &t2 );

	fp2_sub( &t1, &(Alpha->X), &(Alpha->Z), p );
	fp2_mult( &t1, &t1, &t1, p );

	fp2_add( &t2, &(Alpha->X), &(Alpha->Z), p );
	fp2_mult( &t2, &t2, &t2, p );
	fp2_sub( &t2, &t1, &t2, p );

	fp2_set( &(F->A), &t1 );
	fp2_set( &(F->C), &t2 );
	//F is now in the form (A+2C : 4C).

	fp2_clear(&t1);
	fp2_clear(&t2);
}

/* 
HEY ! LISTEN : the prime p has to be 3 mod 4, to avoid
having -1 in Fp. Think about it while testing...
Anyway p = 2^e*3^e' - 1
*/

// Main is only for testing ! 
// Uncomment the following to run the test, or copy/paste it in Main.c. 
//Expected values come from Sage.
/*
int main(int argc, char const *argv[])
{	
	printf("   With curve x^3 + x : \n");
	mpz_t p;
	mpz_init( p );
	mpz_set_ui( p, 7 );
	curve E, F;
	curve_init_ui( &E, 0, 0, 1, 0, p );
	curve_init_none( &F );

	point P, Q, R, S, twoP, PplusQ;
	point_init_ui( &P, 1, 0, 1, 0, p ); // point (1 : 3 : 1) on E.
	point_init_ui( &Q, 5, 0, 1, 0, p ); // point (5 : 2 : 1) on E.
	point_init_ui( &R, 3, 0, 1, 0, p ); // point (3 : 3 : 1) on E. R = P-Q.
	point_init_none( &S);
	point_init_none( &PplusQ );
	point_init_none( &twoP );

	printf("criss_cross P Q = (6 : 3) :\n");
	criss_cross( &S, &P, &Q, p );
	point_print( &S );

	printf("xADD PplusQ = (5 : 1) :\n");
	xADD( &PplusQ, &P, &Q, &R, p );
	point_normalize( &PplusQ, &PplusQ, p );
	point_print( &PplusQ );

	printf("xDBL twoP = (0 : 1) :\n");
	curve_Aplus_C( &F, &E, p );
	xDBL( &twoP, &P, &F, p) ;
	point_normalize( &twoP, &twoP, p );
	point_print( &twoP );

	printf("xDBLe 4Q = (0 : 1) :\n");
	curve_Aplus_C( &F, &E, p );
	xDBLe( &S, &Q, &F, p, 2 );
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("xDBLADD 2P = (0 : 1) PplusQ = (5 : 1) :\n");
	curve_aplus( &F, &E, p );
	xDBLADD( &twoP, &PplusQ, &P, &Q, &R, &F, p );
	point_normalize( &twoP, &twoP, p );
	point_normalize( &PplusQ, &PplusQ, p );
	point_print( &twoP );
	point_print( &PplusQ );

	printf("xTPL 3P = (1 : 1) :\n");
	curve_Aplus_Aminus( &F, &E, p );
	xTPL( &S, &P, &F, p ) ;
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("xTPL 9P = (1 : 1) :\n");
	curve_Aplus_Aminus( &F, &E, p );
	xTPLe( &S, &P, &F, p, 2 ) ;
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("Ladder 7P = (1 : 1) :\n");
	curve_Aplus_C( &F, &E, p );
	Ladder( &S, &P, &F, p, p ) ;
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("Ladder3pt Q+7P = (3: 1) :\n");
	Ladder3pt( &S, p, &(P.X), &(Q.X), &(R.X), &E, p ) ;
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("J_invariant from E = (0 : 1) : 6\n");
	fp2 j;
	fp2_init_none( &j );
	jInvariant( &j, &E, p );
	fp2_print( &j );

	printf("\nget_A from E = (0 : 1) : 0\n");
	fp2 A;
	fp2_init_none( &A );
	get_A( &A, &(P.X), &(Q.X), &(R.X), p );
	fp2_print( &A );

	point Alpha;
	point_init_ui( &Alpha, 0, 1, 1, 0, p );
	point_print(&Alpha);

	printf("\ncurve_from_Alpha  = (5i : 3i) = (2 : 4) \n");
	curve_from_Alpha( &F, &Alpha, p );
	curve_print( &F );	

	printf("\ncurve_from_alpha  = (5i : 3i) = (2 : 4) \n");
	curve_from_alpha( &F, &(Alpha.X), p );
	curve_print( &F );

	fp2_clear(&A);
	fp2_clear(&j);
	point_clear(&P);
	point_clear(&Q);
	point_clear(&R);
	point_clear(&S);
	point_clear(&twoP);
	point_clear(&PplusQ);
	point_clear(&Alpha);
	curve_clear(&E);
	curve_clear(&F);

	printf("\n\n   With curve x^3 + 3x^2 + x : \n");
	mpz_init( p );
	mpz_set_ui( p, 19 );
	curve_init_ui( &E, 3, 0, 1, 0, p );
	curve_init_none( &F );

	point_init_ui( &P, 18, 0, 1, 0, p ); // point (18 : 1 : 1) on E.
	point_init_ui( &Q, 6, 0, 1, 0, p ); // point (6 : 8 : 1) on E.
	// point out of two torsion to test DBL in case of :
	//point_init_ui( &Q, 2, 0, 1, 0, p ); // point (2 : 4i : 1) on E.
	point_init_ui( &R, 8, 0, 1, 0, p ); // point (8 : 16 : 1) on E. R = P-Q.
	point_init_none( &S);
	point_init_none( &PplusQ );
	point_init_none( &twoP );

	printf("criss_cross P Q = (2 : 15) :\n");
	criss_cross( &S, &P, &Q, p );
	point_print( &S );

	printf("xADD PplusQ = (12 : 1) :\n");
	xADD( &PplusQ, &P, &Q, &R, p );
	point_normalize( &PplusQ, &PplusQ, p );
	point_print( &PplusQ );

	printf("xDBL twoP = (0 : 1) :\n");
	curve_Aplus_C( &F, &E, p );
	xDBL( &twoP, &P, &F, p) ;
	point_normalize( &twoP, &twoP, p );
	point_print( &twoP );

	printf("xDBLe 4Q = (0 : 1) :\n");
	curve_Aplus_C( &F, &E, p );
	xDBLe( &S, &Q, &F, p, 2 );
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("xDBLe 8Q = (0 : 0) :\n");
	curve_Aplus_C( &F, &E, p );
	xDBLe( &S, &Q, &F, p, 23 );
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("xDBLADD 2P = (0 : 1) PplusQ = (12 : 1) :\n");
	curve_aplus( &F, &E, p );
	xDBLADD( &twoP, &PplusQ, &P, &Q, &R, &F, p );
	point_normalize( &twoP, &twoP, p );
	point_normalize( &PplusQ, &PplusQ, p );
	point_print( &twoP );
	point_print( &PplusQ );

	printf("xTPL 3Q = (16 : 1) :\n");
	curve_Aplus_Aminus( &F, &E, p );
	xTPL( &S, &Q, &F, p ) ;
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("xTPL 9Q = (6 : 1) :\n");
	curve_Aplus_Aminus( &F, &E, p );
	xTPLe( &S, &Q, &F, p, 2 ) ;
	point_normalize( &S, &S, p );
	point_print( &S );

	// 19 !
	printf("Ladder 7Q = (16 : 1) :\n");
	curve_Aplus_C( &F, &E, p );
	Ladder( &S, &Q, &F, p, p ) ;
	point_normalize( &S, &S, p );
	point_print( &S );

	printf("Ladder3pt Q+14P = (16 : 1) :\n");
	mpz_t t;
	mpz_init(t);
	mpz_set_ui(t, 14);
	Ladder3pt( &S, t, &(Q.X), &(P.X), &(R.X), &E, p ) ;
	point_normalize( &S, &S, p );
	point_print( &S );
	mpz_clear(t);

	printf("J_invariant from E = (0 : 1) : 5\n");
	fp2_init_none( &j );
	jInvariant( &j, &E, p );
	fp2_print( &j );

	printf("\nget_A from E = (3 : 0) : 3\n");
	fp2_init_none( &A );
	get_A( &A, &(P.X), &(Q.X), &(R.X), p );
	fp2_print( &A );

	fp2_clear(&A);
	fp2_clear(&j);
	point_clear(&P);
	point_clear(&Q);
	point_clear(&R);
	point_clear(&S);
	point_clear(&twoP);
	point_clear(&PplusQ);
	curve_clear(&E);
	curve_clear(&F);

}
*/

