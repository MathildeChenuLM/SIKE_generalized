/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
* Contains definitions and methods to set parameters and 
* secret key, to generate public key, and compute shared 
* secret.
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
 					PUBLIC KEY
----------------------------------------------------------*/

void pk_init_none( pk *PK ) {
	/* Initializes public key PK with value 0 for every fields.
	*/
	fp2_init_none( &(PK->beta) );
	fp2_init_none( &(PK->x1) );
	fp2_init_none( &(PK->x2) );
	fp2_init_none( &(PK->x3) );
}

void pk_clear( pk *PK ) {
	/* Clears the fields of public key PK.
	*/
	fp2_clear( &(PK->beta) );
	fp2_clear( &(PK->x1) );
	fp2_clear( &(PK->x2) );
	fp2_clear( &(PK->x3) );
}

void pk_print( pk *PK ) {
	printf("2 torsion point (optional) :\n");
	fp2_print( &(PK->beta) );
	printf("x1 :\n");
	fp2_print( &(PK->x1) );
	printf("x2 :\n");
	fp2_print( &(PK->x2) );
	printf("x3 :\n");
	fp2_print( &(PK->x3) );
	printf("\n\n");
}

/* --------------------------------------------------------
                    SECRET KEY
----------------------------------------------------------*/

void sk_init( sk *SK, mpz_t secret ) {
	/* Initializes secret key SK with secret value.
	*/
	mpz_init( SK->secret );
	mpz_set( SK->secret, secret );
}

void sk_clear( sk *SK ) {
	mpz_clear( SK->secret );
}

/* --------------------------------------------------------
            	PARAMETERS SETTING
----------------------------------------------------------*/

void init_parameters( parameters *param ) {
	/* Initializes public parameters param with 
	value 0 for every fields.
	*/
	mpz_init( param->p );

	fp2_init_none( &(param->alpha) );
	fp2_init_none( &(param->xPA) );
	fp2_init_none( &(param->xQA) );
	fp2_init_none( &(param->xRA) );
	fp2_init_none( &(param->xPB) );
	fp2_init_none( &(param->xQB) );
	fp2_init_none( &(param->xRB) );
}

void set_prime( parameters *param, mpz_t p, int pA, int pB, int eA, int eB ) {
	mpz_set( param->p, p );
	param->pA = pA ;
	param->pB = pB ;
	param->eA = eA;
	param->eB = eB;
}

void set_points( parameters *param, 
	mpz_t alpha0, mpz_t alpha1, 
	mpz_t xPA0, mpz_t xPA1, mpz_t xQA0, mpz_t xQA1, mpz_t xRA0, mpz_t xRA1,
	mpz_t xPB0, mpz_t xPB1, mpz_t xQB0, mpz_t xQB1, mpz_t xRB0, mpz_t xRB1 ) {
	/*
	*/
	//fp2_init( &(param->alpha), alpha0, alpha1, param->p );
	fp2_init( &(param->alpha), alpha0, alpha1, param->p );

	fp2_init( &(param->xPA), xPA0, xPA1, param->p );
	fp2_init( &(param->xQA), xQA0, xQA1, param->p );
	fp2_init( &(param->xRA), xRA0, xRA1, param->p );

	fp2_init( &(param->xPB), xPB0, xPB1, param->p );
	fp2_init( &(param->xQB), xQB0, xQB1, param->p );
	fp2_init( &(param->xRB), xRB0, xRB1, param->p );
}

void clear_parameters( parameters *param ) {
	/*
	*/
	mpz_clear( param->p );

	fp2_clear( &(param->alpha) );
	fp2_clear( &(param->xPA) );
	fp2_clear( &(param->xQA) );
	fp2_clear( &(param->xRA) );
	fp2_clear( &(param->xPB) );
	fp2_clear( &(param->xQB) );
	fp2_clear( &(param->xRB) );
}

/* --------------------------------------------------------
               	KEY GENERATION REGULAR
----------------------------------------------------------*/

void isogen_two( pk *PK2, sk *SK2, parameters *param ) {
	/* Key generation for the "two-isogeny side".
	Sets the public key according to the secret key and the parameters.
	*/
	curve E, F;
	curve_init_ui( &E, 0, 0, 1, 0, param->p ); // (A : C) = (0 : 1)
	curve_init_ui( &F, 1, 0, 2, 0, param->p ); // (A+2C : 4C)

	point P1, P2, P3;
	point_init_normalize( &P1, &(param->xPB), param->p );
	point_init_normalize( &P2, &(param->xQB), param->p );
	point_init_normalize( &P3, &(param->xRB), param->p );

	point S;
	point_init_none(&S);
	Ladder3pt( &S, 
	 	(SK2->secret), &(param->xPA), &(param->xQA), &(param->xRA), &E, (param->p) );
	// Now S = [secret]P2 + Q2.

	two_e_iso( &F, &P1, &P2, &P3, 
		&F, &S, &P1, &P2, &P3, param->eA, param->p );
	// Now F = Ea, the image of E by the 2^e2 isogeny with kernel <S>.

	point_normalize( &P1, &P1, param->p );
	point_normalize( &P2, &P2, param->p );
	point_normalize( &P3, &P3, param->p );

	fp2_set( &(PK2->x1), &(P1.X) );
	fp2_set( &(PK2->x2), &(P2.X) );
	fp2_set( &(PK2->x3), &(P3.X) );
	// F can be described with abscisses x1, x2, and x3.

	point_clear(&S);
	point_clear(&P1);
	point_clear(&P2);
	point_clear(&P3);
	curve_clear(&E);
	curve_clear(&F);
}

void isogen_three( pk *PK3, sk *SK3, parameters *param ) {
	/* Key generation for the "three-isogeny side".
	Sets the public key according to the secret key and the parameters.
	*/
	curve E, F;
	curve_init_ui( &E, 0, 0, 1, 0, param->p );
	curve_init_ui( &F, 2, 0, -2, 0, param->p );

	point P1, P2, P3;
	point_init_normalize( &P1, &(param->xPA), param->p );
	point_init_normalize( &P2, &(param->xQA), param->p );
	point_init_normalize( &P3, &(param->xRA), param->p );

	point S;
	point_init_none(&S);
	Ladder3pt( &S, 
		(SK3->secret), &(param->xPB), &(param->xQB), &(param->xRB), &E, (param->p) );
	
	three_e_iso( &F, &P1, &P2, &P3, 
		&F, &S, &P1, &P2, &P3, param->eB, param->p );

	point_normalize( &P1, &P1, param->p );
	point_normalize( &P2, &P2, param->p );
	point_normalize( &P3, &P3, param->p );

	fp2_set( &(PK3->x1), &(P1.X) );
	fp2_set( &(PK3->x2), &(P2.X) );
	fp2_set( &(PK3->x3), &(P3.X) );

	point_clear(&S);
	point_clear(&P1);
	point_clear(&P2);
	point_clear(&P3);
	curve_clear(&E);
	curve_clear(&F);
}

/* --------------------------------------------------------
               	KEY EXCHANGE REGULAR
----------------------------------------------------------*/

void isoex_two( fp2 *j, sk *SK2, pk *PK3, parameters *param ) {
	/* Key exchange for the "two-isogeny side".
	Sets j as the j invariant of the common curve.
	*/
	curve E;
	curve_init_none(&E);
	fp2 A;
	fp2_init_none(&A);

	get_A( &A, &(PK3->x1), &(PK3->x2), &(PK3->x3), param->p );
	curve_init_normalize( &E, &A, param->p );

	point S;
	point_init_none(&S);
	Ladder3pt( &S, 
		SK2->secret, &(PK3->x1), &(PK3->x2), &(PK3->x3), &E, param->p );

	// (Aplus, C) = (A+2, 4)
	curve F;
	curve_init_none(&F);
	curve_Aplus_C( &F, &E, param->p );

	point T;
	point_init_none(&T);
	two_e_iso( &F, &T, &T, &T, &F, &S, &T, &T, &T, param->eA, param->p );

	// (A, C) = (4Aplus - 2C, C)
	fp2_set( &(E.C), &(F.C) );
	fp2 tmp;
	fp2_init_none(&tmp);
	fp2_add( &tmp, &(F.C), &(F.C), param->p );
	fp2_add( &(E.A), &(F.A), &(F.A), param->p );
	fp2_add( &(E.A), &(E.A), &(E.A), param->p );
	fp2_sub( &(E.A), &(E.A), &tmp, param->p );

	jInvariant( j, &E, param->p );

	fp2_clear(&tmp);
	point_clear(&T);
	curve_clear(&F);
	point_clear(&S);
	fp2_clear(&A);
	curve_clear(&E);
}

void isoex_three( fp2 *j, sk *SK3, pk *PK2, parameters *param ) {
	/* Key exchange for the "three-isogeny side".
	Sets j as the j invariant of the common curve.
	*/
	curve E;
	curve_init_none(&E);
	fp2 A;
	fp2_init_none(&A);

	get_A( &A, &(PK2->x1), &(PK2->x2), &(PK2->x3), param->p );
	curve_init_normalize( &E, &A, param->p );

	point S;
	point_init_none(&S);
	Ladder3pt( &S, 
		SK3->secret, &(PK2->x1), &(PK2->x2), &(PK2->x3), &E, param->p );

	// (Aplus, C) = (A+2, A-2)
	curve F;
	curve_init_none(&F);
	curve_Aplus_Aminus( &F, &E, param->p );

	point T;
	point_init_none(&T);
	three_e_iso( &F, &T, &T, &T, &F, &S, &T, &T, &T, param->eB, param->p );

	// (E.A, E.C) = (2*( F.A + F.C ), ( F.A - F.C ) )
	fp2 tmp;
	fp2_init_none(&tmp);
	fp2_sub( &(E.C), &(F.A), &(F.C), param->p );
	fp2_add( &(E.A), &(F.A), &(F.C), param->p );
	fp2_add( &(E.A), &(E.A), &(E.A), param->p );

	jInvariant( j, &E, param->p );

	fp2_clear(&tmp);
	point_clear(&T);
	curve_clear(&F);
	point_clear(&S);
	fp2_clear(&A);
	curve_clear(&E);
}

/* --------------------------------------------------------
               	KEY GENERATION GENERAL
----------------------------------------------------------*/

void isogen_Alice( pk *PKA, sk *SKA, parameters *param ) {
	/* Key generation for any odd isogeny.
	Sets the public key according to the secret key and the parameters.

	Variable written with a capital letter are points. Variable written
	with a small letter are elements from Fp2.
	*/
	point Beta, P1, P2, P3; // Beta is a two torsion point.
	point_init_normalize( &Beta, &(param->alpha), param->p );
	point_init_normalize( &P1, &(param->xPB), param->p );
	point_init_normalize( &P2, &(param->xQB), param->p );
	point_init_normalize( &P3, &(param->xRB), param->p );

	curve E;
	curve_init_none( &E );
	curve_from_alpha( &E, &(param->alpha), param->p );
	curve_normalize( &E, &E, param->p );

	point R;
	point_init_none( &R );
	// Note that E is already in the form aPlus = (A+2C / 4C : 1).
	Ladder3pt_without_conversion( &R, 
		(SKA->secret), &(param->xPA), &(param->xQA), &(param->xRA), &E, (param->p) );
	// Now R = [secret]PA + QA

	curve F;
	curve_init_none( &F );
	point S;
	point_init_none( &S );
	mpz_t q;
	mpz_init(q);
	
	for( int i = ( (param->eA) - 1); i>= 0; i-- ) {
		curve_from_Alpha( &F, &Beta, param->p ); 
	 	// WATCH OUT ! F will be in the form ( A+2C, 4C ) !	
	 	mpz_ui_pow_ui( q, (param->pA), i );

	 	Ladder( &S, &R, &F, q, param->p );

	 	simultaneous_odd_isogeny( &R, &Beta, 
	 		&P1, &P2, &P3, 
	 		&S, &F, &P1, &P2, &P3, 
	 		param->p, ( (param->pA)-1)/2 ); 
	}

	point_normalize_X( &(PKA->beta), &Beta, param->p );
	point_normalize_X( &(PKA->x1), &P1, param->p );
	point_normalize_X( &(PKA->x2), &P2, param->p );
	point_normalize_X( &(PKA->x3), &P3, param->p );

	mpz_clear(q);
	point_clear(&S);
	curve_clear(&F);
	point_clear(&R);
	curve_clear(&E);
	point_clear(&Beta);
	point_clear(&P1);
	point_clear(&P2);
	point_clear(&P3);

}

void isogen_Bob( pk *PKB, sk *SKB, parameters *param ) {
	/* Key generation for any odd isogeny.
	Sets the public key according to the secret key and the parameters.

	Variable written with a capital letter are points. Variable written
	with a small letter are elements from Fp2.
	*/
	point Beta, P1, P2, P3; // Beta is a two torsion point.
	point_init_normalize( &Beta, &(param->alpha), param->p );
	point_init_normalize( &P1, &(param->xPA), param->p );
	point_init_normalize( &P2, &(param->xQA), param->p );
	point_init_normalize( &P3, &(param->xRA), param->p );

	curve E;
	curve_init_none( &E );
	curve_from_alpha( &E, &(param->alpha), param->p );
	curve_normalize( &E, &E, param->p );
	
	point R;
	point_init_none( &R );
	Ladder3pt_without_conversion( &R, 
	 	(SKB->secret), &(param->xPB), &(param->xQB), &(param->xRB), &E, (param->p) );
	// Now R = PB + [secret]QB.

	curve F;
	curve_init_none( &F );
	point S;
	point_init_none( &S );
	mpz_t q;
	mpz_init(q);
	
	for( int i = ( (param->eB) - 1 ); i>= 0; i-- ) {
		curve_from_Alpha( &F, &Beta, param->p ); 
		// WATCH OUT ! F will be in the form ( A+2C, 4C ) !

		mpz_ui_pow_ui( q, (param->pB), i );

		Ladder( &S, &R, &F, q, param->p );

		simultaneous_odd_isogeny( &R, &Beta, 
			&P1, &P2, &P3, 
			&S, &F, &P1, &P2, &P3, 
			param->p, ( (param->pB)-1)/2 ); 
	}
	point_normalize_X( &(PKB->beta), &Beta, param->p );
	point_normalize_X( &(PKB->x1), &P1, param->p );
	point_normalize_X( &(PKB->x2), &P2, param->p );
	point_normalize_X( &(PKB->x3), &P3, param->p );

	mpz_clear(q);
	point_clear(&S);
	curve_clear(&F);
	point_clear(&R);
	point_clear(&Beta);
	point_clear(&P1);
	point_clear(&P2);
	point_clear(&P3);

}

/* --------------------------------------------------------
               	KEY EXCHANGE GENERAL
----------------------------------------------------------*/

void isoex_Alice( fp2 *jA, pk *PKB, sk *SKA, parameters *param ) {
	/* Key exchange for any odd isogeny.

	Variable written with a capital letter are points. Variable written
	with a small letter are elements from Fp2.
	*/
	point Beta; 
	point_init_normalize( &Beta, &(PKB->beta), param->p );
	// Beta is a two torsion point on EB given by Alice's public key.

	curve EB;
	curve_init_none( &EB );
	curve_from_Alpha( &EB, &Beta, param->p );
	curve_normalize( &EB, &EB, param->p );

	point R;
	point_init_none( &R );
	Ladder3pt_without_conversion( &R, 
	 	(SKA->secret), &(PKB->x1), &(PKB->x2), &(PKB->x3), &EB, (param->p) );
	// Now R = phiB(PA) + [secret]phiB(QA).

	curve F;
	curve_init_none( &F );
	point S;
	point_init_none( &S );
	mpz_t q;
	mpz_init(q);
	
	for( int i = ( (param->eA) - 1 ); i>= 0; i-- ){
		mpz_ui_pow_ui( q, (param->pA), i );
		curve_from_Alpha( &F, &Beta, param->p ); 
		Ladder( &S, &R, &F, q, param->p );

		simultaneous_odd_isogeny_without_points( &R, &Beta, 
			&S, &F, 
			param->p, ( (param->pA)-1)/2 ); 
	}
	
	curve EAB;
	curve_init_none( &EAB );
	curve_from_Alpha( &EAB, &Beta, param->p );
	jInvariant( jA, &EAB, param->p );

	curve_clear(&EAB);
	mpz_clear(q);
	point_clear(&S);
	curve_clear(&F);
	point_clear(&R);
	point_clear(&Beta);

}

void isoex_Bob( fp2 *jB, pk *PKA, sk *SKB, parameters *param ) {
	/* Key exchange for any odd isogeny.

	Variable written with a capital letter are points. Variable written
	with a small letter are elements from Fp2.
	*/
	point Beta; 
	point_init_normalize( &Beta, &(PKA->beta), param->p );
	// Beta is a two torsion point on EA given by Alice's public key.

	curve EA;
	curve_init_none( &EA );
	curve_from_Alpha( &EA, &Beta, param->p );
	curve_normalize( &EA, &EA, param->p );

	point R;
	point_init_none( &R );
	Ladder3pt_without_conversion( &R, 
	 	(SKB->secret), &(PKA->x1), &(PKA->x2), &(PKA->x3), &EA, (param->p) );
	// Now R = phiA(PB) + [secret]phiA(QB).

	curve F;
	curve_init_none( &F );
	point S;
	point_init_none( &S );
	mpz_t q;
	mpz_init(q);
	
	for( int i = ( (param->eB) - 1 ); i>= 0; i-- ){
		mpz_ui_pow_ui( q, (param->pB), i );
		curve_from_Alpha( &F, &Beta, param->p ); 
		Ladder( &S, &R, &F, q, param->p );

		simultaneous_odd_isogeny_without_points( &R, &Beta, 
			&S, &F,
			param->p, ( (param->pB)-1)/2 ); 
	}

	curve EAB;
	curve_init_none( &EAB );
	curve_from_Alpha( &EAB, &Beta, param->p );
	jInvariant( jB, &EAB, param->p );

	curve_clear(&EAB);
	mpz_clear(q);
	point_clear(&S);
	curve_clear(&F);
	point_clear(&R);
	point_clear(&Beta);

}






