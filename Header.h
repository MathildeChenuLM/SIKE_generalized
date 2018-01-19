/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*
* Contains declaration and short description for each method
* in the poject, list in bottom/up approach.
*
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

#ifndef HEADER_H
#define HEADER_H

/* ---------------------------------------------------------
                       Fp PART
-----------------------------------------------------------*/

typedef struct fp2 {
   mpz_t s0;
   mpz_t s1;
} fp2;

void fp_inv( mpz_t *res, mpz_t a, mpz_t p );
	/* Sets res with the inverse of a mod p.
	*/
void fp2_init_none( fp2 *res );
	/* Initializes res with value 0 for every fields.
	*/
void fp2_init( fp2 *res, mpz_t x, mpz_t y, mpz_t p );
	/* Initializes res with fields s0 = x and s1 = y.
	*/
void fp2_clear( fp2 *res );
	/* Clears the fields of res.
	*/
void fp2_print( fp2 *my_fp2 );
	/* Prints my_fp2 fields. 
	*/
void fp2_add( fp2 *res, fp2 *a, fp2 *b, mpz_t p);
	/* Adds a and b and writes the result in res.
	*/
void fp2_sub( fp2 *res, fp2 *a, fp2 *b, mpz_t p);
	/* Substractes a and b and writes the result in res.
	*/
void fp2_mult( fp2 *res, fp2 *a, fp2 *b, mpz_t p );
	/* Multiplies a and b and writes the result in res.
	*/
void fp2_addinv( fp2 *res, fp2 *a, mpz_t p );
	/* Fills res with the additive inverse of a. 
	Note : The gmp_mod function always returns a positive value.
	*/
void fp2_multinv( fp2 *res, fp2 *a, mpz_t p );
	/* Fills res with the multiplicative inverse of a, ie res*a = 1.
	*/
int fp2_are_equal( fp2 *a, fp2 *b);
	/* Returns 1 if a and b are equal, and 0 otherwise.
	*/
int fp2_is_zero( fp2 *a, mpz_t p );
	/* Returns 1 if a is (0 : 0), and 0 otherwise.
	*/
void fp2_set( fp2 *res, fp2 *a );
	/* Sets fields from res with fields from a.
	*/

/* ---------------------------------------------------------
                    CURVE & POINT PART
-----------------------------------------------------------*/

typedef struct curve {
   fp2 A;
   fp2 C;
} curve;

typedef struct point {
   fp2 X;
   fp2 Z;
} point;

void curve_init_none( curve *E );
	/* Initializes a curve E with parameters
	(A, C) = ( 0 + i*0, 0 + i*0).
	*/
void curve_init_fp( curve *E, fp2 *A, fp2 *C );
	/* Initializes a curve P with parameters (A, C);
	*/
void curve_init_mpz( curve *E, mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t p);
	/* Initializes a curve E with 
	parameters ( a + ib, c + id).
	*/
void curve_init_ui( curve *E, int a, int b, int c, int d, mpz_t p );
	/* Initializes a curve E with 
	parameters ( a + ib, c + id), with a,b,c,d of type int.
	*/
void curve_init_normalize( curve *E, fp2 *A, mpz_t p );
	/* Initializes a curve E with 
	parameters ( A, 1 ).
	*/
void curve_init_fp2_int( curve *E, fp2 *A, int a, mpz_t p );
	/* Initializes a curve E with 
	parameters ( A, a ), with A of type fp2 and a of type int.
	*/
void curve_set( curve *E, curve *F );
	/* Sets curve E as curve F.
	*/ 
void curve_clear( curve *E );
	/* Clears fields of curve E.
	*/
void curve_print( curve *E );
	/* Prints E on two lines, 
	respectively E.A and E.C with fp representation.
	*/
void curve_Aplus_C( curve *F, curve *E, mpz_t p );
	/* Sets the curve F such that F = (A+2C, 4C),
	where E = (A, C).
	*/
void curve_Aplus_Aminus( curve *F, curve *E, mpz_t p );
	/* Sets the curve F such that F = (A+2C, A-2C),
	where E = (A, C).
	*/
void curve_aplus( curve * F, curve *E, mpz_t p );
	/* Sets the curve F such that F = (A+2C / 4C, 1),
	where E = (A, C);
	*/
void curve_normalize( curve *F, curve *E, mpz_t p );
	/* Sets F as ( A/C, 1 ), where E = ( A, C ).
	*/

void point_init_none( point *P );
	/* Initializes a point P as ( 0 + i*0, 0 + i*0).
	*/
void point_init_fp( point *P, fp2 *X, fp2 *Z );
	/* Initializes a point P as (X, Z).
	*/
void point_init_mpz( point *P, mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t p);
	/* Initializes a point P as ( a + ib, c + id).
	*/
void point_init_ui( point *P, int a, int b, int c, int d, mpz_t p );
	/* Initializes a point P with 
	parameters ( a + ib, c + id), with a,b,c,d of type int.
	*/
void point_init_normalize( point *P, fp2 * xP, mpz_t p );
	/* Initializes a point P as ( xP, 1 ).
	*/
void point_set( point *R, point *P );
	/* Sets point R as point P.
	*/ 
void point_clear( point *P );
	/* Clears fields of point P.
	*/
void point_print( point *P );
	/*Prints P on two lines, 
	respectively P->X and P->y with fp representation.
	*/
int points_are_equal( point *P, point *Q );
	/* Returns 1 if points are equal, and 0 otherwise.
	Note : Might need some normalization before.
	*/
void point_normalize( point *R, point *P, mpz_t p );
	/* Sets R as ( X/Z, 1 ), where P = ( X, Z ), with Z non zero.
	WATCH OUT ! If Z is 0, sets R as (0, 0).
	*/
void point_normalize_X( fp2 *X, point *P, mpz_t p );
	/* Sets X as  X/Z, where P = ( X, Z ), with Z non zero.
	WATCH OUT ! If Z is 0, sets R as (0, 0).
	*/

/* ---------------------------------------------------------
                     MONTGOMERY PART
-----------------------------------------------------------*/

void xADD( point *A, 
	point *P, point *Q, point *R, mpz_t p );
	/* Sets A = P+Q on group F, using Montgomery differential addition.
	R = P-Q.
	*/
void xDBL( point *R, 
	point *P, curve *F, mpz_t p );
	/* Sets R = 2P on group F.
	WATCH OUT ! F has to be in the form AplusC = (A+2C : 4C).
	*/
void xDBLe( point *R, 
	point *P, curve *F, mpz_t p, int e );
	/* Sets R = 2^eP on group F.
	WATCH OUT ! F has to be in the form AplusC = (A+2C : 4C).
	*/
void xDBLADD( point * twoP, point *PplusQ,
	point *P, point *Q, point *QminusP, curve *F, mpz_t p );
	/* Sets twoP as P+P, and PplusQ as P+Q.
	WATCH OUT ! F has to be in the form aplus.
	*/
void xTPL( point *R, 
	point *P, curve *F, mpz_t p );
	/*
	Sets R = 3P on group F.
	WATCH OUT ! F has to be in the form AplusAminus.
	*/
void xTPLe( point *R, 
	point *P, curve *F, mpz_t p, int e );
	/* Sets R = 3eP on group F.
	WATCH OUT ! F has to be in the form AplusAminus.
	*/
void Ladder( point *R, 
	point *P, curve *F, mpz_t m, mpz_t p );
	/* 
	Sets R as [m]P, with addition on curve F.
	WATCH OUT ! F has to be in the form ( A+2C : 4C).
	*/
void Ladder3pt( point *QplusmP, 
	mpz_t m, fp2 *xP, fp2 *xQ, fp2 *xQminusP, curve *F, mpz_t p);
	/* Sets QplusmP as Q + mP, with xP, xQ and x(P-Q) as inputs,
	using the Montgomery ladder on curve F = (A : C).
	WATCH OUT : the spec suggest that F is in the form aplus, but
	uses Ladder3pt on (A : C) form. Conversion is thus added HERE.
	*/
void Ladder3pt_without_conversion( point *QplusmP, mpz_t m, 
	fp2 *xP, fp2 *xQ, fp2 *xQminusP, curve *F, mpz_t p);
	/* Sets QplusmP as Q + mP, with xP, xQ and x(P-Q) as inputs,
	using the Montgomery ladder on curve F = (A : C).
	WATCH OUT : assert that F is in the form (A+2C : 4C) ! 
	*/
void jInvariant( fp2 *j, 
	curve *E, mpz_t p );
	/* Returns j invariant from curve E.
	*/
void get_A( fp2 *A, 
	fp2 *xP, fp2 *xQ, fp2 *xQminusP, mpz_t p );
	/* Sets A as the first parameter of a curve passing through
	xP, xQ and x(Q-P).
	*/
void criss_cross( point *S, 
	point *P, point *Q, mpz_t p );
	/* Sets S as ( XP*ZQ + ZP*XQ : XP*ZQ - ZP*XQ ).
	*/
void curve_from_alpha( curve *F, 
	fp2 *alpha, mpz_t p );
	/* Sets F as the curve having (alpha : 1) as two torsion point.
	WATCH OUT ! F will be in the form ( A+2C, 4C ) !
	NOTE : needs alpha as an element from fp2. 
	*/
void curve_from_Alpha( curve *F, 
	point *Alpha, mpz_t p );
	/* Sets F as the curve having Alpha as two torsion point.
	WATCH OUT ! F will be in the form ( A+2C, 4C) !
	Note : doesn't follows the paper which gives ( A-2C/4, C)
	*/

/* ---------------------------------------------------------
                    STACK_ELEM & STACK PART
-----------------------------------------------------------*/

typedef struct stack_elem {
	int e;
	point P;
} stack_elem;

typedef struct stack {
	int capacity;
	int top;
	stack_elem* array;
} stack;

void stack_elem_init( stack_elem* elem, int e, point* P );
void stack_elem_init_none( stack_elem* elem );
void stack_elem_print( stack_elem* elem );
void stack_elem_set( stack_elem* target, stack_elem* source );
void stack_elem_quick_set( stack_elem* target, int e, point* P );
void stack_elem_clear( stack_elem* elem );

stack* stack_init( int capacity );
void stack_clear( stack* S );
int stack_is_full( struct stack* S );
int stack_is_empty( struct stack* S );
void stack_push( stack* S, stack_elem* elem );
void stack_pop(stack_elem* elem, stack* S);
void stack_set( stack* S_target, stack* S_source );


/* ---------------------------------------------------------
                     ISOGENY PART
-----------------------------------------------------------*/
void four_iso_curve( curve *F, fp2 *K1, fp2 *K2, fp2 *K3, 
	point *P4, mpz_t p );
	/* Sets F as the four-isogenous curve, ie F = E/<P4>,
	where P4 has order 4.
	WATCH OUT : F will be in the form AplusC = (A+2C : 4C)
	*/
void four_iso_eval( point* Q, 
	fp2 *K1, fp2 *K2, fp2 *K3, point *P, mpz_t p );
	/* Sets Q as phi(P), where phi is a four-isogeny,
	described by K1, K2 and K3.
	*/
void three_iso_curve( curve *F, fp2 *K1, fp2 *K2, 
	point* P, mpz_t p );
	/* Sets F as the three-isogenous curve, ie F = E/<P>,
	where P has order 3.
	WATCH OUT : curve will be in the form (A+ : A-) = (A+2C : A-2C)
	*/
void three_iso_eval( point *phiP, 
	fp2 *K1, fp2* K2, point *P, mpz_t p );
	/* Sets phiP as phi(P), where phi is a three-isogeny,
	described by K1 and K2.
	Note : Does not follows exactly the specification, error on line 3.
	*/
void two_e_iso( curve *F, point *phiP1, point *phiP2, point *phiP3,
	curve *E, point *S, point *P1, point *P2, point *P3, 
	int e2, mpz_t p );
	/* Sets F as the (2^e2)-isogenous curve, ie F = E/<S>,
	where S has order 2^e2 in E.
	Also sets phiPi as phi(Pi).
	WATCH OUT ! E and F are in AplusC form.
	*/
void three_e_iso( curve *F, point *phiP1, point *phiP2, point *phiP3,
	curve *E, point *S, point *P1, point *P2, point *P3, 
	int e3, mpz_t p );
	/* Sets F as the (3^e3)-isogenous curve, ie F = E/<S>,
	where S has order 3^e3 in E.
	WATCH OUT E and F are in the form AplusAmoins !
	*/
void three_e_iso_with_strategy( curve *F, point *phiP1, point *phiP2, point *phiP3,
	curve *E, point *G, point *P1, point *P2, point *P3, 
	int e3, mpz_t p, int strategy[e3-1] );
	/* Sets F as the (3^e3)-isogenous curve, ie F = E/<S>,
	where S has order 3^e3 in E. Uses strategy.
	WATCH OUT E and F are in the form AplusAmoins !
	*/

void kernel_point( int d, point ker[d], 
	point *G, curve *F, mpz_t p );
	/*
	Computes the d points from the subgroup of the curve F 
	generated by G, with G a [2d+1] torsion point. Sets ker as <G>.
	WATCH OUT ! F has to be in the form (A+2C : 4C) !
	*/
void kernel_reshape( int d, point ker[d], mpz_t p );
	/*
	Reshapes the points in the kernel as (X+Z : X-Z).
	The kernel is then ready for odd_isogeny.
	*/
void odd_isogeny( point *S, 
	int d, point ker[d], point *P, mpz_t p );
	/* Sets S as phi(P), where phi is the isogeny of kernel ker.
	Uses the formula from Costello and Hisil, to compute the images
	from the point efficiently.
	WATCH OUT ! The kernel points have to be in reshaped form (X+Z, X-Z).
	*/
void simultaneous_odd_isogeny( point *R, point *Beta, 
	point *S1, point *S2, point *S3,
	point *G, curve *F, point *P, point *Q, point *QminusP, 
	mpz_t p, int d );
	/* Sets R as phi(R), where R is the 'local' image from the 'global' kernel generator, 
	sets the two torsion point Beta as phi(Beta),
	sets S1 as phi(P), S2 as phi(Q), S3 as phi(QminusP),
	where phi is an odd isogeny with degree 2d+1 and kernel <G>.
	WATCH OUT ! F is in the form ( A+2C : 4C ) !
	*/
void simultaneous_odd_isogeny_without_points( point *R, point *Beta, 
	point *G, curve *F,
	mpz_t p, int d );
	/* Sets R as phi(R), the two torsion point Beta as phi(Beta),
	where phi is an odd isogeny with degree 2d+1.
	WATCH OUT ! F is in the form ( A+2C : 4C ) !
	*/

/* ---------------------------------------------------------
        			PK, SK, PARAMETERS
-----------------------------------------------------------*/

typedef struct pk {
	fp2 beta;
	fp2 x1;
   	fp2 x2;
   	fp2 x3;
} pk;

typedef struct sk {
	mpz_t secret;
} sk;

typedef struct parameters {
   	int eA;
   	int eB;
   	int pA;
   	int pB;
   	mpz_t p;
   	fp2 alpha;
   	fp2 xPA;
   	fp2 xQA;
   	fp2 xRA;
   	fp2 xPB;
   	fp2 xQB;
   	fp2 xRB;
} parameters;

void pk_init_none( pk *PK );
	/* Initializes public key PK with value 0 for every fields.
	*/
void pk_clear( pk *PK );
	/* Clears the fields of public key PK.
	*/
void pk_print( pk *PK );
	/* Prints the fields of public key PK.
	*/

void sk_init( sk *SK, mpz_t secret );
	/* Initializes secret key SK with secret value.
	*/
void sk_clear( sk *SK );
	/* Clears the field of secret key SK.
	*/

void init_parameters( parameters *param );
	/* Initializes public key PK with value 0 for every fields.
	*/
void set_prime( parameters *param, 
	mpz_t p, int pA, int pB, int eA, int eB );
	/* Sets the public parameters in param as given inputs.
	*/
void set_points( parameters *param, 
	mpz_t alpha0, mpz_t alpha1, 
	mpz_t xPA0, mpz_t xPA1, mpz_t xQA0, mpz_t xQA1, mpz_t xRA0, mpz_t xRA1,
	mpz_t xPB0, mpz_t xPB1, mpz_t xQB0, mpz_t xQB1, mpz_t xRB0, mpz_t xRB1 );
	/* Sets the public parameters in param as given inputs.
	*/
void clear_parameters( parameters *param );
	/* Clears the fields of public parameters param.
	*/

void isogen_two( pk *PK2, sk *SK2, parameters *param );
	/* Key generation for the 2-isogeny side.
	Sets the public key according to the secret key and the parameters.
	*/
void isogen_three( pk *PK3, sk *SK3, parameters *param );
	/* Key generation for the 3-isogeny side.
	Sets the public key according to the secret key and the parameters.
	*/
void isoex_two( fp2 *j, sk *SK2, pk *PK3, parameters *param );
	/* Key exchange for the 2-isogeny side.
	Sets jA as the shared secret.
	*/
void isoex_three( fp2 *j, sk *SK3, pk *PK2, parameters *param );
	/* Key exchange for the 3-isogeny side.
	Sets jB as the shared secret.
	*/

void isogen_Alice( pk *PKA, sk *SKA, parameters *param );
	/* Key generation for any odd isogeny.
	Sets the public key according to the secret key and the parameters.
	*/
void isogen_Bob( pk *PKB, sk *SKB, parameters *param );
	/* Key generation for any odd isogeny.
	Sets the public key according to the secret key and the parameters.
	*/
void isoex_Alice( fp2 *jA, pk *PKB, sk *SKA, parameters *param );
	/* Key exchange for any odd isogeny.
	Sets jA as the shared secret.
	*/
void isoex_Bob( fp2 *jB, pk *PKA, sk *SKB, parameters *param );
	/* Key exchange for any odd isogeny.
	Sets jB as the shared secret.
	*/

/* ---------------------------------------------------------
        					MAIN
-----------------------------------------------------------*/

void key_gen_Alice( pk *PKA, sk *SKA, parameters *param, int regular );
	/* Key_generation for Alice. Given a secret key SKA and public parameters param,
	sets PKA to the associated public key.
	Uses algorithms for 2/3 only if regular is set to 1.
	*/
void key_gen_Bob( pk *PKB, sk *SKB, parameters *param, int regular );
	/* Key_generation for Bob. Given a secret key SKB and public parameters param,
	sets PKB to the associated public key.
	Uses algorithms for 2/3 only if regular is set to 1.
	*/
void key_exchange_Alice( fp2 *jA, sk *SKA, pk *PKB, parameters *param, int regular );
	/* Key exchange for Alice. Given Bob's public key PKB, Alice's own secret key SKA and
	public parameters param, sets jA to the shared j-invariant jA.
	Uses algorithms for 2/3 only if regular is set to 1.
	*/
void key_exchange_Bob( fp2 *jB, sk *SKB, pk *PKA, parameters *param, int regular );
	/* Key exchange for Bob. Given Alice's public key PKA, Bob's own secret key SKB and
	public parameters param, sets jB to the shared j-invariant jB.
	Uses algorithms for 2/3 only if regular is set to 1.
	*/
void verification( fp2 *jA, fp2 *jB );
	/* Checks if jA and jB are equal over Fp2.
	*/

#endif