#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include "Header.h"

typedef struct stack_elem {
	int e;
	point P;
} stack_elem;

typedef struct stack {
	int capacity;
	int top;
	stack_elem* array;
} stack;

/* ---------------------------------------------------------
                    STACK_ELEM PART
-----------------------------------------------------------*/

void stack_elem_init( stack_elem* elem, int e, point* P ) {
	elem->e = e;
	point_init_none( &(elem->P) );
	point_set( &(elem->P), P ) ;
}

void stack_elem_init_none( stack_elem* elem ) {
	elem->e = 0;
	point_init_none( &(elem->P) );
}

void stack_elem_print( stack_elem* elem ) {
	printf("%d\n", elem->e);
	point_print( &(elem->P) );
}

void stack_elem_set( stack_elem* target, stack_elem* source ){
	target->e = source->e;
	point_set( &(target->P), &(source->P) );
}

void stack_elem_quick_set( stack_elem* target, int e, point* P ){
	target->e = e;
	point_set( &(target->P), P );
}

void stack_elem_clear( stack_elem* elem ) {
	point_clear( &(elem->P) );
}

/* ---------------------------------------------------------
                    	STACK PART
-----------------------------------------------------------*/

stack* stack_init( int capacity ) {
	stack* S = (struct stack*) malloc(sizeof(struct stack));
	S->capacity = capacity;
	S->top = -1;
	S->array = (stack_elem*) malloc(S->capacity * sizeof(stack_elem));
	return S;
}

void stack_clear( stack* S ) {
	free((S->array));
	free(S);
}

int stack_is_full( struct stack* S ) {   
	return (S->top == S->capacity - 1); 
}

int stack_is_empty( struct stack* S ){   
	return (S->top == -1);  
}

void stack_push( stack* S, stack_elem* elem ) {
	if( stack_is_full(S) ) {
		printf("Stack is full !!");
	}
	else {
		(S->top) ++;
		(S->array)[S->top] = *elem;
	}
}

void stack_pop(stack_elem* elem, stack* S)
{
    if (stack_is_empty(S))
        printf("Stack is empty !!");
    else{
    	stack_elem_set( elem, &( S->array[ S->top ] ) );
    	(S->top) --;
    }
}

void stack_set( stack* S_target, stack* S_source ) {
	S_target->capacity = S_source->capacity;
	S_target->top = S_source->top;
	S_target->array = S_source->array;
}

/* ---------------------------------------------------------
                    	STATEGY PART
-----------------------------------------------------------*/

void three_e_iso_with_strategy( curve *F, point *phiP1, point *phiP2, point *phiP3,
	curve *E, point *G, point *P1, point *P2, point *P3, 
	int e3, mpz_t p, int strategy[e3-1] ) {
	/* Sets F as the (3^e3)-isogenous curve, ie F = E/<G>,
	where S has order 3^e3 in E.
	WATCH OUT E and F are in the form AplusAmoins !
	*/
	// Some preliminaries.
	fp2 K1, K2;
	fp2_init_none(&K1);
	fp2_init_none(&K2);

	curve_set(F, E);

	point_set(phiP1, P1);
	point_set(phiP2, P2);
	point_set(phiP3, P3);

	//1
	stack* S;
	stack* S_prime;
	S = stack_init(200);

	//2
	stack_elem elem, elem_bis;
	stack_elem_init( &elem, e3, G );
	stack_push( S, &elem );

	//3
	int i = 1;

	//4
	while( !( stack_is_empty(S) ) ) {
		//5
		stack_pop( &elem, S ); // h = elem.e
		//6
		if( elem.e == 1 ) {
			//7
			three_iso_curve( F, &K1, &K2, &(elem.P), p );
			// 8
			S_prime = stack_init(200);
			// 9
			while( !( stack_is_empty(S) ) ) {
				// 10
				stack_pop( &elem, S );
				// 11
				three_iso_eval( &(elem.P), &K1, &K2, &(elem.P), p );
				// 12
				stack_elem_quick_set( &elem_bis, ( (elem.e) - 1), &(elem.P) );
				stack_push( S_prime, &elem_bis );
			}
			// 13
			stack_set( S, S_prime );

			//14 - 15
			three_iso_eval( phiP1, &K1, &K2, phiP1, p );
			three_iso_eval( phiP2, &K1, &K2, phiP2, p );
			three_iso_eval( phiP3, &K1, &K2, phiP3, p );

			// Release stack S'
			stack_clear(S_prime);
		}
		// 16
		else if( (strategy[i] > 0) && (strategy[i] < elem.e) ) { // h = elem.e
			// 17
			stack_push( S, &elem );
			// 18
			xTPL( &(elem.P), &(elem.P), F, p );
			// 19
			stack_elem_quick_set( &elem_bis, (elem.e - strategy[i]), &(elem.P) );
			stack_push( S, &elem_bis );
			// 20
			i++;
		}
		else {
			printf("Sorry, this is an invalid strategy... :'(\n");
		}
	} 

	fp2_clear(&K2);
	fp2_clear(&K1);
	stack_elem_clear(&elem);
	stack_elem_clear(&elem_bis);
	stack_clear(S);
}

/*
int main(int argc, char const *argv[])
{	
	mpz_t p;
	mpz_init(p);
	mpz_set_ui(p, 7);
	point P;
	point_init_ui( &P, 1, 0, 3, 0, p); // P = (1 : 3)

	stack_elem E, F;
	stack_elem_init( &E, 4, &P );
	stack_elem_print( &E );
	stack_elem_init_none( &F );

	stack* S;
	S = stack_init(4);
	printf("%d\n", S->top);
	printf( "%d\n", stack_is_empty(S) );
	printf( "%d\n", stack_is_full(S) );
	stack_push( S, &E );
	stack_pop( &F, S );
	stack_elem_print( &F );
	printf( "%d\n", stack_is_empty(S) );


	stack_clear(S);
	stack_elem_clear( &E );
	stack_elem_clear( &F );
	point_clear(&P);
	mpz_clear(p);
	printf("here !\n");
}
*/