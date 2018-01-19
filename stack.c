#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include "Header.h"


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

void stack_elem_set( stack_elem* target, stack_elem* source ) {
	target->e = source->e;
	point_set( &(target->P), &(source->P) );
}

void stack_elem_quick_set( stack_elem* target, int e, point* P ) {
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

int stack_is_empty( struct stack* S ) {   
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

void stack_pop(stack_elem* elem, stack* S) {
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