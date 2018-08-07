/*   This file is part of the program Snaffour.
 *
 *   Copyright (C) 2018 by Marc Culler, Nathan Dunfield, Matthias Görner
 *   and others.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Project homepage: https://bitbucket.org/t3m/snaffour/
 *   Author homepage: https://marc-culler.info
 *   Author homepage: http://dunfield.info
 *   Author homepage: http://www.unhyperbolic.org/
 */

#ifndef F_FOUR_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

/** A Term16_t is an array of 16 chars representing exponents.
 *
 * Assigning the vector_size attribute enables gcc to use the 128-bit SSE
 * MMX registers to do arithmetic operations on all 16 bytes in a single
 * instruction.  This is meant to optimize operations such as multiplying
 * or dividing Terms.
 */

typedef char Term16_t __attribute__ ((vector_size (16)));

/** Terms
 *
 * We use two Term16_t vectors to allow up to 32 variables.  For computing
 * Ptolemy varieties, 16 variables would not be enough.  The code uses a
 * for loop to process the two 16 byte vectors, but those loops will be
 * unrolled by the optimizer.
 *
 * Wrapping the array of two Term16_t types in a struct makes Cython easier
 * to deal with.
 */
typedef struct Term_s {
  Term16_t degree[2];
} Term_t;

void Term_print(Term_t *t, int rank);
bool Term_equals(Term_t* t, Term_t* s);
int  Term_total_degree(Term_t *t, int rank);
bool Term_divides(Term_t* t, Term_t* s);
bool Term_divide(Term_t *t, Term_t *s, Term_t *answer);
void Term_multiply(Term_t *t, Term_t *s, Term_t *answer);
void Term_lcm(Term_t *t, Term_t *s, Term_t *answer);
void Term_gcd(Term_t *t, Term_t *s, Term_t *answer);
int  Term_revlex_diff(Term_t *t, Term_t *s, int rank);
long Term_hash(Term_t *t);
bool Term_merge(Term_t* s, Term_t* t, int s_size, int t_size,
		Term_t** answer, int* answer_size, int rank);

/** Coefficients in a Polynomial.
 *
 * The coefficient also stores a column index, which is used in the echelon
 * reduction to speed up comparisons of Terms.  The index should be set to
 * INDEX_UNSET when a coefficient is created and is reset to INDEX_UNSET after
 * computing the reduced echelon form.
 */

/* Negative values of the column index are used as flags.*/
#define INDEX_UNSET  -1
#define DUPLICATE  -2

typedef struct coeff_s {
  int column_index;
  int value;
} coeff_t;

typedef struct monomial_s {
  int total_degree;
  Term_t* term;
  coeff_t* coefficient;
} monomial_t;

typedef struct monomial_array_s {
  int size;
  monomial_t* monomials;
} monomial_array_t;

int inverse_mod(int p, int x);

/** Polynomials
 *
 * Polynomials come in two flavors.  The standard flavor has two arrays,
 * one containing terms and one containing coefficients.  Using two arrays
 * saves memory since a term must be stored in memory which is aligned to
 * 16 bytes while a coefficient only occupies 8 bytes.
 *
 * The second "compact" flavor saves even more memory by using an external table
 * of terms.  The term_order element of the coefficient is an index into the
 * external table, and the internal term pointer is NULL.  The external table
 * can be shared among several polynomials. This is done when computing the
 * echelon form of a "matrix" whose "rows" are Polynomials.  This saves both
 * space and time.  The basic row operation consists primarily of copying
 * coefficients from one row to another, based on a comparison of the associated
 * terms. Terms can be compared by comparing their indices, which is much faster
 * than comparing them as vectors, even with MMX instructions.  And copying 8
 * bytes is much faster than copying 40 bytes.  The compact flavor also saves
 * time by storing the values of the coefficients in their Montgomery
 * representation, which allows fast reduction of products mod P.  This is
 * important because the other major component of the basic row operation
 * row_i -> row_i + a*row_j is that every non-zero coefficient in row_j must
 * be multiplied by a.
 *
 * A Polynomial in compact form should have its terms element set to NULL and
 * each coefficient should have a non-negative column_index.  A Polynomial in
 * standard form should have its table element set to NULL and each column_index
 * should be set to INDEX_UNSET.
 *
 * The rank indicates the number of variables, i.e. the rank of the parent ring.
 * Loops which deal with the exponents as separate bytes use this to avoid
 * iterating through the unused exponents which, incidentally, are expected to
 * all be zero.
 *
 * The array terms of Term_t types, or the external table must be maintained in
 * sorted order by descending grevlex.  This means that adding or subtracting
 * two Polynomials is done by interleaving the lists of terms and coefficients,
 * or just the coefficients in the compact case, and occasionally combining two
 * coefficients in the relatively rare case where a term appears in both
 * Polynomials.  Finding the coefficient of a Term in a Polynomial is done by
 * bisection, which assumes that the term arrays be ordered.
 *
 * The zero polynomial typically has num_terms=0 and both pointers equal to NULL.
 * However, any polynomial with num_terms=0 is equal to 0.  The meaning of
 * num_terms is algebraic -- it does not indicate how much memory has been
 * allocated for the Polynomial.
 *
 * In our Cython implementaton, a Polynomial object holds a Polynomial_t within
 * the object data structure.  So the task of allocating and deallocating
 * Polynomials is almost always handled by Python.  The task of allocating and
 * deallocating the arrays of terms and coefficients is handled by this C
 * library.  (The Cython __cinit__ and __dealloc__ methods call Poly_alloc
 * and Poly_free.)
 */

typedef struct Polynomial_s {
  int num_terms;          /* How many terms in this Polynomial. */
  int max_size;           /* How many terms and coefficients the allocated memory will hold. */
  int rank;               /* The number of variables in the parent polynomial ring. */
  coeff_t* coefficients;
  Term_t* terms;
  Term_t* table;
} Polynomial_t;

bool Poly_alloc(Polynomial_t *P, int size, int rank);
void Poly_free(Polynomial_t* P);
void Poly_print(Polynomial_t* P, int rank);
void Poly_copy(Polynomial_t* src, Polynomial_t* dest);
void Poly_new_term(Polynomial_t* P, Term_t* term, coeff_t coefficient, int rank);
bool Poly_equals(Polynomial_t* P, Polynomial_t *Q);
bool Poly_add(Polynomial_t* P, Polynomial_t* Q, Polynomial_t* answer, int prime, int rank);
bool Poly_sub(Polynomial_t* P, Polynomial_t* Q, Polynomial_t* answer, int prime, int rank);
int  Poly_coeff(Polynomial_t* P, Term_t* t, int rank);
bool Poly_make_row(Polynomial_t* P, Term_t* t, Polynomial_t* answer, int prime, int rank);
void Poly_make_monic(Polynomial_t* P, int prime, int rank, int mu, int R_cubed);
bool Poly_echelon(Polynomial_t** P, Polynomial_t *answer, int num_rows, int *num_columns,
                  int prime, int rank);
bool Poly_times_term(Polynomial_t *P, Term_t* t, Polynomial_t* answer, int prime, int rank);
bool Poly_times_int(Polynomial_t* P, int a, Polynomial_t* answer, int prime, int rank);
void Poly_sort(Polynomial_t* P, int num_polys, bool increasing);
bool Poly_terms(Polynomial_t* P, int num_polys, Term_t** answer, int* answer_size, int rank);
int Poly_column_index(Polynomial_t* P, Term_t* t, int rank);


/* Use Python's memory allocation */
#include "pymem.h"

#define malloc PyMem_Malloc
#define realloc PyMem_Realloc
#define free PyMem_Free

#define F_FOUR_H
#endif
