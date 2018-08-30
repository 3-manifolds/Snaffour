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

/** Coefficients.
 *
 * Coefficients are used in both Polynomials and Rows.
 * The coefficient stores a column index, which is used in the echelon
 * reduction to speed up comparisons of Terms.  The index should be set to
 * NO_INDEX if the coefficient appears in a Polynomial.
 */

/* Negative values of the column index are used as flags.*/
#define NO_INDEX  -1
#define DUPLICATE  -2

typedef struct coeff_s {
  int column_index;  /* Used when converting Polynomials to Rows. */
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

/** Polynomials
 *
 * A Polynomials manages two arrays, one containing terms and one containing
 * coefficients.  Using two arrays saves memory since a term must be stored in
 * memory which is aligned to 16 bytes while a coefficient only occupies 8
 * bytes.
 *
 * The rank indicates the number of variables, i.e. the rank of the parent ring.
 * Loops which deal with the exponents as separate bytes use this to avoid
 * iterating through the unused exponents which, incidentally, are expected to
 * all be zero.
 *
 * The array terms of Term_t types must be maintained in sorted order by
 * descending grevlex.  This means that adding or subtracting two Polynomials is
 * done by interleaving the lists of terms and coefficients, or just the
 * coefficients in the compact case, and occasionally combining two coefficients
 * in the relatively rare case where a term appears in both Polynomials.
 * Finding the coefficient of a Term in a Polynomial is done by bisection, which
 * assumes that the term arrays be ordered.
 *
 * The zero polynomial typically has num_terms=0 and both pointers equal to
 * NULL.  However, any polynomial with num_terms=0 is equal to 0.  The meaning
 * of num_terms is algebraic -- it does not indicate how much memory has been
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
  int num_terms;   /* How many terms are used by this Polynomial. */
  int max_size;    /* How many terms (and coefficients) have been allocated. */
  int rank;        /* The rank of the parent polynomial ring. */
  coeff_t* coefficients;
  Term_t* terms;
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
void Poly_make_monic(Polynomial_t* P, int prime, int rank);
bool Poly_echelon(Polynomial_t** P, Polynomial_t *answer, int num_rows, int *num_columns,
                  int prime, int rank);
bool Poly_times_term(Polynomial_t *P, Term_t* t, Polynomial_t* answer, int prime, int rank);
bool Poly_times_int(Polynomial_t* P, int a, Polynomial_t* answer, int prime, int rank);
void Poly_sort(Polynomial_t* P, int num_polys, bool increasing);
bool Poly_terms(Polynomial_t* P, int num_polys, Term_t** answer, int* answer_size, int rank);

/** Rows in a Polynomial Matrix
 *
 * A Row is a compact representation of a Polynomial used when
 * computing echelon forms. A row saves memory by using a sorted
 * external table of terms.  The term_order element of the coefficient
 * is an index into the external table.  The external table is shared
 * among all rows in a matrix.  This saves both space and time. The
 * basic row operation consists primarily of copying coefficients from
 * one row to another, based on a comparison of the associated
 * terms. Terms can be compared by comparing their indices, which is
 * much faster than comparing them as vectors, even with MMX
 * instructions.  And copying 8 bytes is much faster than copying 40
 * bytes.  The row also saves time by storing the values of the
 * coefficients in their Montgomery representation, which allows fast
 * reduction of products mod P.  This is important because the other
 * major component of the basic row operation row_i -> row_i + a*row_j
 * is that every non-zero coefficient in row_j must be multiplied by
 * a.
 */

//typedef coeff_t row_coeff_t;
typedef int64_t row_coeff_t;

typedef struct Row_s {
  int num_terms;   /* How many coefficients are used by this Row. */
  int max_size;    /* How many coefficients have been allocated. */
  row_coeff_t* coefficients;
  Term_t* term_table;
} Row_t;

/** Arithmetic
 *
 * When computing echelon forms over Fp, we use the Montgomery representation
 * of a conjugacy class mod p.  Given a class X, its Montgomery representative
 * M(X) is an element of [0, p) such that M(X) is congruent to RX mod p, where
 * the Montgomery radix R will always be 2^31 for us.
 */

int inverse_mod(int p, int x);

/** Constants
 *
 * An MConstant struct holds the constants needed for computing with
 * Montgomery representations, and for converting between standard
 * and Montgomery representatives.
 */

typedef struct MConstants_s {
  int prime;
  int mu;
  int R_mod_p;
  int R_squared;
  int R_cubed;
} MConstants_t;

#define M_RADIX (1 << 31)
#define M_RADIX64 ((int64_t)(1) << 31)

/* Reducing mod R is equivalent to anding with this constant. */

#define MOD_R 0x7fffffff

/**
 * The distributive law implies M(X) + M(Y) = M(X+Y).  However, M(X)*M(Y) is
 * not congruent to M(X*Y) mod p.  In fact M(X)*M(Y) = R*M(X*Y) mod p.  So
 * computing the Montgomery representative of a product requires being able to
 * divide by R mod p.  The point of choosing R to be a power of 2 is that this
 * can be done without using hardware division.
 *
 * Given X in [0, p^2) the following macro computes an element of [0, 2p)
 * representing X/R in Fp. It requires a precomputed constant mu in [0, p)
 * such that mu represents -1/R in Fp.  (Hardware division can by used to
 * compute mu, since it is only done once per echelon form.)  The macro uses 2
 * multiplies, 2 AND operations 1 addition and a shift.  At most one extra
 * subtraction is needed to compute the standard representative of X/R in Fp.
 *
 * NOTE: the arguments for the macro must be 64 bits wide to avoid overflow!!!
 *
 */

#define M_REDUCE(X, p, mu) ((X + (((X & MOD_R)*mu) & MOD_R)*p) >> 31)

#define F_FOUR_H
#endif
