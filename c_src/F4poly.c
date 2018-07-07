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

#include "F4.h"
#include <stdlib.h>

/** Struct to hold the result of an extended gcd computation.
 *
 * For inputs x and y, this struct is meant to hold d, a and b where
 * d is a greatest common divisor (with arbitrary sign) and d = ax + by. 
 */

typedef struct gcd_s {
  int d;
  int a;
  int b;
} gcd_t;

/** Extended Greatest Common Divisor
 *
 * Compute d = gcd(x,y), a and b such that a*x + b*y = d.
 * The sign of d may or may not be positive.  If x and y are
 * both 0, the result is undefined.
 */

static void x_gcd(int x, int y, gcd_t* answer) {
  static int d, a, b;
  int q, r, temp;
  if (y == 0) {
    if (x == 0) {
      return;
    } else {
      d = x; a = 1; b = 0;
    }
  } else {
    q = x / y; r = x % y;
    x_gcd(y, r, NULL);
    temp = b; b = a - q*b; a = temp;
    if (answer != NULL) {
      answer->d = d; answer->a = a; answer->b = b;
    }
  }
}

/** \page modP Arithmetic mod p.
 *
 * All computations return an int value in the interval [0, p).
 * Multiplications must be done with 64-bit arithmetic to avoid overflow.
 */


/** Compute the inverse of a number modulo p
 */

static inline int inverse_mod(int p, int x) {
  int a;
  gcd_t xgcd = {.d=0, .a=0, .b=0};
  x_gcd(x, p, &xgcd);
  a = xgcd.d < 0 ? (-xgcd.a) % p : xgcd.a % p;
  return a < 0 ? p + a : a;
}

/** Multiply two numbers mod [
 *
 * Be careful about overflow!
 */

static inline int multiply_mod(int p, int x, int y) {
  int64_t P = (int64_t)p, X = (int64_t)x, Y = (int64_t)y, answer;
  if (x == 1) {
    answer = Y;
  } else if (x == -1) {
    answer = -Y;
  } else {
    answer = X*Y;
  }
  answer = answer % P;
  return (answer < 0 ? (int)(P + answer) : (int)answer);
}

/** Compute x + ay mod p
 *
 * This is the scalar version of a row operation.
 */

static inline int x_plus_ay_mod(int p, int x, int a, int y) {
  int64_t P = p, A = a, X = x, Y = y, answer;
  if (a == 1) {
    answer = X + Y;
  } else if (a == p - 1) {
    answer = X - Y;
  } else {
    answer = X + A*Y;
  }
  answer = answer % P;
  if (answer < 0) {
    answer += P;
  }
  return (int)answer;
}

/** \page modP Polynomials
 */

/** Allocate memory for a Polynomial with specified number of terms.
 *
 * Sets num_terms to 0, but allocates memory for terms and coefficients.
 * Note that allocation and deallocation of Polynomial_s structs is not
 * handled by this library.  That shoud be done by Python.
 */

bool Poly_alloc(Polynomial_t* P, size_t size, int rank) {
  Poly_free(P);
  P->num_terms = 0;
  P->rank = rank;
  P->terms = malloc(sizeof(Term_t)*size);
  if (P->terms == NULL) {
    return false;
  }
  P->coefficients = malloc(sizeof(coeff_t)*size);
  if (P->coefficients == NULL) {
    free(P->terms);
    return false;
  }
  return true;
}

/** Free the terms and coefficients of a Polynomial, if not NULL.
 *
 * The pointers are set to NULL after freeing the arrays and num_terms is
 * set to 0.
 */

void Poly_free(Polynomial_t* P) {
  if (P->terms != NULL) {
    free(P->terms);
  }
  if (P->coefficients != NULL) {
    free(P->coefficients);
  }
  P->num_terms = 0;
  P->rank = 0;
  P->terms = NULL;
  P->coefficients = NULL;
}

/** Print a Polynomial in a crude way, without variable names.
 */

void Poly_print(Polynomial_t* P, int rank) {
  if (P->num_terms == 0) {
    printf("0\n");
  } else {
    for (int i=0; i<P->num_terms; i++) {
      printf("%d*", P->coefficients[i].value);
      Term_print(&P->terms[i], rank);
    }
    printf("\n");
  }
}

/** Initialize a Polynomial from arrays of terms and coefficients.
 */

void Poly_init(Polynomial_t* P, size_t size, Term_t* terms, coeff_t* coefficients, int rank) {
  int i;
  P->num_terms = size;
  P->rank = rank;
  for (i=0; i < size; i++) {
    P->terms[i] = terms[i];
    P->coefficients[i] = coefficients[i];
  }
}

/** Static function to compute P + a*Q for an int a.
 *
 * This is the core of the row operation used to reduce a matrix to echelon form
 * and also handles addition and subtraction (by taking a=1 or a=-1).
 *
 * NOTE: P and Q must point to different polynomials for this to work.
 */

static bool Poly_p_plus_aq(Polynomial_t* P, int a, Polynomial_t* Q, Polynomial_t* answer,
                  int prime, int rank) {
  int size = P->num_terms + Q->num_terms, p = 0, q = 0, N = 0, p_td, q_td;
  Term_t term;
  coeff_t p_coeff, q_coeff;
  if (! Poly_alloc(answer, size, rank)) {
    return false;
  }
  while (p < P->num_terms && q < Q->num_terms) {
    p_td = P->coefficients[p].total_degree;
    q_td = Q->coefficients[q].total_degree;
    int td_cmp = p_td - q_td;
    int revlex_cmp = td_cmp == 0 ? Term_revlex_diff(Q->terms + q, P->terms + p, rank) : 0;
    if (td_cmp > 0 || revlex_cmp > 0) {
      /* deg P > deg Q */
      answer->terms[N] = P->terms[p];
      answer->coefficients[N++] = P->coefficients[p++];
    } else if (td_cmp < 0 || revlex_cmp < 0) {
      /* deg P < deg Q */
      answer->terms[N] = Q->terms[q];
      q_coeff = Q->coefficients[q++];
      q_coeff.value = multiply_mod(prime, a, q_coeff.value);
      answer->coefficients[N++] = q_coeff;
    } else {
      /* deg P == deg Q */
      term = P->terms[p];
      p_coeff = P->coefficients[p++];
      q_coeff = Q->coefficients[q++];
      p_coeff.value = x_plus_ay_mod(prime, p_coeff.value, a, q_coeff.value);
      if (p_coeff.value != 0) {
        answer->terms[N] = term;
        answer->coefficients[N++] = p_coeff;
      }
    }
  }
  /* At most one of these two loops will be non-trivial. */
  for (; q < Q->num_terms; q++) {
    answer->terms[N] = Q->terms[q];
    q_coeff = Q->coefficients[q];
    q_coeff.value = multiply_mod(prime, a, q_coeff.value);
    answer->coefficients[N++] = q_coeff;
  }
  for (; p < P->num_terms; p++) {
    answer->terms[N] = P->terms[p];
    answer->coefficients[N++] = P->coefficients[p];
  }
  answer->num_terms = N;
  answer->rank = rank;
  /*
   * It's not clear if this is a good idea. It often won't save much memory and
   * could conceivably cause fragmentation or waste time.  But it does free
   * everything when the answer is 0, at least on linux.
   */
  answer->terms = realloc((void*)answer->terms, sizeof(Term_t)*N);
  answer->coefficients = realloc((void*)answer->coefficients, sizeof(coeff_t)*N);
  if (N != 0 && (answer->terms == NULL || answer->coefficients == NULL)) {
      Poly_free(answer);
      return false;
    }
  return true;
}

/** Add Polynomials P and Q and store the result in answer.
 *
 * The work is done by Poly_p_plus_aq, but we need to deal with the special case
 * where P and Q are the same Polynomial.
 */

bool Poly_add(Polynomial_t* P, Polynomial_t* Q, Polynomial_t* answer, int prime, int rank) {
  if (P == Q) {
    int N = 0;
    coeff_t p_coeff;
    if (! Poly_alloc(answer, P->num_terms, rank)) {
      return false;
    }
    for (N=0; N < P->num_terms; N++) {
      answer->terms[N] = P->terms[N];
      p_coeff = P->coefficients[N];
      p_coeff.value = x_plus_ay_mod(prime, p_coeff.value, 1, p_coeff.value);
      answer->coefficients[N] = p_coeff;
    }
    answer->num_terms = N;
    answer->rank = rank;
    return true;
  }
  return Poly_p_plus_aq(P, 1, Q, answer, prime, rank);
}

/** Subtract Polynomials P and Q and store the result P - Q in answer.
 *
 * The work is done by Poly_p_plus_aq, but we need to deal with the special case
 * where P and Q are the same Polynomial (by returning a 0 Polynomial).
 */

bool Poly_sub(Polynomial_t* P, Polynomial_t* Q, Polynomial_t* answer, int prime, int rank) {
  if (P == Q) {
    Poly_free(answer);
    return true;
  }
  return Poly_p_plus_aq(P, prime - 1, Q, answer, prime, rank);
}

/** Use bisection to find the index of the given term in P->terms for a Polynomial P.
 *
 * Return false if the term is not in P->terms.
 */

static bool find_index(Polynomial_t* P, Term_t* t, int t_td, int rank,
		      int bottom, int top, int* index) {
  int middle;
  int64_t td;
  if (top - bottom == 1) {
    if (Term_equals(t, P->terms + bottom)) {
	*index = bottom;
	return true;
      } else {
	return false;
      }
  }
  middle = (top + bottom) >> 1;
  td = P->coefficients[middle].total_degree;
  if (td < t_td || (td == t_td && Term_revlex_diff(P->terms + middle, t, rank) > 0)) {
    return find_index(P, t, t_td, rank, bottom, middle, index);
  } else {
    return find_index(P, t, t_td, rank, middle, top, index);
  }
}

/** Return the coefficient of the given term in P.
 *
 * If P->terms does not contain the term, return 0.  This uses find_index to
 * find the coefficient or determine that it does not exist.
 */

int Poly_coeff(Polynomial_t* P, Term_t* t, int rank) {
  int index, t_td = Term_total_degree(t, rank);
  if (P->num_terms == 0) {
    return 0;
  }
  if (find_index(P, t, t_td, rank, 0, P->num_terms, &index)) {
    return P->coefficients[index].value;
  }
  return 0;
}

/** Multiply a Polynomial by a Term.
 *
 * Much of the F4 algorithm works with "unevaluated products" (t, f) but eventually
 * they need to be evaluated.  This function does the evaluation.
 */

bool Poly_times_term(Polynomial_t *P, Term_t *t, Polynomial_t *answer, int prime, int rank) {
  int N = 0;
  coeff_t p_coeff;
  if (! Poly_alloc(answer, P->num_terms, rank)) {
    return false;
  }
  for (N=0; N < P->num_terms; N++) {
    Term_multiply(t, P->terms + N, answer->terms + N);
    p_coeff = P->coefficients[N];
    p_coeff.total_degree = Term_total_degree(answer->terms + N, rank);
    answer->coefficients[N] = p_coeff;
  }
  answer->num_terms = N;
  answer->rank = rank;
  return true;
}

/** Copy a Polynomial, while dividing by the head coefficient.
 */

bool Poly_make_monic(Polynomial_t *P, Polynomial_t *answer, int prime, int rank) {
  int N = 0, a_inverse;
  coeff_t p_coeff;
  if (! Poly_alloc(answer, P->num_terms, rank)) {
    return false;
  }
  a_inverse = inverse_mod(prime, P->coefficients[0].value);
  for (N=0; N < P->num_terms; N++) {
    answer->terms[N] = P->terms[N];
    p_coeff = P->coefficients[N];
    p_coeff.value = multiply_mod(prime, a_inverse, p_coeff.value);
    answer->coefficients[N] = p_coeff;
  }
  answer->num_terms = N;
  answer->rank = rank;
  return true;
}

/** A global zero polynomial.
 */

Polynomial_t zero = {.num_terms = 0, .terms = NULL, .coefficients = NULL};

/** The basic row operation
 *
 * Assume that f and g are both non-zero and have the same leading term.
 * Kill a term of g by subtracting an integer multiple of f.
 *
 * First subtract a*f from g where a is the coefficient in g of the head term of
 * f.  We assume that f is monic here, so after this row op, the new value of g
 * will not have a term equal to the head term of f.  In the case where f and g
 * have the same head term, this operation may not produce a monic result, so
 * we then divide g - a*f by its leading coefficient.
 */
static inline bool row_op(Polynomial_t *f, Polynomial_t *g, int g_coeff, int prime, int rank) {
  Polynomial_t new_row = zero;
  /* The coefficient should have been normalized to lie in [0, p).*/
  int a = prime - g_coeff;
  if (! Poly_p_plus_aq(g, a, f, &new_row, prime, rank)) {
    return false;
  }
  if (new_row.coefficients != NULL && new_row.coefficients[0].value != 1) {
    int a_inverse = inverse_mod(prime, new_row.coefficients[0].value);
    for (int n = 0; n < new_row.num_terms; n++) {
      int coeff = new_row.coefficients[n].value;
      new_row.coefficients[n].value = multiply_mod(prime, a_inverse, coeff);
    }
  }
  Poly_free(g);
  *g = new_row;
  return true;
}

/** Echelon Form
 * 
 * Input an array P of Polynomial pointers, and an array answer of unitialized
 * polynomials.  Both arrays should have size num_rows.  Each Polynomial gets
 * divided by its leading coefficient while being copied into the answer
 * array.  Then row operations are performed on the answer array until each
 * leading term occurs in exactly one row.
 */

bool Poly_echelon(Polynomial_t **P, Polynomial_t *answer, int prime, int rank,
                  size_t num_rows) {
  int i, j, coeff;
  Polynomial_t *row_i, *row_j;
  Term_t* head;
  for(i = 0; i < num_rows; i++) {
    *(answer + i) = zero;
    if (! Poly_make_monic(*(P+i), answer+i, prime, rank)) {
      return false;
    }
  }
  for (i = 0; i < num_rows; i++) {
    row_i = answer + i;
    if (row_i->num_terms == 0) continue;
    head = row_i->terms;
    for (j = 0; j < num_rows; j++) {
      if (i == j) continue;
      row_j = answer + j;
      if (row_j->num_terms == 0) continue;
      coeff = Poly_coeff(row_j, head, rank);
      if (coeff != 0) {
        if (! row_op(row_i, answer+j, coeff, prime, rank)) {
          return false;
        }
      }
    }
  }
  return true;
}
