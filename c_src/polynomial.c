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

/** Basic Polynomial operations
 * 
 * This file contains support for the Cython extension class Polynomial.  All
 * functions in this file expect a Polynomial_t of standard flavor.  Functions
 * to support the echelon form computation are in the file echelon.c.
 */

#include "snaffour.h"
#include <assert.h>

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

static inline void x_gcd(int x, int y, gcd_t* answer) {
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
 * We assume that p is odd.
 */


/** Compute the standard representative of the inverse of x modulo p.
 */

int inverse_mod(int p, int x) {
  int a;
  gcd_t xgcd = {.d=0, .a=0, .b=0};
  x_gcd(x, p, &xgcd);
  a = xgcd.d < 0 ? (-xgcd.a) % p : xgcd.a % p;
  return a < 0 ? p + a : a;
}

/** Multiply two elements of Z/pZ.
 *
 * Arguments should be ints in the interval [0,p).
 */

static inline int multiply_mod(int prime, int x, int y) {
  int64_t prime64 = (int64_t)prime, x64 = (int64_t)x, y64 = (int64_t)y, answer64;
  if (x == 1) {
    return (int)y64;
  } else if (x64 == prime64 - 1) {
    return (int)(prime64 - y64);
  }
  answer64 = x64*y64;
  answer64 = answer64 % prime64;
  return (int)answer64;
}

/** Compute x + ay mod p
 *
 * This is the scalar version of a row operation.
 */

static inline int x_plus_ay_mod(int prime, int x, int a, int y) {
  int64_t prime64 = prime, x64 = x, y64 = y, a64 = a, answer64;
  if (a == 1) {
    answer64 = x64 + y64;
  } else if (a == prime - 1) {
    answer64 = x64 - y64;
  } else {
    answer64 = a64*y64;
    answer64 = answer64 % prime64;
    answer64 += x64;
  }
  if (answer64 < 0) {
    answer64 += prime64;
  } else if (answer64 >= prime64) {
    answer64 -= prime64;
  }
  return (int)answer64;
}

/** \page modP Polynomials
 */

/** Allocate or reallocate memory for a Polynomial with specified maximum size.
 *
 * Sets num_terms to 0 and ensures that there is enough memory for up to "size"
 * terms and coefficients.  These arrays are enlarged with realloc if necessary.
 *
 * Note that allocation and deallocation of Polynomial_s structs is not
 * handled by this library.  That should be done by Cython.
 */

bool Poly_alloc(Polynomial_t* P, int size, int rank) {
  int old_size = P->max_size;
  if (size > old_size) {
    if (NULL == (P->terms = realloc(P->terms, sizeof(Term_t)*size))) {
      goto oom;
    }
    if (NULL == (P->coefficients = realloc(P->coefficients, sizeof(coeff_t)*size))) {
      goto oom;
    }
    P->max_size = size;
  }
  P->num_terms = 0;
  P->rank = rank;
  return true;

 oom:
  free(P->terms);
  free(P->coefficients);
  return false;
}

/** Free the terms and coefficients of a Polynomial, if not NULL.
 *
 * The pointers are set to NULL after freeing the arrays and num_terms is
 * set to 0.  The flavor does not change!
 */

void Poly_free(Polynomial_t* P) {
  free(P->terms);
  free(P->coefficients);
  P->num_terms = 0;
  P->rank = 0;
  P->terms = NULL;
  P->coefficients = NULL;
}

/** Copy a Polynomial's data into another Polynomial.
 *
 */
void Poly_copy(Polynomial_t* src, Polynomial_t* dest) {
  int i;
  /* First make sure there is enough room in the destination. */
  Poly_alloc(dest, src->num_terms, src->rank);
  for (i=0; i < src->num_terms; i++) {
    dest->terms[i] = src->terms[i];
  }
  for (i=0; i < src->num_terms; i++) {
    dest->coefficients[i] = src->coefficients[i];
  }
  dest->num_terms = src->num_terms;
}

/** Print a Polynomial in a crude way, without variable names.
 */

void Poly_print(Polynomial_t* P, int rank) {
  int i;
  if (P->num_terms == 0) {
    printf("0\n");
  } else {
    for (i=0; i<P->num_terms; i++) {
      printf("%d*", P->coefficients[i].value);
      Term_print(&P->terms[i], rank);
      printf("\n");
    }
  }
}

/** Static function to compare two Terms in two Polynomials
 *
 * Compare the pth term of P to the qth term of Q and return an integer which is
 * < 0 if the first one is smaller, > 0 if the second one is smaller and 0 if
 * they are equal.  The 0 term is considered larger than any non-zero term.
 */

static inline int Poly_compare_terms(Polynomial_t *P, int p, Polynomial_t *Q, int q) {
  if (P->num_terms > 0 && Q->num_terms > 0) {
    int P_index = P->coefficients[p].column_index;
    int Q_index = Q->coefficients[q].column_index;
    if (P_index >= 0 && Q_index >= 0) {
      /* If both column indexes are non-negative, use them for the comparison. */
      return P_index - Q_index;
    } else {
      /*
       * Otherwise do the comparison from scratch.  This only happens for the
       * standard flavor
       */
      Term_t *P_term = P->terms + p, *Q_term = Q->terms + q;
      int td1, td2, td_cmp;
      td1 = Term_total_degree(P_term, P->rank);
      td2 = Term_total_degree(Q_term, Q->rank);
      td_cmp = td1 - td2;
      return td_cmp == 0 ? Term_revlex_diff(Q_term, P_term, P->rank) : td_cmp;
    }
  } else if (P->num_terms == 0 && Q->num_terms == 0) {
    return 0;
  } else if (P->num_terms == 0) {
    return 1;
  } else {
      return -1;
    }
}

/** Static function to compute P + a*Q for an element a of the coefficient field.
 *
 * This is the core of the row operation used to reduce a matrix to echelon form
 * and also handles addition and subtraction (by taking a=1 or a=-1).
 *
 * The operands must have the standard flavor.
 *
 * NOTE: P and Q must point to different polynomials for this to work.
 */

static inline bool Poly_p_plus_aq(Polynomial_t* P, int a, Polynomial_t* Q,
                                  Polynomial_t* answer, int prime, int rank) {
  int size = P->num_terms + Q->num_terms, p = 0, q = 0, N = 0, cmp;
  coeff_t p_coeff, q_coeff;
  int new_value;
  if (! Poly_alloc(answer, size, rank)) {
    return false;
  }
  p_coeff = P->coefficients[0];
  q_coeff = Q->coefficients[0];
  while (p < P->num_terms && q < Q->num_terms) {
    cmp = Poly_compare_terms(P, p, Q, q);
    if (cmp > 0) { /* deg P > deg Q */
      answer->terms[N] = P->terms[p];
      answer->coefficients[N++] = p_coeff;
      p_coeff = P->coefficients[++p];
    } else if (cmp < 0) { /* deg P < deg Q */
      answer->terms[N] = Q->terms[q];
      new_value = multiply_mod(prime, a, q_coeff.value);
      answer->coefficients[N].column_index = q_coeff.column_index;
      answer->coefficients[N++].value = new_value;
      q_coeff = Q->coefficients[++q];
    } else { /* deg P == deg Q */
      new_value = x_plus_ay_mod(prime, p_coeff.value, a, q_coeff.value);
      if (new_value != 0) {
	answer->terms[N] = P->terms[p];
	answer->coefficients[N].column_index = p_coeff.column_index;
	answer->coefficients[N++].value = new_value;
      }
      p_coeff = P->coefficients[++p];
      q_coeff = Q->coefficients[++q];
    }
  }
  /* At most one of these two loops will be non-trivial. */
  for (; q < Q->num_terms; q++, N++) {
    answer->terms[N] = Q->terms[q];
    q_coeff = Q->coefficients[q];
    new_value = multiply_mod(prime, a, q_coeff.value);
    answer->coefficients[N].column_index = q_coeff.column_index;
    answer->coefficients[N].value = new_value;
  }
  for (; p < P->num_terms; p++, N++) {
    answer->terms[N] = P->terms[p];
    answer->coefficients[N] = P->coefficients[p];
  }
  answer->num_terms = N;
  answer->rank = rank;
  return true;
}


/** Add Polynomials P and Q and store the result in answer.
 *
 * The operands must have normal flavor.
 *
 * The work is done by Poly_p_plus_aq, but we need to deal with the special case
 * where P and Q are the same Polynomial.
 */

bool Poly_add(Polynomial_t* P, Polynomial_t* Q, Polynomial_t* answer,
	      int prime, int rank) {
  if (P == Q) {
    int N;
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
 * The operands must have normal flavor.
 *
 * The work is done by Poly_p_plus_aq, but we need to deal with the special case
 * where P and Q are the same Polynomial (by returning a 0 Polynomial).
 */

bool Poly_sub(Polynomial_t* P, Polynomial_t* Q, Polynomial_t* answer,
	      int prime, int rank) {
  if (P == Q) {
    answer->num_terms = 0;
    return true;
  }
  return Poly_p_plus_aq(P, prime - 1, Q, answer, prime, rank);
}

/** Determine if two Polynomials are equal.
 *
 */

bool Poly_equals(Polynomial_t* P, Polynomial_t *Q) {
  if (P->num_terms != Q->num_terms) {
    return false;
  }
  for (int N = 0; N < P->num_terms; N++) {
    if (! Term_equals(P->terms + N, Q->terms + N)) {
      return false;
    }
    if (P->coefficients[N].value != Q->coefficients[N].value) {
      return false;
    }
  }
  return true;
}

/** Use bisection to find the index of the given term in P->terms for a Polynomial P.
 *
 * Return false if the term is not in P->terms.
 */

static inline bool find_index(Polynomial_t* P, Term_t* t, int t_td, int rank,
		      int bottom, int top, int* index) {
  int middle;
  int td;
  Term_t *terms = P->terms;
  if (top - bottom == 1) {
    if (Term_equals(t, terms + bottom)) {
	*index = bottom;
	return true;
      } else {
	return false;
      }
  }
  middle = (top + bottom) >> 1;
  td = Term_total_degree(terms + middle, rank);
  if (td < t_td || (td == t_td && Term_revlex_diff(terms + middle, t, rank) > 0)) {
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
 * The Polynomial must have the standard flavor.
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
    answer->coefficients[N] = p_coeff;
  }
  answer->num_terms = N;
  answer->rank = rank;
  return true;
}

/** Multiply a Polynomial by an int.
 *
 *  The Polynomial must have the standard flavor.
 */

bool Poly_times_int(Polynomial_t *P, int a, Polynomial_t *answer, int prime, int rank) {
  int N = 0;
  coeff_t p_coeff;
  if (! Poly_alloc(answer, P->num_terms, rank)) {
    return false;
  }
  for (N=0; N < P->num_terms; N++) {
    answer->terms[N] = P->terms[N];
    p_coeff = P->coefficients[N];
    p_coeff.value = multiply_mod(prime, a, p_coeff.value);
    answer->coefficients[N] = p_coeff;
  }
  answer->num_terms = N;
  answer->rank = rank;
  return true;
}

/** Divide a polynomial by its head coefficient.
 *
 * The coefficients of P are modified in place.
 */

void Poly_make_monic(Polynomial_t *P, int prime, int rank) {
  int N = 0, a_inverse = inverse_mod(prime, P->coefficients[0].value);
  for (N=0; N < P->num_terms; N++) {
    P->coefficients[N].value = multiply_mod(prime, a_inverse,
                                            P->coefficients[N].value);
  }
}

/*
 * Sort an array of type Polynomial_t, in place, by head term.  Zero polynomials
 * go to the end;
 *
 */

static inline int compare_heads(const void* p1, const void* p2) {
  Polynomial_t* P1 = (Polynomial_t*)p1;
  Polynomial_t* P2 = (Polynomial_t*)p2;
  return Poly_compare_terms(P1, 0, P2, 0);
}

static inline int compare_heads_dec(const void* p1, const void* p2) {
  Polynomial_t* P1 = (Polynomial_t*)p1;
  Polynomial_t* P2 = (Polynomial_t*)p2;
  return Poly_compare_terms(P2, 0, P1, 0);
}

void Poly_sort(Polynomial_t *P, int num_polys, bool increasing) {
  if (increasing) {
    qsort(P, num_polys, sizeof(Polynomial_t), compare_heads);
  } else {
    qsort(P, num_polys, sizeof(Polynomial_t), compare_heads_dec);
  }
}

/*
 * Allocate an array and fill it with all of the distinct terms, in descending
 * order, that appear in the input array of Polynomials.
 */

bool Poly_terms(Polynomial_t *P, int num_polys, Term_t **answer, int* answer_size,
		int rank) {
  if (num_polys == 1) {
    *answer = (Term_t*)malloc(P[0].num_terms*sizeof(Term_t));
    if (*answer == NULL) {
      return false;
    }
    memcpy((void*)*answer, P[0].terms, P[0].num_terms*sizeof(Term_t));
    *answer_size = P[0].num_terms;
    return true;
  }
  if (num_polys == 2) {
      return Term_merge(P[0].terms, P[1].terms, P[0].num_terms, P[1].num_terms,
			answer, answer_size, rank);
  }
  int left = num_polys / 2, right = num_polys - left;
  Term_t *left_answer, *right_answer;
  int left_size, right_size;
  bool result;
  if (! Poly_terms(P, left, &left_answer, &left_size, rank)) {
    return false;
  }
  if (! Poly_terms(P + left, right, &right_answer, &right_size, rank)) {
    free(left_answer);
    return false;
  }
  result = Term_merge(left_answer, right_answer, left_size, right_size,
  		      answer, answer_size, rank);
  free(left_answer);
  free(right_answer);
  return result;
}
