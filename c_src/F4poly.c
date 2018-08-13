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

/**
 * When computing echelon forms over Fp, we use the Montgomery representation
 * of a conjugacy class mod p.  Given a class X, its Montgomery representative
 * M(X) is an element of [0, p) such that M(X) is congruent to RX mod p, where
 * the Montgomery radix R will always be 2^31 for us.
 */

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

#define M_REDUCE(X, mu, p) ((X + (((X & MOD_R)*mu) & MOD_R)*p) >> 31)

/** Multiply two elements of Z/pZ.
 *
 * Arguments should be ints in the interval [0,p).
 *
 * Supplying a value of 0 for the reduction constant mu signals that hardware
 * division should be used to compute the standard representative of x + ay.
 * Otherwise, Montgomery reduction is used to compute the Montgomery
 * representative.  The compact version assumes mu is non-zero and expects
 * 64 bit values for p and mu.
 *
 * This function can also be used to convert between standard and Montgomery
 * representations.  Multiplying by the standard representative of R^2 mod P
 * converts a standard element to its Montgomery representation, and
 * multiplying by 1 converts back.
 */

static inline int multiply_mod(int prime, int x, int y, int mu) {
  int64_t prime64 = (int64_t)prime, x64 = (int64_t)x, y64 = (int64_t)y, mu64, answer64;
  if (mu != 0) {
    mu64 = mu;
    answer64 = M_REDUCE(x64*y64, mu64, prime64);
    /* Note: My tests indicate that this "if" is faster than:
     * answer64 -= (~((answer64 - prime64) >> 63)) & prime64;
     */
    if (answer64 >= prime64) {
      answer64 -= prime64;
    }
  } else { /* Not using a Montgomery representation. */
    if (x == 1) {
      return (int)y64;
    } else if (x64 == prime64 - 1) {
      return (int)(prime64 - y64);
    }
    answer64 = x64*y64;
    answer64 = answer64 % prime64;
  }
  return (int)answer64;
}

static inline int compact_multiply_mod(int64_t prime64, int x, int y,
				       int64_t mu64) {
  int64_t x64 = (int64_t)x, y64 = (int64_t)y, answer64;
    answer64 = M_REDUCE(x64*y64, mu64, prime64);
    if (answer64 >= prime64) {
      answer64 -= prime64;
    }
    return (int)answer64;
}

/** Compute x + ay mod p
 *
 * This is the scalar version of a row operation.  Supplying a value of 0 for
 * the reduction constant mu signals that ordinary division should be used to
 * compute the standard representative of x + ay.  Otherwise, Montgomery
 * reduction is used to compute the Montgomery representative.  The compact
 * version assumes mu is non-zero and expects 64 bit values for p and mu.
 */

static inline int x_plus_ay_mod(int prime, int x, int a, int y, int mu) {
  int64_t prime64 = prime, x64 = x, y64 = y, a64 = a, mu64, answer64;
    if (mu != 0) {
      mu64 = mu;
      answer64 = M_REDUCE(a64*y64, mu64, prime64);
      if (answer64 >= prime64) {
      	answer64 -= prime64;
      }
      answer64 += x64;
    } else {  /* Not using a Montgomery representation. */
      if (a == 1) {
	answer64 = x64 + y64;
      } else if (a == prime - 1) {
	answer64 = x64 - y64;
      } else {
	answer64 = a64*y64;
	answer64 = answer64 % prime64;
	answer64 += x64;
      }
    }
  if (answer64 < 0) {
    answer64 += prime64;
  } else if (answer64 >= prime64) {
    answer64 -= prime64;
  }
  return (int)answer64;
}

static inline int compact_x_plus_ay_mod(int64_t prime64, int x, int a, int y,
					int64_t mu64) {
  int64_t x64 = x, y64 = y, a64 = a, answer64;
  answer64 = M_REDUCE(a64*y64, mu64, prime64);
  if (answer64 >= prime64) {
    answer64 -= prime64;
  }
  answer64 += x64;
  if (answer64 >= prime64) {
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
    if (P->table == NULL) {
      if (NULL == (P->terms = realloc(P->terms, sizeof(Term_t)*size))) {
	goto oom;
      }
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
  free(P->table);
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

/** Decompress a compact Polynomial, converting it to normal flavor.
 */

static inline bool Poly_decompress(Polynomial_t* P, int prime, int mu) {
  int i, value;
  int64_t prime64 = prime, mu64 = mu;
  if (P->num_terms > 0) {
    P->terms = (Term_t*)malloc(P->num_terms*sizeof(Term_t));
    if (P->terms == NULL){
      return false;
    }
  }
  for (i = 0; i < P->num_terms; i++) {
    P->terms[i] = P->table[P->coefficients[i].column_index];
    value = compact_multiply_mod(prime64, P->coefficients[i].value, 1, mu64);
    P->coefficients[i].column_index = INDEX_UNSET;
    P->coefficients[i].value = value;
  }
  P->table = NULL;
  return true;
}

/** Copy a Polynomial's data into another Polynomial.
 *
 * The src and dest must have the same flavor and, if compact, must share the
 * same table.
 */
void Poly_copy(Polynomial_t* src, Polynomial_t* dest) {
  int i;
  /* First make sure there is enough room in the destination. */
  Poly_alloc(dest, src->num_terms, src->rank);
  if (src->table == NULL) {
    for (i=0; i < src->num_terms; i++) {
      dest->terms[i] = src->terms[i];
    }
  }
  for (i=0; i < src->num_terms; i++) {
    dest->coefficients[i] = src->coefficients[i];
  }
  dest->num_terms = src->num_terms;
}

/** Compress a Polynomial which is in transitional state.
 *
 * This assumes that the source is almost a Polynomial of normal flavor, but
 * its column indexes have been set to indexes into a table containing all of
 * its terms.  The destination should have been initialized to zero, and will
 * be changed into a Polynomial of compact flavor.
 */
static inline bool Poly_compress(Polynomial_t* src, Polynomial_t* dest,
                          Term_t* table, int prime, int mu, int R_squared) {
  int i, value;
  dest->table = table;
  /* First make sure there is enough room in the destination. */
  if (!Poly_alloc(dest, src->num_terms, src->rank)) {
    return false;
  }
  for (i=0; i < src->num_terms; i++) {
    value = src->coefficients[i].value;
    value = multiply_mod(prime, value, R_squared, mu);
    dest->coefficients[i].column_index = src->coefficients[i].column_index;
    dest->coefficients[i].value = value;
  }
  dest->num_terms = src->num_terms;
  return true;
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
      if (P->table != NULL) {
	Term_print(P->table + P->coefficients[i].column_index, rank);
      } else {
	Term_print(&P->terms[i], rank);
      }
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
 * The operands must have the same flavor and, if compact, must share the same
 * table.
 *
 * NOTE: P and Q must point to different polynomials for this to work.
 */

static inline bool Poly_p_plus_aq_normal(Polynomial_t* P, int a, Polynomial_t* Q,
                                  Polynomial_t* answer, int prime, int rank, int mu) {
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
      new_value = multiply_mod(prime, a, q_coeff.value, mu);
      answer->coefficients[N].column_index = q_coeff.column_index;
      answer->coefficients[N++].value = new_value;
      q_coeff = Q->coefficients[++q];
    } else { /* deg P == deg Q */
      new_value = x_plus_ay_mod(prime, p_coeff.value, a, q_coeff.value, mu);
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
    new_value = multiply_mod(prime, a, q_coeff.value, mu);
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

static inline bool Poly_p_plus_aq_compact(Polynomial_t* P, int a, Polynomial_t* Q,
                                  Polynomial_t* answer, int prime, int rank, int mu) {
  int size = P->num_terms + Q->num_terms, p = 0, q = 0, N = 0, cmp, new_value;
  coeff_t p_coeff, q_coeff;
  int64_t prime64 = prime, mu64 = mu;
  if (! Poly_alloc(answer, size, rank)) {
    return false;
  }
  answer->table = P->table;
  p_coeff = P->coefficients[0];
  q_coeff = Q->coefficients[0];
  while (p < P->num_terms && q < Q->num_terms) {
    cmp = p_coeff.column_index - q_coeff.column_index;
    if (cmp > 0) { /* deg P > deg Q */
      answer->coefficients[N++] = p_coeff;
      p_coeff = P->coefficients[++p];
    } else if (cmp < 0) { /* deg P < deg Q */
      new_value = compact_multiply_mod(prime64, a, q_coeff.value, mu64);
      answer->coefficients[N].column_index = q_coeff.column_index;
      answer->coefficients[N++].value = new_value;
      q_coeff = Q->coefficients[++q];
    } else { /* deg P == deg Q */
      new_value = compact_x_plus_ay_mod(prime64, p_coeff.value, a, q_coeff.value, mu64);
      if (new_value != 0) {
	answer->coefficients[N].column_index = p_coeff.column_index;
	answer->coefficients[N++].value = new_value;
      }
      p_coeff = P->coefficients[++p];
      q_coeff = Q->coefficients[++q];
    }
  }
  /* At most one of these two loops will be non-trivial. */
  for (; q < Q->num_terms; q++, N++) {
    q_coeff = Q->coefficients[q];
    new_value = compact_multiply_mod(prime64, a, q_coeff.value, mu64);
    answer->coefficients[N].column_index = q_coeff.column_index;
    answer->coefficients[N].value = new_value;
  }
  for (; p < P->num_terms; p++, N++) {
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
      p_coeff.value = x_plus_ay_mod(prime, p_coeff.value, 1, p_coeff.value, 0);
      answer->coefficients[N] = p_coeff;
    }
    answer->num_terms = N;
    answer->rank = rank;
    return true;
  }
  return Poly_p_plus_aq_normal(P, 1, Q, answer, prime, rank, 0);
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
  return Poly_p_plus_aq_normal(P, prime - 1, Q, answer, prime, rank, 0);
}

/** Determine if two Polynomials are equal.
 *
 * The operands must have the same flavor and, if compact, must share
 * the same table.
 */

bool Poly_equals(Polynomial_t* P, Polynomial_t *Q) {
  if (P->num_terms != Q->num_terms) {
    return false;
  }
  if (P->table != NULL) {
    for (int N = 0; N < P->num_terms; N++) {
      if (P->coefficients[N].column_index != Q->coefficients[N].column_index ||
	  P->coefficients[N].value != Q->coefficients[N].value) {
	return false;
      }
    }
  } else {
    for (int N = 0; N < P->num_terms; N++) {
      if (! Term_equals(P->terms + N, Q->terms + N)) {
	return false;
      }
      if (P->coefficients[N].value != Q->coefficients[N].value) {
	return false;
      }
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
  Term_t *terms = P->table == NULL ? P->terms : P->table;
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
    p_coeff.value = multiply_mod(prime, a, p_coeff.value, 0);
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

void Poly_make_monic(Polynomial_t *P, int prime, int rank, int mu, int R_cubed) {
  int N = 0, a_inverse = inverse_mod(prime, P->coefficients[0].value);
  if (mu != 0) {
    /* Use the Montgomery inverse. */
    a_inverse = multiply_mod(prime, a_inverse, R_cubed, mu);
  }
  for (N=0; N < P->num_terms; N++) {
    P->coefficients[N].value = multiply_mod(prime, a_inverse,
                                            P->coefficients[N].value, mu);
  }
}

/** A global zero polynomial.
 */

Polynomial_t zero = {.num_terms = 0,
		     .max_size = 0,
		     .terms = NULL,
		     .coefficients = NULL};

/** The basic row operation
 *
 * Assume that f and g are both non-zero, that f is monic, and that the head term of
 * f appears in g.  Kill that term of g by subtracting a scalar multiple of f.
 *
 * First subtract a*f from g where a is the coefficient in g of the head term of
 * f.  Since f is monic, after this subtraction, the result will not have a term
 * equal to the head term of f.  However, in the case where the cancelled term
 * is the head term of g, i.e. f and g have the same head term, this operation
 * may not produce a monic result. So in that case we then divide g - a*f by its
 * leading coefficient.
 */

static inline bool row_op(Polynomial_t *f, Polynomial_t *g, Polynomial_t *answer,
                          int g_coeff, int prime, int rank, int mu, int R_cubed,
                          int R_mod_p) {
  /* The coefficient should have been normalized to lie in [0, p).*/
  /* Note that p - M(X) = M(p - X) = M(-X). */
  int a = prime - g_coeff;
  if (! Poly_p_plus_aq_compact(g, a, f, answer, prime, rank, mu)) {
    return false;
  }
  /* Make sure the new row is monic. This is always called with non-ero mu!*/
  if (answer->coefficients != NULL && answer->coefficients[0].value != R_mod_p) {
    /* Use the Montgomery inverse.
     * This does several divisions for each row op -- think about that!
     */
    int a_inverse = inverse_mod(prime, answer->coefficients[0].value);
    a_inverse = multiply_mod(prime, a_inverse, R_cubed, mu);
    for (int n = 0; n < answer->num_terms; n++) {
      int coeff = answer->coefficients[n].value;
      answer->coefficients[n].value = multiply_mod(prime, a_inverse, coeff, mu);
    }
  }
  return true;
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

/** Merge two arrays of monomials.
 */

static bool monomial_merge_two(monomial_array_t M0, monomial_array_t M1,
			       monomial_array_t *answer, int rank) {
  int size = M0.size + M1.size, p = 0, q = 0, N = 0, td_cmp, revlex_cmp;
  monomial_array_t merged;
  merged.monomials = (monomial_t*)malloc(sizeof(monomial_t)*size);
  if (merged.monomials == NULL) {
    answer->size = 0;
    answer->monomials = NULL;
    return false;
  }
  while (p < M0.size && q < M1.size) {
    td_cmp = M0.monomials[p].total_degree - M1.monomials[q].total_degree;
    if (td_cmp == 0) {
        Term_t *T0 = M0.monomials[p].term, *T1 = M1.monomials [q].term;
        revlex_cmp = Term_revlex_diff(T1, T0, rank);
    } else {
      revlex_cmp = 0;
    }
    if (td_cmp > 0 || revlex_cmp > 0) {
      /* M0[p] > M1[q] */
      merged.monomials[N++] = M0.monomials[p++];
    } else if (td_cmp < 0 || revlex_cmp < 0) {
      /* M1[q] > M0[p]*/
      merged.monomials[N++] = M1.monomials[q++];
    } else {
      /* M0[p] == M1[q] */
      merged.monomials[N++] = M0.monomials[p++];
      merged.monomials[N] = M1.monomials[q++];
      /* Mark the duplicate by setting its column_index to a special value. */
      merged.monomials[N++].coefficient->column_index = DUPLICATE;
    }
  }
  /* At most one of these two loops will be non-trivial. */
  while (p < M0.size) {
    merged.monomials[N++] = M0.monomials[p++];
  }
  while (q < M1.size) {
    merged.monomials[N++] = M1.monomials[q++];
  }
  merged.size = N;
  *answer = merged;
  return true;
}

/** Merge many arrays of monomials by divide and conquer.
 */

static bool monomial_merge(monomial_array_t* M, int num_arrays,
                           monomial_array_t *answer, int rank) {
  if (num_arrays == 1) {
    answer->monomials = (monomial_t*)malloc(M[0].size*sizeof(monomial_t));
    if (answer->monomials == NULL) {
      return false;
    }
    memcpy((void*)answer->monomials, M[0].monomials, M[0].size*sizeof(monomial_t));
    answer->size = M[0].size;
    return true;
  }
  if (num_arrays == 2) {
      return monomial_merge_two(M[0], M[1], answer, rank);
  }
  int left = num_arrays / 2, right = num_arrays - left;
  monomial_array_t left_answer, right_answer;
  left_answer.monomials = NULL;
  right_answer.monomials = NULL;
  bool result;
  if (! monomial_merge(M, left, &left_answer, rank)) {
    return false;
  }
  if (! monomial_merge(M + left, right, &right_answer, rank)) {
    free(left_answer.monomials);
    return false;
  }
  result = monomial_merge_two(left_answer, right_answer, answer, rank);
  free(left_answer.monomials);
  free(right_answer.monomials);
  return result;
}

/* Initialize a matrix
 *
 * Inputs an array P of polynomial pointers.  Sorts all of the terms which
 * appear in any of the polynomials and uses their sort position to set the
 * column_index field of each coefficient.  This prepares the polynomials for
 * conversion to the compact flavor, which is then carried out.
 *
 * The number of columns in the matrix having the input polynomials as its rows
 * is stored in the int referenced by the num_columns input, and the table of
 * all terms is stored in the term_table input.
 */

static bool Poly_matrix_init(Polynomial_t **P, int num_rows, int *num_columns,
                             Term_t **term_table, Polynomial_t *matrix,
			     int prime, int rank, int mu, int R_squared) {
  int num_monomials = 0, i, j, index;
  monomial_t *pool;
  monomial_array_t monomial_arrays[num_rows], previous, merged;
  Term_t *table = NULL;
  for (i = 0; i < num_rows; i++) {
    num_monomials += P[i]->num_terms;
  }
  if (NULL == (pool = (monomial_t*)malloc(num_monomials*sizeof(monomial_t)))) {
    goto oom;
  }
  previous.size = 0; previous.monomials = pool;
  for (i = 0; i < num_rows; i++) {
    monomial_arrays[i].size = P[i]->num_terms;
    monomial_arrays[i].monomials = previous.monomials + previous.size;
    previous = monomial_arrays[i];
  }
  for (i = 0; i < num_rows; i ++) {
    for (j = 0; j < P[i]->num_terms; j++){
      monomial_arrays[i].monomials[j].total_degree = Term_total_degree(P[i]->terms + j, rank);
      monomial_arrays[i].monomials[j].term = P[i]->terms + j;
      monomial_arrays[i].monomials[j].coefficient = P[i]->coefficients + j;
    }
  }
  merged.monomials = NULL;
  if (! monomial_merge(monomial_arrays, num_rows, &merged, rank)) {
    goto oom;
  }
  if (NULL == (table = (Term_t*)malloc(merged.size*sizeof(Term_t)))) {
    goto oom;
  }
  index = 0;
  for (i = merged.size - 1; i >= 0; i--) {
    int saved_index = merged.monomials[i].coefficient->column_index;
    merged.monomials[i].coefficient->column_index = index;
    if (saved_index != DUPLICATE) {
      table[index] = *merged.monomials[i].term;
      index++;
    }
  }
  *num_columns = index;
  if (NULL == (table = realloc(table, index*sizeof(Term_t)))) {
    goto oom;
  }
  *term_table = table;
  free(merged.monomials);
  free(pool);
  for (i = 0; i < num_rows; i++) {
    matrix[i] = zero;
    if (!Poly_compress(P[i], matrix + i, table, prime, mu, R_squared)) {
      for (j = 0; j < i; j++) {
        Poly_free(matrix + i);
      }
      goto oom;
    }
  }
  return true;

 oom:
  free(merged.monomials);
  free(pool);
  free(table);
  return false;
}

/** Use bisection to find the coefficient of P with a given column index.
 *
 * Return false if no non-zero coefficient of P has the column index.
 * When the column indexes are available this is quite a bit faster.
 */

static inline bool coeff_in_column(Polynomial_t* P, int column, int bottom,
                                   int top, int* coefficient) {
  int middle;
  if (top - bottom == 1) {
    if (column == P->coefficients[bottom].column_index) {
	*coefficient = P->coefficients[bottom].value;
	return true;
      } else {
	return false;
      }
  }
  middle = (top + bottom) >> 1;
  if (P->coefficients[middle].column_index < column) {
    return coeff_in_column(P, column, bottom, middle, coefficient);
  } else {
    return coeff_in_column(P, column, middle, top, coefficient);
  }
}

/** Echelon Form
 *
 * Input an array P of Polynomial pointers, and an array answer of unitialized
 * polynomials.  Both arrays should have size num_rows.  Each Polynomial gets
 * divided by its leading coefficient while being copied into the answer
 * array.  Then row operations are performed on the answer array until each
 * leading term occurs in exactly one row.
 */

bool Poly_echelon(Polynomial_t** P, Polynomial_t* answer, int num_rows,
                  int* num_columns, int prime, int rank) {
  int i, j, coeff, mu, R_squared, R_cubed, R_mod_p;
  int64_t prime64 = prime, radix64 = M_RADIX64, R_squared64, R_cubed64;
  Term_t* term_table = NULL;
  Polynomial_t *row_i, *row_j;
  Polynomial_t buffer = zero, tmp;
  int head;
  /* Compute constants needed for the Montgomery representation. */
  /* First, mu, the negative inverse of p mod R.
   * This seems to work, even though M_RADIX is negative as a signed int.
   */
  mu = inverse_mod(M_RADIX, (int)(radix64 - prime64));
  assert(0 <= mu && mu & ~MOD_R == 0); 
  assert((1 + ((int64_t)mu * radix64) % radix64) == 0);
  /* R mod p, i.e. M(1) (used to check if a row is monic). */
  R_mod_p = (int)(radix64 % prime64);
  /* R^2 mod p (used to compute Montgomery representatives). */
  R_squared64 = (radix64 * radix64) % prime64;
  R_squared = (int)R_squared64;
  assert(0 <= R_squared && R_squared < prime);
  assert((R_squared - M_RADIX * M_RADIX) % prime == 0);
  /* R^3 mod p (used to compute Montgomery inverses). */
  R_cubed64 = (R_squared64 * radix64) % prime64;
  R_cubed = (int)R_cubed64;
  assert(0 <= R_cubed && R_cubed < prime);
  assert((R_cubed - M_RADIX * M_RADIX * M_RADIX) % prime == 0);
  
  if (!Poly_matrix_init(P, num_rows, num_columns, &term_table,
                        answer, prime, rank, mu, R_squared)) {
    goto oom;
  }
  buffer.table = term_table;
  if (!Poly_alloc(&buffer, *num_columns, rank)) {
    goto oom;
  }
  for(i = 0; i < num_rows; i++) {
    Poly_make_monic(P[i], prime, rank, mu, R_cubed);
  }
  /*
   * The best pivoting strategy I have found is to sort by increasing head term
   * once, at the beginning.  Doing subsequent sorts seems to cost more than it
   * saves.
   */
  qsort(answer, num_rows, sizeof(Polynomial_t), compare_heads);
  for (i = 0; i < num_rows; i++) {
    row_i = answer + i;
    if (row_i->num_terms == 0) continue;
    head = row_i->coefficients->column_index;
    for (j = 0; j < num_rows; j++) {
      if (i == j) continue;
      row_j = answer + j;
      if (row_j->num_terms == 0) continue;
      if (coeff_in_column(row_j, head, 0, row_j->num_terms, &coeff)) {
        if (! row_op(row_i, row_j, &buffer, coeff, prime, rank, mu, R_cubed, R_mod_p)) {
          return false;
        }
        tmp = answer[j];
        answer[j] = buffer;
        buffer = tmp;
      }
    }
  }
  Poly_free(&buffer);
  /*
   * While we are here in C land, let's sort the result by decreasing head term.
   */
  qsort(answer, num_rows, sizeof(Polynomial_t), compare_heads_dec);
  /*
   * Free unneeded memory and convert the row polynomials back to standard flavor.
   */
  for (i = 0; i < num_rows; i++) {
    row_i = answer + i;
    if (row_i->num_terms == 0) {
      Poly_free(row_i);
      continue;
    }
    if (!Poly_decompress(row_i, prime, mu)) {
      for (j = 0; j <= i; j++) {
	Poly_free(answer + j);
      }
      goto oom;
    }
  }
  free(term_table);
  return true;

 oom:
  Poly_free(&buffer);
  free(term_table);
  return false;
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
