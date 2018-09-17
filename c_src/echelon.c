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

/** Echelon forms
 * 
 * This file contains functions used by the echelon form computation.
 */

#include "snaffour.h"

/* Getter and setter macros for typedef coeff_t row_coeff_t */
/*
#define SET_COEFF(x, new_coeff) do {x.value = new_coeff;} while (0)
#define GET_COEFF(x) (x.value)
#define SET_COLUMN(x, new_column) do {x.column_index = new_column;} while (0)
#define GET_COLUMN(x) (x.column_index)
#define ZERO_COEFF ({ .value : 0, .column_index : 0 }) 
*/

/* Getter and setter macros for typedef int64_t row_coeff_t */
#define SET_COEFF(x, new_coeff) do {                            \
    x = ((int64_t)new_coeff << 32) | (x & 0xffffffffL);         \
  } while (0)
#define GET_COEFF(x) ((int)(x >> 32))
#define SET_COLUMN(x, new_column) do {                  \
    x = (x & ~0xffffffffL) | ((int64_t)new_column);     \
  } while (0)
#define GET_COLUMN(x) ((int)(x & 0xffffffffL))
#define ZERO_COEFF (0)

/** Global zero values for initialization.
 */

Polynomial_t zero_poly = {.num_terms = 0,
                          .max_size = 0,
                          .terms = NULL,
                          .coefficients = NULL};

Row_t zero_row = {.num_terms = 0,
                  .max_size = 0,
                  .coefficients = NULL,
                  .term_table = NULL};

/** Memory allocation.
 *
 */

static bool Row_alloc(Row_t* P, int size) {
  int old_size = P->max_size;
  if (size > old_size) {
    if (NULL == (P->coefficients = realloc(P->coefficients,
                                           sizeof(row_coeff_t)*size))) {
      goto oom;
    }
    P->max_size = size;
  }
  P->num_terms = 0;
  return true;

 oom:
  free(P->term_table);
  free(P->coefficients);
  return false;
}

static void Row_free(Row_t* P) {
  free(P->coefficients);
  P->num_terms = 0;
  P->coefficients = NULL;
}

/** Montgomery arithmetic.
 *
 * Convention: for low-level arithmetic operations we pass the
 * required auxiliary constants as 32 bit int values, rather
 * than passing the entire MConstants struct.  The intent is
 * to make life easier for the optimizer when inlining these
 * operations, although it is unclear whether it does that.
 */


static inline MConstants_t montgomery_init(int prime) {
  int64_t prime64 = prime, radix64 = M_RADIX64;
  int64_t R_squared64 = (radix64 * radix64) % prime64;
  MConstants_t C;
  C.prime = prime;
  /* mu = -1/p mod R */
  C.mu = inverse_mod(M_RADIX, (int)(radix64 - prime64));
  /* R mod p, i.e. M(1) (used to check if a row is monic). */
  C.R_mod_p = (int)(radix64 % prime64);
  /* R^2 mod p (used to compute Montgomery representatives). */
  C.R_squared = (int)R_squared64;
  /* R^3 mod p (used to compute Montgomery inverses). */
  int64_t R_cubed64 = (R_squared64 * radix64) % prime64;
  C.R_cubed = (int)R_cubed64;
  return C;
}

static inline int montgomery_multiply(int x, int y, int prime, int mu) {
  int64_t x64 = x, y64 = y, prime64 = prime, mu64 = mu, answer64;
  answer64 = M_REDUCE(x64*y64, prime64, mu64);
  if (answer64 >= prime64) {
    answer64 -= prime64;
  }
  return (int)answer64;
}

static inline int montgomery_inverse(int x, int prime, int mu, int R_cubed) {
  int x_inverse = inverse_mod(prime, x);
  int64_t x_inverse64 = x_inverse, R_cubed64 = R_cubed, prime64 = prime, mu64 = mu;
  /* Montgomery multiply by R^3 to get the Montgomery inverse. */
  int64_t answer64 = M_REDUCE(x_inverse64*R_cubed64, prime64, mu64);
  if (answer64 >= prime64) {
    answer64 -= prime64;
  }
  return (int)answer64;
}

static inline int montgomery_x_plus_ay(int x, int a, int y, int prime, int mu) {
  int64_t x64 = x, y64 = y, a64 = a, prime64 = prime, mu64 = mu, answer64;
  answer64 = M_REDUCE(a64*y64, prime64, mu64);
  if (answer64 >= prime64) {
    answer64 -= prime64;
  }
  answer64 += x64;
  if (answer64 >= prime64) {
    answer64 -= prime64;
  }
  return (int)answer64;
}

/** Divide the coefficients of a row by the leading coefficient.
 *
 * The coefficients of P are modified in place.
 */

static void Row_make_monic(Row_t *P, MConstants_t C) {
  register row_coeff_t coeff = P->coefficients[0];
  int N = 0;
  int factor = montgomery_inverse(GET_COEFF(coeff), C.prime, C.mu, C.R_cubed);
  for (N = 0; N < P->num_terms; N++) {
    coeff = P->coefficients[N];
    SET_COEFF(coeff, montgomery_multiply(factor, GET_COEFF(coeff), C.prime, C.mu));
    P->coefficients[N] = coeff;
  }
}

/** Convert a Polynomial which is in transitional state to a Row.
 *
 * This assumes that the source's column indexes have been set to indexes into a
 * table containing all of its terms.  The destination should have been
 * initialized to zero.
 */
static inline bool Poly_to_Row(Polynomial_t* src, Row_t* dest, Term_t* table,
                                 MConstants_t C) {
  register row_coeff_t coeff = ZERO_COEFF;
  int i;
  dest->term_table = table;
  /* First make sure there is enough room in the destination. */
  if (!Row_alloc(dest, src->num_terms)) {
    return false;
  }
  for (i = 0; i < src->num_terms; i++) {
    /* Montgomery multiplication by R^2 converts standard to Montgomery. */
    SET_COEFF(coeff, montgomery_multiply(src->coefficients[i].value,
                                         C.R_squared, C.prime, C.mu));
    SET_COLUMN(coeff, src->coefficients[i].column_index);
    dest->coefficients[i] = coeff;
  }
  dest->num_terms = src->num_terms;
  return true;
}

/** Convert a Row into a Polynomial
 *
 * The Polynomial is expanded if necessary and the Row is freed.
 */

static inline bool Row_to_Poly(Row_t* src, Polynomial_t* dest, int rank,
                               MConstants_t C) {
  register row_coeff_t coeff;
  int i;
  
  *dest = zero_poly;
  if (!Poly_alloc(dest, src->num_terms, rank)) {
    return false;
  }
  dest->num_terms = src->num_terms;
  for (i = 0; i < dest->num_terms; i++) {
    coeff = src->coefficients[i];
    dest->terms[i] = src->term_table[GET_COLUMN(coeff)];
    /* Montgomery multiplication by 1 converts Montgomery to standard. */
    int value = montgomery_multiply(GET_COEFF(coeff), 1, C.prime, C.mu);
    dest->coefficients[i].value = value;
    dest->coefficients[i].column_index = NO_INDEX;
  }
  return true;
}

/** The basic row operation
 *
 * Assume that Q is non-zero, and that the head term of Q has nonzero
 * coefficient in P.  Kill that term of P by subtracting a scalar multiple of
 * Q.  If the head coefficient of Q is a and the corresponding coefficient of
 * P is b, then the computed answer is P - (b/a)*Q.
 *
 * This is the most performance critical step in this implementation.  In the
 * Cyclic-8 example, reducing the largest matrix (3360 x 9326) in step 20
 * (degree 13) executes the core loop of this function about 1367373602
 * times in 7.33 seconds at 3.0GHz. That is an average of about 16.1 cycles
 * per iteration of the loop.  There are many possible way to arrange the
 * loop with the aim of coaxing the gcc optimizer into writing fast code.
 * The results are very unstable: minor rearrangements which are obviously
 * logically equivalent can easily lead to a 50% slowdown.  The version
 * below is the fastest that we have stumbled upon.
 */

static inline bool row_op(Row_t *Q, Row_t *P, Row_t *answer, int P_coeff,
                          MConstants_t C) {
  int inv = montgomery_inverse(GET_COEFF(Q->coefficients[0]),
                               C.prime, C.mu, C.R_cubed);
  /* Multiply by b and negate. Note that p - M(X) = M(p - X) = M(-X). */
  int factor = C.prime - montgomery_multiply(inv, P_coeff, C.prime, C.mu);
  int size = P->num_terms + Q->num_terms + 1; /* see the comment below */
  register int cmp;
  register row_coeff_t p_coeff, q_coeff;
  register row_coeff_t *p_ptr = P->coefficients, *q_ptr = Q->coefficients;
  register row_coeff_t *ans_ptr;
  register row_coeff_t *p_done = P->coefficients + P->num_terms;
  register row_coeff_t *q_done = Q->coefficients + Q->num_terms;
  register int64_t factor64 = factor, prime64 = C.prime, mu64 = C.mu, temp64;
  if (! Row_alloc(answer, size)) {
    return false;
  }
  answer->term_table = P->term_table;
  ans_ptr = answer->coefficients;
  /* It seems fastest to prefetch the coefficients, but this requires loading
   * a coefficient with index num_terms.  This is why we added 1 when computing
   * how much memory to allocate above.
   */
  p_coeff = *p_ptr;
  q_coeff = *q_ptr;
  while (p_ptr < p_done && q_ptr < q_done) {
    cmp = GET_COLUMN(p_coeff) - GET_COLUMN(q_coeff);
    if (cmp > 0) { /* deg p_coeff > deg q_coeff */
      *ans_ptr++ = p_coeff;
      p_coeff = *++p_ptr;
    }
    if (cmp == 0) { /* deg p_coeff == deg q_coeff */
      temp64 = M_REDUCE(factor64*GET_COEFF(q_coeff), prime64, mu64);
      if (temp64 >= prime64) {
        temp64 -= prime64;
      }
      temp64 += GET_COEFF(p_coeff);
      if (temp64 >= prime64) {
        temp64 -= prime64;
      }
      if (temp64 != 0) {
        SET_COEFF(p_coeff, temp64);
        *ans_ptr++ = p_coeff;
      }
      p_coeff = *++p_ptr;
      q_coeff = *++q_ptr;
    }
    if (cmp < 0) { /* deg p_coeff < deg q_coeff */
      temp64 = M_REDUCE(factor64*GET_COEFF(q_coeff), prime64, mu64);
      if (temp64 >= prime64) {
        temp64 -= prime64;
      }
      SET_COEFF(q_coeff, temp64);
      *ans_ptr++ = q_coeff;
      q_coeff = *++q_ptr;
    }
  }
  /* At most one of these two loops will be non-trivial. */
  while (p_ptr < p_done) {
    /* We have already fetched the next p_coeff. */
    *ans_ptr++ = *p_ptr++;
  }
  while (q_ptr < q_done) {
    /* We have already fetched the next q_coeff. */
    q_coeff = *q_ptr++;
    temp64 = M_REDUCE(factor64*GET_COEFF(q_coeff), prime64, mu64);
    if (temp64 >= prime64) {
      temp64 -= prime64;
    }
    SET_COEFF(q_coeff, temp64);
    *ans_ptr++ = q_coeff;
  }
  answer->num_terms = ans_ptr - answer->coefficients;
  return true;
}

/** Static function to compare two Terms in two Rows.
 *
 * Compare the pth term of P to the qth term of Q and return an integer which is
 * < 0 if the first one is smaller, > 0 if the second one is smaller and 0 if
 * they are equal.  The 0 term is considered larger than any non-zero term.
 */

static inline int compare_heads(const void* p1, const void* p2) {
  Row_t *P1 = (Row_t*)p1, *P2 = (Row_t*)p2;
  return GET_COLUMN(P1->coefficients[0]) - GET_COLUMN(P2->coefficients[0]);
}
    
static inline int compare_heads_dec(const void* p1, const void* p2) {
  Row_t *P1 = (Row_t*)p1, *P2 = (Row_t*)p2;
  return GET_COLUMN(P2->coefficients[0]) - GET_COLUMN(P1->coefficients[0]);
}

/** Merge two arrays of monomials.
 *
 * Used in building the table of terms.
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
 *
 * Used in building the table of terms.
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
 * column_index field of each coefficient.  This prepares the Polynomials for
 * conversion to Rows, which is then carried out.
 *
 * The number of columns in the matrix having the input polynomials as its rows
 * is stored in the int referenced by the num_columns input, and the table of
 * all terms is stored in the term_table input.
 */

static bool Poly_matrix_init(Polynomial_t** P, int num_rows, int* num_columns,
                             Term_t** term_table, Row_t* matrix,
			     int rank, MConstants_t C) {
  int num_monomials = 0, i, j, index;
  monomial_t *pool;
  monomial_array_t monomial_arrays[num_rows], previous, merged;
  Term_t *table = NULL;

  /* Construct the shared table of Terms. */
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

  /* Construct the matrix. */
  for (i = 0; i < num_rows; i++) {
    matrix[i] = zero_row;
    if (!Poly_to_Row(P[i], matrix + i, table, C)) {
      for (j = 0; j < i; j++) {
        Row_free(matrix + i);
      }
      goto oom;
    }
  }
  /* This implementation requires sorting the rows by increasing head term.*/
  qsort(matrix, num_rows, sizeof(Row_t), compare_heads);
  return true;

 oom:
  free(merged.monomials);
  free(pool);
  free(table);
  return false;
}

/** Use bisection to find the coefficient of P with a given column index.
 *
 * Return false if P has a zero coefficient in the column.
 */

static inline bool coeff_in_column(Row_t* P, int column, int bottom,
                                   int top, int* coefficient) {
  int middle;
  if (top - bottom == 1) {
    if (column == GET_COLUMN(P->coefficients[bottom])) {
      *coefficient = GET_COEFF(P->coefficients[bottom]);
      return true;
    } else {
      return false;
    }
  }
  middle = (top + bottom) >> 1;
  if (GET_COLUMN(P->coefficients[middle]) < column) {
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
  int i, j, target;
  Term_t* table = NULL;
  Row_t *row_i, *row_j, *matrix = NULL;
  Row_t buffer = zero_row, tmp;
  int head, max_head = -1;
  MConstants_t C = montgomery_init(prime);

  /* Create the matrix. */
  if (NULL == (matrix = (Row_t*)malloc(num_rows*sizeof(Row_t)))) {
    goto oom;;
  }
  if (!Poly_matrix_init(P, num_rows, num_columns, &table, matrix, rank, C)) {
    goto oom;
  }
  
  /* Initialize the buffer row. */
  buffer.term_table = table;
  
  /* Reduce the matrix to echelon form. */
  for (i = 0; i < num_rows; i++) {
    row_i = matrix + i;
    if (row_i->num_terms == 0) continue;
    head = GET_COLUMN(row_i->coefficients[0]);
    /* Clear the column above the head term of row i.  Since head terms are
     * initially non-decreasing, there is nothing to do if the head term
     * of row i is larger that the head term of all earlier rows.
     */
    if (head <= max_head) {
      for (j = 0; j < i; j++) {
        row_j = matrix + j;
        if (row_j->num_terms == 0) continue;
        if (coeff_in_column(row_j, head, 0, row_j->num_terms, &target)) {
          if (! row_op(row_i, row_j, &buffer, target, C)) {
            return false;
          }
        /* Swap row j with the buffer row. */
          tmp = matrix[j];
          matrix[j] = buffer;
          buffer = tmp;
        }
      }
    } else {
      max_head = head;
    }
    /* Clear the column below the head term of row i. */
    for (j = i + 1; j < num_rows; j++) {
      row_j = matrix + j;
      if (row_j->num_terms == 0) continue;
      if (coeff_in_column(row_j, head, 0, row_j->num_terms, &target)) {
        if (! row_op(row_i, row_j, &buffer, target, C)) {
          return false;
        }
        /* Swap row j with the buffer row. */
        tmp = matrix[j];
        matrix[j] = buffer;
        buffer = tmp;
      }
    }
  }
  Row_free(&buffer);
  
  /*
   * While we are here in C land, let's sort the result by decreasing head term.
   */
  qsort(matrix, num_rows, sizeof(Row_t), compare_heads_dec);

  /*
   * Convert the rows to monic polynomials and free the rows.
   */
  for (i = 0; i < num_rows; i++) {
    row_i = matrix + i;
    if (row_i->num_terms == 0) {
      answer[i] = zero_poly;
      Row_free(row_i);
    } else {
      Row_make_monic(row_i, C);
      if (!Row_to_Poly(row_i, answer + i, rank, C)) {
        for (j = 0; j <= i; j++) {
          Row_free(matrix + j);
        }
        goto oom;
      }
      Row_free(row_i);
    }
  }
  free(matrix);
  free(table);
  return true;

 oom:
  Row_free(&buffer);
  free(matrix);
  free(table);
  return false;
}
