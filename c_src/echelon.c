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

/** Global zero values.
 */

Polynomial_t zero_poly = {.num_terms = 0,
                          .max_size = 0,
                          .terms = NULL,
                          .coefficients = NULL};

Row_t zero_row = {.num_terms = 0,
                  .max_size = 0,
                  .coefficients = NULL,
                  .term_table = NULL};

static bool Row_alloc(Row_t* P, int size) {
  int old_size = P->max_size;
  if (size > old_size) {
    if (NULL == (P->coefficients = realloc(P->coefficients, sizeof(coeff_t)*size))) {
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

/** Arithmetic with Montgomery representations.
 *
 * Convention: for low-level arithmetic operations we pass the
 * required auxiliary constants as 32 bit int values, rather
 * than passing the entire MConstants struct.  The intent is
 * to make life easier for the optimizer when inlining these
 * operations, although it is unclear whether it does that.
 */

static inline int montgomery_multiply(int x, int y, int prime, int mu) {
  int64_t x64 = x, y64 = y, prime64 = prime, mu64 = mu, answer64;
  answer64 = M_REDUCE(x64*y64, mu64, prime64);
  if (answer64 >= prime64) {
    answer64 -= prime64;
  }
  return (int)answer64;
}

static inline int montgomery_inverse(int x, int prime, int mu, int R_cubed) {
  int x_inverse = inverse_mod(prime, x);
  int64_t R_cubed64 = R_cubed, prime64 = prime, mu64 = mu, answer64;
  int64_t x_inverse64 = montgomery_multiply(x_inverse, R_cubed, prime, mu);
  answer64 = M_REDUCE(x_inverse64*R_cubed64, mu64, prime64);
  if (answer64 >= prime64) {
    answer64 -= prime64;
  }
  return (int)answer64;
}

static inline int montgomery_x_plus_ay(int x, int a, int y, int prime, int mu) {
  int64_t x64 = x, y64 = y, a64 = a, prime64 = prime, mu64 = mu, answer64;
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

/** Divide the coefficients of a row by the leading coefficient.
 *
 * The coefficients of P are modified in place.
 */

static void Row_make_monic(Row_t *P, MConstants_t C) {
  int N = 0, a = P->coefficients->value;
  int a_inverse = montgomery_inverse(a, C.prime, C.mu, C.R_cubed);
  for (N=0; N < P->num_terms; N++) {
    P->coefficients[N].value = montgomery_multiply(
      a_inverse, P->coefficients[N].value, C.prime, C.mu);
  }
}

/** Convert a Polynomial which is in transitional state to a Row.
 *
 * This assumes that the source's column indexes have been set to indexes into a
 * table containing all of its terms.  The destination should have been
 * initialized to zero.
 */
static inline bool Poly_compress(Polynomial_t* src, Row_t* dest,
                                 Term_t* table, char* pivot_info, int num_pivots,
                                 MConstants_t C) {
  int i, value;
  dest->term_table = table;
  /* First make sure there is enough room in the destination. */
  if (!Row_alloc(dest, src->num_terms)) {
    return false;
  }
  for (i=0; i < src->num_terms; i++) {
    value = src->coefficients[i].value;
    value = montgomery_multiply(value, C.R_squared, C.prime, C.mu);
    dest->coefficients[i].column_index = src->coefficients[i].column_index;
    dest->coefficients[i].value = value;
  }
  dest->num_terms = src->num_terms;
  return true;
}

/** Convert a Row into a Polynomial
 *
 * Allocates a new Polynomial and frees the Row.
 */

static inline bool Poly_decompress(Row_t* P, Polynomial_t* Q, int rank, MConstants_t C) {
  int i, value;
  *Q = zero_poly;
  if (!Poly_alloc(Q, P->num_terms, rank)) {
    return false;
  }
  Q->num_terms = P->num_terms;
  /* Copy the terms and convert the coefficients from the Montgomery
   * representation to the standard representation.
   */
  for (i = 0; i < Q->num_terms; i++) {
    Q->terms[i] = P->term_table[P->coefficients[i].column_index];
    value = montgomery_multiply(P->coefficients[i].value, 1, C.prime, C.mu);
    Q->coefficients[i].column_index = INDEX_UNSET;
    Q->coefficients[i].value = value;
  }
  Row_free(P);
  return true;
}

static inline bool Row_p_plus_aq(Row_t* P, int a, Row_t* Q, Row_t* answer,
                                 MConstants_t C) {
  int size = P->num_terms + Q->num_terms, p = 0, q = 0, N = 0, cmp, new_value;
  coeff_t p_coeff, q_coeff;
  if (! Row_alloc(answer, size)) {
    return false;
  }
  answer->term_table = P->term_table;
  p_coeff = P->coefficients[0];
  q_coeff = Q->coefficients[0];
  while (p < P->num_terms && q < Q->num_terms) {
    cmp = p_coeff.column_index - q_coeff.column_index;
    if (cmp > 0) { /* deg P > deg Q */
      answer->coefficients[N++] = p_coeff;
      p_coeff = P->coefficients[++p];
    } else if (cmp < 0) { /* deg P < deg Q */
      new_value = montgomery_multiply(a, q_coeff.value, C.prime, C.mu);
      answer->coefficients[N].column_index = q_coeff.column_index;
      answer->coefficients[N++].value = new_value;
      q_coeff = Q->coefficients[++q];
    } else { /* deg P == deg Q */
      new_value = montgomery_x_plus_ay(p_coeff.value, a, q_coeff.value,
                                       C.prime, C.mu);
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
    new_value = montgomery_multiply(a, q_coeff.value, C.prime, C.mu);
    answer->coefficients[N].column_index = q_coeff.column_index;
    answer->coefficients[N].value = new_value;
  }
  for (; p < P->num_terms; p++, N++) {
    answer->coefficients[N] = P->coefficients[p];
  }
  answer->num_terms = N;
  return true;
}

/** The basic row operation
 *
 * Assume that f is non-zero, and that the head term of f has nonzero
 * coefficient in g.  Kill that term of g by subtracting a scalar
 * multiple of f.  If the head coefficient of f is a and the
 * corresponding coefficient of g is b, then the computed answer is
 * g - (b/a)*f.
 *
 */

static inline bool row_op(Row_t *f, Row_t *g, Row_t *answer,
                          int g_coeff, int num_pivots, MConstants_t C) {
  int a = f->coefficients->value;
  int a_inverse = montgomery_inverse(a, C.prime, C.mu, C.R_cubed);
  /* Multiply by b and negate. Note that p - M(X) = M(p - X) = M(-X)u. */
  int factor = C.prime - montgomery_multiply(a_inverse, g_coeff, C.prime, C.mu);
  if (! Row_p_plus_aq(g, factor, f, answer, C)) {
    return false;
  }
  return true;
}

/** Static function to compare two Terms in two Rows.
 *
 * Compare the pth term of P to the qth term of Q and return an integer which is
 * < 0 if the first one is smaller, > 0 if the second one is smaller and 0 if
 * they are equal.  The 0 term is considered larger than any non-zero term.
 */

static inline int compare_heads(const void* p1, const void* p2) {
  Row_t* P1 = (Row_t*)p1;
  Row_t* P2 = (Row_t*)p2;
  return P1->coefficients->column_index - P2->coefficients->column_index;
}

static inline int compare_heads_dec(const void* p1, const void* p2) {
  Row_t* P1 = (Row_t*)p1;
  Row_t* P2 = (Row_t*)p2;
  return P2->coefficients->column_index - P1->coefficients->column_index;
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
 * column_index field of each coefficient.  This prepares the Polynomials for
 * conversion to Rows, which is then carried out.
 *
 * The number of columns in the matrix having the input polynomials as its rows
 * is stored in the int referenced by the num_columns input, and the table of
 * all terms is stored in the term_table input.
 */

static bool Poly_matrix_init(Polynomial_t** P, int num_rows, int* num_columns,
                             int* num_pivots,
                             Term_t** term_table, Row_t* matrix,
			     int rank, MConstants_t C) {
  int num_monomials = 0, i, j, index;
  monomial_t *pool;
  monomial_array_t monomial_arrays[num_rows], previous, merged;
  Term_t *table = NULL;
  char *pivot_info = NULL;
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
  if (NULL == (pivot_info = (char*)calloc(*num_columns, sizeof(char)))) {
    goto oom;
  }
  *num_pivots = 0;
  for (i = 0; i < num_rows; i++) {
    if (P[i]->num_terms > 0) {
      int column_index = P[i]->coefficients->column_index;
      if (pivot_info[column_index] == 0) {
        pivot_info[column_index] = 1;
        (*num_pivots)++;
      }
    }
  }
  /* Temporary */
  /* for (i = 0; i < *num_columns; i++) { */
  /*   printf("%d", pivot_info[i]); */
  /* } */
  /* printf(" (%d pivots %dx%d)\n", *num_pivots, num_rows, *num_columns); */
  /**/
  for (i = 0; i < num_rows; i++) {
    matrix[i] = zero_row;
    if (!Poly_compress(P[i], matrix + i, table, pivot_info, *num_pivots, C)) {
      for (j = 0; j < i; j++) {
        Row_free(matrix + i);
      }
      goto oom;
    }
  }
  free(pivot_info);
  return true;

 oom:
  free(merged.monomials);
  free(pool);
  free(table);
  free(pivot_info);
  return false;
}

/** Use bisection to find the coefficient of P with a given column index.
 *
 * Return false if no non-zero coefficient of P has the column index.
 * When the column indexes are available this is quite a bit faster.
 */

static inline bool coeff_in_column(Row_t* P, int column, int bottom,
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
  int i, j, coeff;
  Term_t* table = NULL;
  Row_t *row_i, *row_j, *matrix = NULL;
  Row_t buffer = zero_row, tmp;
  int head, num_pivots, last_pivot = -1;
  MConstants_t constants = montgomery_init(prime);
  
  /* Create the matrix. */
  if (NULL == (matrix = (Row_t*)malloc(num_rows*sizeof(Row_t)))) {
    goto oom;;
  }
  if (!Poly_matrix_init(P, num_rows, num_columns, &num_pivots, &table, 
                        matrix, rank, constants)) {
    goto oom;
  }
  
  /* Allocate one extra row to use as a buffer. */
  buffer.term_table = table;
  if (!Row_alloc(&buffer, *num_columns)) {
    goto oom;
  }
  
  /* This implementation requires sorting the rows by increasing head term.*/
  qsort(matrix, num_rows, sizeof(Row_t), compare_heads);

  /* Reduce the matrix to echelon form. */
  for (i = 0; i < num_rows; i++) {
    row_i = matrix + i;
    if (row_i->num_terms == 0) continue;
    head = row_i->coefficients->column_index;
    /* Clear above.  Since head terms are initially non-decreasing, we can skip
     * this step the first time that a new head term is seen because we know
     * that the column is already clear above.
     */
    if (head <= last_pivot) {
      for (j = 0; j < i; j++) {
        row_j = matrix + j;
        if (row_j->num_terms == 0) continue;
        if (coeff_in_column(row_j, head, 0, row_j->num_terms, &coeff)) {
          if (! row_op(row_i, row_j, &buffer, coeff, num_pivots, constants)) {
            return false;
          }
          tmp = matrix[j];
          matrix[j] = buffer;
          buffer = tmp;
        }
      }
    } else {
      SET_IS_PIVOT_ROW(row_i);
      last_pivot = head;
    }
    /* Clear below. */
    for (j = i + 1; j < num_rows; j++) {
      row_j = matrix + j;
      if (row_j->num_terms == 0) continue;
      if (coeff_in_column(row_j, head, 0, row_j->num_terms, &coeff)) {
        if (! row_op(row_i, row_j, &buffer, coeff, num_pivots, constants)) {
          return false;
        }
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
   * Free unneeded memory and convert the row polynomials to monic polynomials
   * of standard flavor.
   */
  for (i = 0; i < num_rows; i++) {
    row_i = matrix + i;
    if (row_i->num_terms == 0) {
      answer[i] = zero_poly;
      Row_free(row_i);
      continue;
    }
    Row_make_monic(row_i, constants);
    if (!Poly_decompress(row_i, answer + i, rank, constants)) {
      for (j = 0; j <= i; j++) {
	Row_free(matrix + j);
      }
      goto oom;
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