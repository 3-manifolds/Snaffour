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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include <immintrin.h>

/*
gcc -o termtest -mavx2 -flax-vector-conversions -O3 termAVX.c
*/

/** A Term32_t is an array of 32 chars representing exponents.
 *
 * Assigning the vector_size attribute enables gcc to use the 256-bit AVX
 * registers to do arithmetic operations on all 32 bytes in a single
 * instruction.  This is meant to optimize operations such as multiplying
 * or dividing Terms.
 */

typedef char Term32_t __attribute__ ((vector_size (32)));

/** Terms
 *
 * We use a Term32_t vector to allow up to 32 variables.  For computing
 * Ptolemy varieties, 16 variables would not be enough.
 *
 * Wrapping the Term32_t in a struct makes Cython easier to deal with.
 */
typedef struct Term_s {
  Term32_t degree;
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

void Term_print(Term_t *t, int rank) {
  char *c = (char *)&t->degree;
  printf("[ ");
  for (int i=0; i<rank; i++) {
    printf("%hhd ", *c++);
  }
  printf("] ");
}


/* Borrowed from Python's tuple hash function. */
long Term_hash(Term_t *t)
{
  uint64_t* w = (uint64_t*)&t->degree;
  uint64_t x;
  uint64_t mult = 1000003UL  /* 0xf4243 */;
  x = 0x345678UL;
  for (int i=0; i<4; i++) {
        x = (x ^ w[i]) * mult;
        mult += 82528UL;
  }
  x += 97531UL;
  if (x == ~0) {
    x = -2;
  }
  return (long)x;
}

static inline bool is_zero_vector(Term32_t T) {
  Term32_t ones = (ones == ones);
  if (!_mm256_testz_si256((__m256i)T, ones)) {
    return true;
  }
  return false;
}

/* Does t equal s? */
bool Term_equals(Term_t *t, Term_t *s) {
  Term32_t result = (t->degree - s->degree);
  return !is_zero_vector(result);
}

/* stack overflow suggestion
  //  Term32_t ones = _mm256_set1_epi8(-1);


uint16_t sum_32(const uint8_t a[32])
{
    __m128i zero = _mm_xor_si128(zero, zero);
    __m128i sum0 = _mm_sad_epu8(
                        zero,
                        _mm_load_si128(reinterpret_cast<const __m128i*>(a)));
    __m128i sum1 = _mm_sad_epu8(
                        zero,
                        _mm_load_si128(reinterpret_cast<const __m128i*>(&a[16])));
    __m128i sum2 = _mm_add_epi16(sum0, sum1);
    __m128i totalsum = _mm_add_epi16(sum2, _mm_shuffle_epi32(sum2, 2));
    return totalsum.m128i_u16[0];
}
*/

int Term_total_degree(Term_t *t, int rank) {
  char *c = (char *)&t->degree;
  int i, degree = 0;
  for (i = 0; i < rank; i++) {
    degree += c[i];
  }
  return degree;
}

/* Does t divide s? */
bool Term_divides(Term_t *t, Term_t *s) {
    /* Return false if any power in t is larger than the corresponding
     * power in s
    */
  Term32_t result = (t->degree > s->degree);
  int64_t *L = (int64_t *) &result;
  if ( L[0] != 0 || L[1] != 0 || L[2] != 0 || L[3] != 0) {
    return false;
  }
  return true;
}

/* Divide t by s, if possible, and report success. */
bool Term_divide(Term_t *t, Term_t *s, Term_t *answer) {
  Term32_t failure;
  answer->degree = t->degree - s->degree;
  failure = (answer->degree < 0);
  int64_t *L = (int64_t *) &failure;
  if ( L[0] != 0 || L[1] != 0 || L[2] != 0 || L[3] != 0) {
    return false;
  }
  return true;
}

/* Multiply t by s */
void Term_multiply(Term_t *t, Term_t *s, Term_t *answer) {
    answer->degree = t->degree + s->degree;
}

/* Compute the least common multiple of t and s. */
void Term_lcm(Term_t *t, Term_t *s, Term_t *answer) {
  Term32_t t_bigger = (t->degree >= s->degree);
  Term32_t s_bigger = (t->degree < s->degree);
  answer->degree = (t->degree & t_bigger) | (s->degree & s_bigger);
}

/* Compute the greatest common divisor of t and s. */
void Term_gcd(Term_t *t, Term_t *s, Term_t *answer) {
  for (int i = 0; i < 2; i++) {
    Term32_t t_smaller = (t->degree <= s->degree);
    Term32_t s_smaller = (t->degree > s->degree);
    answer->degree = (t->degree & t_smaller) | (s->degree & s_smaller);
  }
}

/*
 * Return the difference T[i] - S[i] where i is the first term where they differ,
 * or return 0 if the terms are equal.  Thus a positive value means that t < s,
 * since the lex order is *reversed* in grevlex!
 */

int Term_revlex_diff(Term_t *t, Term_t *s, int rank) {
  int i = rank;
  Term32_t termdiff;
  char* diff = (char*)&termdiff;
  termdiff = t->degree - s->degree;
  while (i >= 0 && diff[i] == 0) {i--;}
  if (i >= 0) {
    return diff[i];
  }
  return 0;
}

/* Merge two strictly decreasing arrays of Terms
 *
 * Duplicates are removed and the resulting strictly decreasing list is stored
 * in answer.  Memory will be allocated for the answer. The caller is responsible
 * for freeing the memory.
 */

bool Term_merge(Term_t* s, Term_t* t, int s_size, int t_size,
                       Term_t** answer, int* answer_size, int rank) {
  int size = s_size + t_size, p = 0, q = 0, N = 0, s_td, t_td;
  Term_t* merged = (Term_t*)aligned_alloc(32, sizeof(Term_t)*size);
  if (merged == NULL) {
    *answer_size = 0;
    *answer = NULL;
    return false;
  }
  while (p < s_size && q < t_size) {
    s_td = Term_total_degree(s + p, rank);
    t_td = Term_total_degree(t + q, rank);
    int td_cmp = s_td - t_td;
    int revlex_cmp = td_cmp == 0 ? Term_revlex_diff(t + q, s + p, rank) : 0;
    if (td_cmp > 0 || revlex_cmp > 0) {
      /* s[p] > t[q] */
      merged[N++] = s[p++];
    } else if (td_cmp < 0 || revlex_cmp < 0) {
      /* t[q] > s[p]*/
      merged[N++] = t[q++];
    } else {
      /* s[p] == t[q] */
      merged[N++] = s[p++];
      q++;
    }
  }
  /* At most one of these two loops will be non-trivial. */
  while (p < s_size) {
    merged[N++] = s[p++];
  }
  while (q < t_size) {
    merged[N++] = t[q++];
  }
  *answer_size = N;
  *answer = merged;
  return true;
}

int main(int argc, char **argv) {
  Term_t X, Y, Z;
  int rank = 5;
  Term32_t a = {1,2,3,4,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Term32_t b = {1,2,3,4,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Term32_t c = {1,2,3,4,127,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  X.degree = a;
  Y.degree = b;
  Z.degree = c;
  Term_print(&X, rank);
  Term_print(&Y, rank);
  Term_print(&Z, rank);
  printf("\nX %s Y\n", Term_equals(&X, &Y)? "==" : "!=");
  printf("\nX %s Z\n", Term_equals(&X, &Z)? "==" : "!=");
  
  return 0;
}
