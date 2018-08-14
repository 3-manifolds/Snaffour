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
  unsigned long* w = (unsigned long*)&t->degree;
  long x;
  long mult = 1000003UL  /* 0xf4243 */;
  x = 0x345678UL;
  for (int i=0; i<4; i++) {
        x = (x ^ w[i]) * mult;
        mult += 82528UL;
  }
  x += 97531UL;
  if (x == ~0) {
    x = -2;
  }
  return x;
}

/* Does t equal s? */
bool Term_equals(Term_t *t, Term_t *s) {
  for (int i = 0; i < 2; i++) {
    Term16_t result = (t->degree[i] - s->degree[i]);
    int64_t *L = (int64_t *) &result;
    if ( L[0] != 0 || L[1] != 0) {
      return false;
    }
  }
  return true;
}

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
  for (int i = 0; i < 2; i++) {
    Term16_t result = (t->degree[i] > s->degree[i]);
    int64_t *L = (int64_t *) &result;
    if ( L[0] != 0 || L[1] != 0) {
      return false;
    }
  }
  return true;
}

/* Divide t by s, if possible, and report success. */
bool Term_divide(Term_t *t, Term_t *s, Term_t *answer) {
  Term16_t failure;
  for (int i = 0; i < 2; i++) {
    answer->degree[i] = t->degree[i] - s->degree[i];
    failure = (answer->degree[i] < 0);
    int64_t *L = (int64_t *) &failure;
    if ( L[0] != 0 || L[1] != 0) {
      return false;
    }
  }
  return true;
}

/* Multiply t by s */
void Term_multiply(Term_t *t, Term_t *s, Term_t *answer) {
  for (int i = 0; i < 2; i++) {
    answer->degree[i] = t->degree[i] + s->degree[i];
  }
}

/* Compute the least common multiple of t and s. */
void Term_lcm(Term_t *t, Term_t *s, Term_t *answer) {
  for (int i = 0; i < 2; i++) {
    Term16_t t_bigger = (t->degree[i] >= s->degree[i]);
    Term16_t s_bigger = (t->degree[i] < s->degree[i]);
    answer->degree[i] = (t->degree[i] & t_bigger) | (s->degree[i] & s_bigger);
  }
}

/* Compute the greatest common divisor of t and s. */
void Term_gcd(Term_t *t, Term_t *s, Term_t *answer) {
  for (int i = 0; i < 2; i++) {
    Term16_t t_smaller = (t->degree[i] <= s->degree[i]);
    Term16_t s_smaller = (t->degree[i] > s->degree[i]);
    answer->degree[i] = (t->degree[i] & t_smaller) | (s->degree[i] & s_smaller);
  }
}

/*
 * Return the difference T[i] - S[i] where i is the first term where they differ,
 * or return 0 if the terms are equal.  Thus a positive value means that t < s,
 * since the lex order is *reversed* in grevlex!
 */

int Term_revlex_diff(Term_t *t, Term_t *s, int rank) {
  int i = rank;
  Term16_t termdiff;
  char* diff = (char*)&termdiff;
  if (i > 16) {
    termdiff = t->degree[1] - s->degree[1];
    while (i >= 16 && diff[i - 16] == 0) {i--;}
    if (i >= 16) {
      return diff[i - 16];
    }
  }
  termdiff = t->degree[0] - s->degree[0];
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
  Term_t* merged = (Term_t*)malloc(sizeof(Term_t)*size);
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
