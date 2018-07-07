#include "F4.h"

void Term_print(Term_t *t, int rank) {
  char *c = (char *)&t->degree;
  printf("[ ");
  for (int i=0; i<rank; i++) {
    printf("%hhd ", *c++);
  }
  printf("] ");
}

int Term_hash(Term_t *t) {
  int *w = (int *)&t->degree;
  int hash = 0;
  for (int i=0; i<8; i++) {
    hash ^= w[i];
  }
  return hash;
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
  int degree = 0;
  for (int i=0; i<rank; i++) {
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
    failure = answer->degree[i] < 0;
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
 * or return 0 if the terms are equal.  Thus a positive value means that t < s.
 */
int Term_revlex_diff(Term_t *t, Term_t *s, int rank) {
  int i = rank;
  char *T = (char*)&t->degree;
  char *S = (char*)&s->degree;
  while (i >= 0 && T[i] == S[i]) {
      i--;
    }
  if (i < 0) {
    return 0;  /* t == s */
  }
  return T[i] - S[i];
}
