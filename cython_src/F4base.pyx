# -*- coding: utf-8 -*-
#   This file is part of the program Snaffour.
#
#   Copyright (C) 2018 by Marc Culler, Nathan Dunfield, Matthias Görner
#   and others.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#   Project homepage: https://bitbucket.org/t3m/snaffour/
#   Author homepage: https://marc-culler.info
#   Author homepage: http://dunfield.info
#   Author homepage: http://www.unhyperbolic.org/

# change the next line to use cProfile
#cython: profile=True

from __future__ import print_function
from collections import Iterable, Mapping
from itertools import combinations, chain
from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t
import re

cdef extern from "F4.h":
    ctypedef int bool
    cdef struct Term_s:
        pass
    ctypedef Term_s Term_t
    cdef bool Term_equals(Term_t *t, Term_t *s)
    cdef int  Term_total_degree(Term_t *t, int rank)
    cdef bool Term_divides(Term_t *t, Term_t *s)
    cdef bool Term_divide(Term_t *t, Term_t *s, Term_t *answer)
    cdef int  Term_multiply(Term_t *t, Term_t *s, Term_t *answer)
    cdef int  Term_gcd(Term_t *t, Term_t *s, Term_t *answer)
    cdef int  Term_lcm(Term_t *t, Term_t *s, Term_t *answer)
    cdef int  Term_revlex_diff(Term_t *t, Term_t *s, int rank)
    cdef void Term_print(Term_t *t)
    cdef long Term_hash(Term_t *t)
    cdef bool Term_merge(Term_t* s, Term_t* t, int s_size, int t_size,
                         Term_t** answer, int* answer_size, int rank)

    cdef int inverse_mod(int p, int x);

    cdef struct coeff_s:
        int total_degree
        int value
    ctypedef coeff_s coeff_t

    cdef struct Polynomial_s:
        int num_terms
        int rank
        coeff_t* coefficients
        Term_t* terms
    ctypedef Polynomial_s Polynomial_t
    cdef bool Poly_alloc(Polynomial_t* P, size_t size, int tank)
    cdef void Poly_free(Polynomial_t* P)
    cdef void Poly_print(Polynomial_t* P, int rank)
    cdef void Poly_init(Polynomial_t* P, size_t size, Term_t* terms, coeff_t* coefficients,
                        int rank)
    cdef void Poly_new_term(Polynomial_t* P, Term_t* term, coeff_t coefficient, int rank)
    cdef bool Poly_equals(Polynomial_t* P, Polynomial_t *Q)
    cdef bool Poly_add(Polynomial_t* P, Polynomial_t* Q, Polynomial_t* answer, int prime, int rank)
    cdef bool Poly_sub(Polynomial_t* P, Polynomial_t* Q, Polynomial_t* answer, int prime, int rank)
    cdef int  Poly_coeff(Polynomial_t* P, Term_t* t, int rank)
    cdef bool Poly_make_monic(Polynomial_t *P, Polynomial_t *answer, int prime, int rank)
    cdef bool Poly_echelon(Polynomial_t **P, Polynomial_t *answer, int prime, int rank,
                           size_t num_rows)
    cdef bool Poly_times_term(Polynomial_t *P, Term_t *t, Polynomial_t *answer, int prime, int rank)
    cdef bool Poly_times_int(Polynomial_t *P, int a, Polynomial_t *answer, int prime, int rank)
    cdef void Poly_sort(Polynomial_t *P, int num_polys, bool increasing)
    cdef bool Poly_terms(Polynomial_t *P, int num_polys, Term_t **answer, int* answer_size,
                         int rank);

cdef extern from "Python.h":
    pass

cdef class PolyRing(object):
    """
    A polynomial ring over a field with a fixed set of named (string)
    variables. Currently the coefficient field is the finite field of order
    P = 2^31 - 1.

    A PolyRing is required when instantiating a Term, a Monomial, a Polynomial
    or an Ideal.  It specifies the number and names of variables and provides
    methods for normalizing the coefficients.  The best way to instantiate a
    a Term, Monomial, Polynomial or Ideal is to use the eponymous method of
    a PolyRing.
    """
    cdef public int rank
    cdef public variables

    BIG_PRIME = 2**31 - 1
    CENTER = BIG_PRIME // 2

    def __init__(self, *variables):
        self.variables = tuple(variables)
        self.rank = len(self.variables)

    cpdef normalize_coeff(self, int c):
        """
        A normalized representative is in the interval [-P//2, P//2).
        This usually looks better when printing coefficients.
        """
        # We assume c is reduced mod BIG_PRIME so it fits in an int
        return c if c < self.CENTER else c - self.BIG_PRIME

    cpdef reduce_coeff(self, c):
        """Reduce mod P."""
        return c%self.BIG_PRIME

    cpdef negate_coeff(self, int c):
        """Negate mod P."""
        return self.BIG_PRIME - c

    cpdef invert_coeff(self, int c):
        """Compute the reciprocal mod P."""
        return inverse_mod(self.BIG_PRIME, c);

    def Term(self, *args):
        return Term(*args, ring=self)

    def Monomial(self, *args, **kwargs):
        return Monomial(*args, **kwargs, ring=self)

    def Polynomial(self, *args):
        return Polynomial(*args, ring=self)

    def Ideal(self, *args):
        return Ideal(*args, ring=self)

no_ring = PolyRing()

cdef class Term(object):
    """
    A Term over a PolyRing is a formal product of powers of the variables of the
    PolyRing.  Its degree is the tuple of powers.  Its total_degree is the sum
    of the powers.  Term degrees are compared using the graded reverse lexicographic
    ordering.

    The C structure underlying a Term is a vector of up to 32 signed chars,
    stored in an array of two 16-byte vector types suitable for use in a 128-bit
    SSE2 MMX register.

    >>> R = PolyRing('x', 'y', 'z')
    >>> t1 = Term(1,2,3,ring=R)
    >>> t2 = Term(2,2,1,ring=R)
    >>> t1, t2
    (x*y^2*z^3, x^2*y^2*z)
    >>> t1 * t2
    x^3*y^4*z^4
    >>> t1 | t2
    x^2*y^2*z^3
    >>> t1 ^ t2
    x*y^2*z
    >>> x = Term(1,0,0,ring=R)
    >>> y = Term(0,1,0,ring=R)
    >>> z = Term(0,0,1,ring=R)
    >>> x*x > x*y > y*y > x*z > y*z > z*z
    True
    """
    cdef Term_t *c_term
    cdef public int rank
    cdef public int total_degree
    cdef public ring

    def __cinit__(self):
        self.c_term = <Term_t*>malloc(sizeof(Term_t))

    def __dealloc__(self):
        #print('Term dealloc')
        if self.c_term is not NULL:
            free(<void*>self.c_term)

    def __init__(self, *degree, PolyRing ring=no_ring):
        cdef int i
        cdef char* powers = <char *>self.c_term
        cdef int64_t *int_powers = <int64_t*>self.c_term
        if ring is no_ring:
            raise ValueError('Please provide a PolyRing when instantiating a Term.')
        else:
            self.ring = ring
            self.rank = ring.rank
        self.total_degree = 0
        if degree:
            assert self.rank == len(degree) <= 32
            for i in range(self.rank):
                self.total_degree += degree[i]
                powers[i] = degree[i]
            for i in range(self.rank, 32):
                powers[i] = 0
        else:
            for i in range(4):
                int_powers[i] = 0

    def __repr__(self):
        powers = []
        for var, n in zip(self.ring.variables, self.degree):
            if n == 0:
                continue
            elif n == 1:
                powers.append(var)
            else:
                powers.append('%s^%d'%(var, n))
        if powers:
            return '*'.join(powers)
        else:
            return '<1>'

    def __hash__(self):
        return Term_hash(self.c_term)

    def __eq__(self, Term other):
        if isinstance(self, Term) and isinstance(other, Term):
            return Term_equals(self.c_term, other.c_term) == 1
        return False

    def __lt__(self, Term other):
        if self.ring is other.ring:
            if self.total_degree < other.total_degree:
                return True
            if self.total_degree == other.total_degree:
                return Term_revlex_diff(self.c_term, other.c_term, self.rank) > 0
            else:
                return False
        else:
            raise ValueError('Terms over different rings are incomparable.')

    def __le__(self, Term other):
        if self == other or self < other:
            return True
        return False

    def __floordiv__(Term self, Term other):
        cdef Term answer = Term(ring=self.ring)
        if Term_divide(self.c_term, other.c_term, answer.c_term) != 1:
            raise RuntimeError("Term division failed for %s // %s!"%(self, other))
        answer.total_degree = Term_total_degree(answer.c_term, self.ring.rank)
        return answer

    def __mul__(self, other):
        cdef Term left, right, answer
        if isinstance(self, Term) and isinstance(other, Term) and self.ring is other.ring:
            left, right = self, other
            answer = Term(ring=left.ring)
            Term_multiply(left.c_term, right.c_term, answer.c_term)
            answer.total_degree = Term_total_degree(answer.c_term, left.ring.rank)
            return answer
        elif isinstance(self, int) and isinstance(other, Term):
            coeff = other.ring.reduce_coeff(self)
            return Monomial(*other.degree, coefficient=coeff, ring=other.ring)
        else:
            return NotImplemented

    def __or__(Term self, Term other):
        """
        The least common multiple of two terms.
        """
        cdef Term answer = Term(ring=self.ring)
        Term_lcm(self.c_term, other.c_term, answer.c_term)
        answer.total_degree = Term_total_degree(answer.c_term, self.ring.rank)
        return answer

    def __xor__(Term self, Term other):
        """
        The greatest common divisor of two terms.
        """
        cdef Term answer = Term(ring=self.ring)
        Term_gcd(self.c_term, other.c_term, answer.c_term)
        answer.total_degree = Term_total_degree(answer.c_term, self.ring.rank)
        return answer

    @property
    def degree(self):
        """
        Return a tuple of integers representing the degree of this term.
        """
        cdef int i
        cdef char* powers = <char *>self.c_term
        return tuple([powers[i] for i in range(self.rank)])

    def divides(self, Term other):
        """
        Return True or False indicating whether this Term divides the other.
        """
        if isinstance(self, Term) and isinstance(other, Term):
            return Term_divides(self.c_term, other.c_term) == 1
        return False

    def mult(self, Polynomial other):
        """
        Return the product of this Term times a Polynomial.
        """
        result = Polynomial(ring=self.ring)
        if not Poly_times_term(&other.c_poly, self.c_term, &result.c_poly,
                        self.ring.BIG_PRIME, self.ring.rank):
            raise RuntimeError('Out of memory!')
        result.decorate()
        return result

cdef class Monomial(object):
    """
    A Monomial is just a Term with a coefficient.  Monomials are not designed to
    be particularly efficient, just convenient for instantiating or printing
    Polynomials.
    """
    cdef public int coefficient
    cdef public term
    cdef public ring

    def __init__(self, *degree, coefficient=1, PolyRing ring=no_ring):
        self.term = Term(*degree, ring=ring)
        self.ring = self.term.ring
        self.coefficient = ring.reduce_coeff(coefficient)

    def __hash__(self):
        return hash(self.coefficient) ^ hash(self.term)

    def __mul__(self, other):
        if isinstance(self, int) and isinstance(other, Monomial):
            coeff = other.ring.reduce_coeff(self*other.coefficient)
            return Monomial(*other.degree, coefficient=coeff, ring=other.ring)
        elif isinstance(self, Term) and isinstance(other, Monomial) and self.ring is other.ring:
            return Monomial(*((self*other.term).degree),
                            coefficient=other.coefficient,
                            ring=self.ring)
        elif isinstance(self, Monomial) and isinstance(other, Term) and self.ring is other.ring:
            return Monomial(*((self.term*other).degree),
                            coefficient=self.coefficient,
                            ring=other.ring)
        elif isinstance(self, Monomial) and isinstance(other, Monomial) and self.ring is other.ring:
            return Monomial(*((self.term*other.term).degree),
                            coefficient=self.ring.reduce_coeff(self.coefficient*other.coefficient),
                            ring=self.ring)
        else:
            return NotImplemented

    def __floordiv__(Monomial self, Term other):
        cdef Term t = Term(ring=self.term.ring)
        cdef Term s = self.term
        assert other.ring is self.ring
        if Term_divide(s.c_term, other.c_term, t.c_term) != 1:
            raise RuntimeError("Term division failed for %s // %s!"%(self, other))
        t.total_degree = Term_total_degree(t.c_term, self.ring.rank)
        return Monomial(*t.degree, coefficient=self.coefficient, ring=self.term.ring)

    def __neg__(self):
        return Monomial(*self.term.degree,
            coefficient=self.ring.negate_coeff(self.coefficient))

    def __eq__(self, other):
        return (self.coefficient == other.coefficient and
                self.term == other.term)

    def __repr__(self):
        return ('%s*%s'%(self.ring.normalize_coeff(self.coefficient), self.term))

cdef class Polynomial(object):
    """
    A polynomial over a PolyRing.  The underlying C strucure of a Polynomial is
    an array of terms, maintained sorted by decreasing degree, and a seperate
    array of coefficients.  (They are kept separate for alignment reasons.)

    Instantiate a polynomial with a mapping from degrees to coefficients or with
    any number of monomial arguments.

    Note that the F4 algorithm only requires adding polynomials and multiplying
    polynomials by scalars or terms.

    >>> R = PolyRing('a', 'b', 'c', 'd')
    >>> f = Polynomial({(1,1,0,0): 1, (0,1,1,0): 1, (0,0,1,1): 1, (1,0,0,1): 1}, ring=R)
    >>> g = Polynomial({(1,1,0,1): 1, (0,1,1,0): 1, (1,0,1,1): 1, (1,0,1,1): 1}, ring=R)
    >>> f
    a*b + b*c + a*d + c*d
    >>> g
    a*b*d + a*c*d + b*c
    >>> f + g
    a*b*d + a*c*d + a*b + 2*b*c + a*d + c*d
    >>> f - g
    -a*b*d - a*c*d + a*b + a*d + c*d
    >>> g - f
    a*b*d + a*c*d - a*b - a*d - c*d
    >>> f - f
    0
    >>> f + f
    2*a*b + 2*b*c + 2*a*d + 2*c*d
    >>> f = Polynomial({(1,0,0,1): 1, (0,1,0,1): 1, (0,0,1,1): 1, (0,0,0,2): 1}, ring=R)
    >>> g = Polynomial({(1,1,0,0): 1, (0,1,1,0): 1, (1,0,0,1): 1, (0,0,1,1): 1}, ring=R)
    >>> f
    a*d + b*d + c*d + d^2
    >>> g
    a*b + b*c + a*d + c*d
    >>> f -g
    -a*b - b*c + b*d + d^2
    >>> g- f
    a*b + b*c - b*d - d^2
    """
    cdef Polynomial_t c_poly
    cdef public PolyRing ring
    cdef public is_nonzero
    cdef public Monomial head_monomial
    cdef public Term head_term
    cdef public _hash

    def __cinit__(self):
        self.c_poly.num_terms = 0
        self.c_poly.terms = NULL
        self.c_poly.coefficients = NULL

    def __dealloc__(self):
        #print('Polynomial dealloc')
        Poly_free(&self.c_poly)

    def __init__(self, *args, PolyRing ring=no_ring):
        cdef int i, size
        self.ring = ring
        self._hash = None
        if not args:
            self.is_nonzero = False
            self.head_monomial = self.head_term = None
            return
        if isinstance(args[0], Mapping):
            coeff_map = args[0]
            monomials = [
                Monomial(*degree,
                         coefficient=self.ring.reduce_coeff(coeff_map[degree]),
                         ring=self.ring)
                for degree in coeff_map]
        else:
            monomials = list(args)
        self.is_nonzero = True
        monomials.sort(key=lambda m : m.term, reverse=True)
        size = len(monomials)
        if not Poly_alloc(&self.c_poly, size, self.ring.rank):
            raise RuntimeError('Out of memory!')
        for m in monomials:
            self.add_term(m)
        self.head_monomial = head = monomials[0]
        self.head_term = head.term

    def __getitem__(self, n):
        return self.monomial(n)

    def __repr__(self):
        if self.is_nonzero:
            result = ' + '.join([str(m) for m in self.monomials])
            result = re.sub(r'\*<1>', '', result)
            result = re.sub(r'\+ -', '- ', result)
            result = re.sub(r' 1\*',' ', result)
            result = re.sub(r'^(-{0,1})1\*', r'\1', result)
            return result
        else:
            return '0'

    def __len__(self):
        return self.c_poly.num_terms

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.monomials)
        return self._hash

    def __eq__(Polynomial self, Polynomial other):
        return Poly_equals(&self.c_poly, &other.c_poly) == 1

    def __add__(Polynomial self, Polynomial other):
        result = Polynomial(ring=self.ring)
        if not Poly_add(&self.c_poly, &other.c_poly, &result.c_poly,
                        self.ring.BIG_PRIME, self.ring.rank):
            raise RuntimeError('Out of memory!')
        result.decorate()
        return result

    def __sub__(Polynomial self, Polynomial other):
        result = Polynomial(ring=self.ring)
        if not Poly_sub(&self.c_poly, &other.c_poly, &result.c_poly,
                        self.ring.BIG_PRIME, self.ring.rank):
            raise RuntimeError('Out of memory!')
        result.decorate()
        return result

    def __mul__(self, Polynomial other):
        cdef Polynomial result = Polynomial(ring=other.ring)
        cdef Polynomial temp
        cdef Term t
        cdef int a
        if isinstance(self, int):
            success = Poly_times_int(&other.c_poly, self, &result.c_poly,
                                other.ring.BIG_PRIME, other.ring.rank)
        elif isinstance(self, Term):
            t = self
            success = Poly_times_term(&other.c_poly, t.c_term, &result.c_poly,
                                other.ring.BIG_PRIME, other.ring.rank)
        elif isinstance(self, Monomial):
            m = self
            t = m.term
            temp = Polynomial(ring=other.ring)
            success = Poly_times_term(&other.c_poly, t.c_term, &temp.c_poly,
                                other.ring.BIG_PRIME, other.ring.rank)
            if not success:
                raise RuntimeError('Out of memory!')
            a = m.coefficient
            success = Poly_times_int(&temp.c_poly, a, &result.c_poly,
                                other.ring.BIG_PRIME, other.ring.rank)
        else:
            return NotImplemented
        if not success:
            raise RuntimeError('Out of memory!')
        result.decorate()
        return result

    cdef decorate(self):
        """
        A Polynomial can be created from a C Polynomial_t by instantiating
        a zero poynomial and then overwriting the pointers in its c_poly.
        After doing so, you must call this method to make the other attributes
        consistent.
        """
        if self.c_poly.num_terms == 0:
            self.is_nonzero=False
            self.head_monomial = self.head_term = None
        else:
            self.is_nonzero=True
            head = self.monomial(0)
            self.head_monomial = head
            self.head_term = head.term

    cdef monomial(self, int n):
        """
        Return this Polynomial's n_th Monomial in reverse grevlex order.
        """
        if n > self.c_poly.num_terms:
            raise IndexError
        cdef char* powers = <char *>&self.c_poly.terms[n]
        cdef int coefficient = self.c_poly.coefficients[n].value
        cdef int i
        degree = [powers[i] for i in range(self.ring.rank)]
        return Monomial(*degree, coefficient=<int>coefficient, ring=self.ring)

    cdef add_term(self, Monomial m):
        cdef Term t = m.term
        cdef int N = self.c_poly.num_terms
        self.c_poly.terms[N] = t.c_term[0]
        self.c_poly.coefficients[N].total_degree = t.total_degree
        self.c_poly.coefficients[N].value = m.coefficient
        self.c_poly.num_terms += 1

    @property
    def monomials(self):
        return (self.monomial(n) for n in range(self.c_poly.num_terms))

    def coefficient(self, Term t):
        """
        Return the coefficient of the given term.  Terms which do not appear
        have coefficient 0.

        >>> R = PolyRing('a', 'b', 'c', 'd')
        >>> f = Polynomial({(1,1,0,1): 1, (0,1,1,0): 2, (1,0,1,1): 3, (0,0,0,0): 4}, ring=R)
        >>> f.coefficient(Term(0,0,0,0, ring=R))
        4
        >>> f.coefficient(Term(1,1,0,1, ring=R))
        1
        >>> f.coefficient(Term(0,1,1,1, ring=R))
        0
        >>> f.coefficient(Term(1,0,1,1, ring=R))
        3
        """
        return Poly_coeff(&self.c_poly, t.c_term, self.ring.rank)

cdef class Pair:
    """
    A "critical pair" consisting of two polynomials which can be combined to
    produce a candidate polynomial for a Gröbner basis.
    >>> R = PolyRing('a', 'b', 'c', 'd')
    >>> f = Polynomial({(1,1,0,0): 1, (0,1,1,0): 1, (0,0,1,1): 1, (1,0,0,1): 1}, ring=R)
    >>> g = Polynomial({(1,0,0,0): 1, (0,1,0,0): 1, (0,0,1,0): 1, (0,0,0,1): 1}, ring=R)
    >>> P = Pair(f, g)
    >>> P.left_poly, P.right_poly, P.lcm
    (a*b + b*c + a*d + c*d, a + b + c + d, a*b)
    """
    cdef public Polynomial left_poly, right_poly
    cdef public Term left_head, right_head, lcm
    cdef public int is_disjoint

    def __init__(self, Polynomial f, Polynomial g):
        cdef Term_t head_product
        assert f.ring is g.ring
        self.lcm = Term(ring=f.ring)
        self.left_poly, self.right_poly = f, g
        self.left_head  = f.head_term
        self.right_head = g.head_term
        Term_lcm(self.left_head.c_term, self.right_head.c_term, self.lcm.c_term)
        self.lcm.total_degree = Term_total_degree(self.lcm.c_term, f.ring.rank)
        Term_multiply(self.left_head.c_term, self.right_head.c_term, &head_product)
        self.is_disjoint = Term_equals(&head_product, self.lcm.c_term)

    def __repr__(self):
        return 'Pair:<%s ... , %s ...>'%(self.left_head, self.right_head)

    def spoly(self):
        ring = self.left_head.ring
        a = ring.invert_coeff(self.left_poly.head_monomial.coefficient)
        b = ring.invert_coeff(self.right_poly.head_monomial.coefficient)
        return (a*(self.lcm // self.left_head)*self.left_poly -
                b*(self.lcm // self.right_head)*self.right_poly)

    def left_prod(self):
        """
        Return (t, f) where f is the left poly and HT(t*f) = p.lcm.
        """
        return (self.lcm // self.left_head, self.left_poly)

    def right_prod(self):
        """
        Return (t, f) where f is the right poly and HT(t*f) = p.lcm.
        """
        return (self.lcm // self.right_head, self.right_poly)

cdef class PolyMatrix(object):
    """
    A PolyMatrix is initialized with a list of polynomials.  During initialization
    it computes and store the reduced echelon form in a way that will make it fast
    to access this data.  It does not hold references to the input polynomials.
    """
    cdef public PolyRing ring
    cdef public int num_rows
    cdef public tuple rows     # The Polynomials as rows of the echelon form
    cdef Term_t** c_heads      # The head terms, accessible from an array

    def __cinit__(self, poly_list):
        cdef int N = len(poly_list)
        self.c_heads = <Term_t**>malloc(N*sizeof(Term_t*))

    def __dealloc__(self):
        free(self.c_heads)

    def __init__(self, poly_list):
        cdef Polynomial_t** polys
        cdef Polynomial_t* answer
        cdef Polynomial_t c_poly
        cdef Term_t** terms
        cdef int N = len(poly_list)
        cdef int m, n
        cdef Polynomial p
        polys = <Polynomial_t **>malloc(N*sizeof(Polynomial_t*))
        answer = <Polynomial_t *>malloc(N*sizeof(Polynomial_t))
        assert isinstance(poly_list[0], Polynomial)
        self.ring = poly_list[0].ring
        for n, p in enumerate(poly_list):
            assert isinstance(p, Polynomial)
            assert p.ring is self.ring
            polys[n] = &p.c_poly
        if not Poly_echelon(polys, answer, self.ring.BIG_PRIME, self.ring.rank, N):
            raise RuntimeError('Out of memory')
        rows = []
        m = 0
        for n in range(N):
            c_poly = answer[n]  # Copies the Polynomial_t, so we can free answer
            if c_poly.num_terms > 0:
                self.c_heads[m] = c_poly.terms
                f = Polynomial(ring=self.ring)
                f.c_poly = c_poly
                f.decorate()
                rows.append(f)
                m += 1
        self.num_rows = m
        self.rows = tuple(rows)
        free(answer)
        free(polys)

cdef class F4State(object):
    """
    Records the state of the F4 algorithm at the beginning of each loop.
    """
    cdef public G
    cdef public P
    cdef public Sel
    cdef public S

    def __init__(self, G, P, Sel):
        self.G, self.P, self.Sel, self.S = list(G), set(P), set(Sel), []

cdef class Ideal(object):
    cdef public generators
    cdef public monic_generators
    cdef public ring
    cdef public _groebner_basis
    cdef public _reduced_groebner_basis
    cdef public echelons
    cdef public matrices
    cdef public select
#    cdef public history

    def __init__(self, *args, PolyRing ring=no_ring):
        if ring is no_ring:
            raise ValueError('A PolyRing must be specified.')
        if args and isinstance(args[0], Iterable):
            self.generators = tuple(args[0])
        else:
            self.generators = args
        self.ring = ring
        self.echelons = []
        self.matrices = []
        self._groebner_basis = None
        self._reduced_groebner_basis = None
        self.select = self.normal_select
#        self.history = []

    cdef mult(self, prod):
        """
        Evaluate an "unevaluated product" (t, f) with t a Term and f a Polynomial.
        """
        cdef Term t
        cdef Polynomial f
        t, f = prod
        result = Polynomial(ring=self.ring)
        if not Poly_times_term(&f.c_poly, t.c_term, &result.c_poly,
                        self.ring.BIG_PRIME, self.ring.rank):
            raise RuntimeError('Out of memory!')
        result.decorate()
        return result

    def set_select(self, name):
        if name == 'id':
            self.select = self.id_select
        elif name == 'buchberger':
            self.select = self.buchberger_select
        elif name == 'normal':
            self.select = self.normal_select
        else:
            raise ValueError('Unknown selector')

    @classmethod
    def is_groebner(cls, G, P=set(), details=False):
        """
        For each pair p = Pair(g1, g2) which is not in P, check whether
        the s-polynomial of p reduces to 0 modulo G.  If details == True,
        return the set of pairs which fail this test.  Otherwise, return
        a Boolean value True if G is Groebner, False if not.

        WARNING: This is very innefficient!
        """
        result = set()
        pairs = {Pair(g1, g2) for g1, g2 in combinations(G, 2)} - P
        for p in pairs:
            s = p.spoly()
            f = cls._normalize(s, G)
            if f.is_nonzero:
                result.add(p)
        if details:
            return result
        else:
            return len(result) == 0

    @classmethod
    def _normalize(self, f, G, return_rep=False):
        """
        If a term of f is divisible by a head term of g in G, then kill it by
        subtracting a multiple of g.  Return f when no further simplification is
        possible.  We assume in this private method that the elements of G are
        monic, since in our application G will be part of a Gröbner basis.
        """
        progress = True
        rep = []
        while progress:
            progress = False
            # One day we might want to make this faster, without the rep.
            for m in f.monomials:
                for i, g in enumerate(G):
                    if g.head_term.divides(m.term):
                        m1 = m // g.head_term
                        f = f - m1*g
                        rep.append((m1, i))
                        progress = True
                        break
        if return_rep:
            return f, rep
        else:
            return f


    cdef make_monic_generators(self):
        cdef Polynomial f, g
        self.monic_generators = []
        for f in self.generators:
             g = Polynomial(ring=self.ring)
             if not Poly_make_monic(&f.c_poly, &g.c_poly, self.ring.BIG_PRIME, self.ring.rank):
                 raise RuntimeError('Out of memory!')
             g.decorate()
             self.monic_generators.append(g)

    def groebner_basis(self):
        """
        Use the F4 algorithm to find a grevlex Gröbner basis of the ideal generated
        by the polynomials in F.  NOTE: This does not compute a reduced groebner
        basis, and hence does not produce a canonical result.
        """
        if self._groebner_basis is not None:
            return self._groebner_basis
        self.make_monic_generators()
        G, P = [], set()
#        self.history.append(F4State(G, P, set()))
        for f in self.monic_generators:
            G, P = self.update(G, P, f)
        while P:
            Sel = self.select(P)
#            self.history.append(F4State(G, P, Sel))
            L = ([p.left_prod() for p in Sel], [p.right_prod() for p in Sel])
            tilde_F_plus = self.reduce(L, G)
            P = P - Sel
            for h in tilde_F_plus:
                G, P = self.update(G, P, h)
        self._groebner_basis = G
        return G

    def id_select(self, pairs):
        """
        The identity selector.
        """
        return set(pairs)

    def buchberger_select(self, pairs):
        """
        Select a single pair.  This converts F4 into the Buchberger algorithm.
        """
        for p in pairs:
            return set((p,))
        # This variant works better, it seems.
        #s = sorted(pairs, key=lambda p: p.lcm.total_degree)
        #return set((s[0],))

    def normal_select(self, pairs):
        """
        The normal selector.
        """
        d = min(p.lcm.total_degree for p in pairs)
        return {p for p in pairs if p.lcm.total_degree == d}

    def reduce(self, L, G):
        r"""
        Start with a pair L of lists containing left and right elements from a
        selected set of critical pairs.  Use preprocessing and row echelon form
        to compute a set of polynomials to be added to the partial grobner basis
        G.

        This is the Reduction subalgorithm of Faugère's F4.

        INPUT:
          * L = (L1, L2), each a list of pairs (t, f), t a term and f a polynomial
          * G, a set of polynomials

        SIDE EFFECTS:
          The echelon form :math:`\tilde F is appended to self.echelons.

        OUTPUT:
          * Returns :math:`\tilde F^plus`, the rows of the echelon form of F
            whose head term was not a head term in F.  In particular, if (t1,g1)
            and (t2, g2) are in L, with HT(t1*g1) = HT(t2*g2) = Pair(g1,  g2).lcm,
            then F includes the (simplified) s-polynomial of Pair(g1, g2) and if this
            s-polynomial has a new head term then it will appear in the reuslt.

        """
        cdef tuple rows
        F = self.preprocess(L, G)
        F_ech = PolyMatrix(F)
        rows = F_ech.rows
        self.echelons.append(F_ech)
        heads = {f.head_term for f in F}
        return [f for f in rows if f.head_term not in heads]

    def terms(self, list polys):
        """
        Return a list of all terms that appear in the input list of Polynomials.
        The terms will be sorted by descending degree.
        """
        for P in polys:
            assert isinstance(P, Polynomial)
            assert P.ring is self.ring, (
                'The Polynomials must belong to the ring containing the ideal.')
        return self.term_list(polys)

    cdef term_list(self, poly_list):
        cdef Term_t *terms
        cdef int i, num_terms
        cdef int rank = self.ring.rank
        cdef Term t
        cdef Polynomial P
        cdef Polynomial_t* polys = <Polynomial_t*>malloc(
            len(poly_list)*sizeof(Polynomial_t))
        for i, P in enumerate(poly_list):
            polys[i] = P.c_poly
        if not Poly_terms(polys, len(poly_list), &terms, &num_terms, rank):
            raise RuntimeError('Out of memory!')
        result = []
        for i in range(num_terms):
            t = Term(ring=self.ring)
            t.c_term[0] = terms[i]
            t.total_degree = Term_total_degree(t.c_term, rank)
            result.append(t)
        free(terms)
        free(polys)
        return result

    cdef tails(self, list F):
        """
        Return a list of all terms which appear in one of the Polynomials in
        F, but do not appear as a head term of any of those Polynomials.
        The list is sorted by descending degree to make divisibility testing
        easier.
        """
        heads = {f.head_term for f in F}
        return tuple(t for t in self.term_list(F) if t not in heads)

    def preprocess(self, L, list G):
        """
        This method returns the matrix (as a list of polynomials) which will be
        reduced to echelon form in the main loop.  The input L consists of the
        left and right products (momomial times polynomial) for the selected
        pairs.  For each pair, the difference between these two products is the
        s-polynomial.  If the s-polynomial does not reduce to 0 modulo G_n then
        its reduction must be added to G_n when constructing G_n+1.  The
        reduction will appear as a row of the echelon form , provided that each
        product t*g which is used in reducing a term of the s-polynomial is
        included in the matrix.  It can be recognized because its head term will
        not appear among the head terms of the rows of the initial matrix.  The
        purpose of the preprocessing step is to add all such multiples of
        elements of G_n to the matrix.

        There is a serious, but easily fixable error in the pseudocode for the
        preprocessing subalgorithm in Faugère's paper, and this error can lead
        to incorrect computations.  He says to begin the preprocessing by
        replacing each product in the list L by its simplification.  He proves
        in Theorem 2.4 that the final G will be a Groebner basis if, for each
        distinct g1, g2 in G, the s-polynomial of their simplifications g1' and
        g2' has a t-representation with respect to G such that t < lcm(g1, g2).
        So it would appear to be sufficient that the simplifications of the
        products in L_n should all reduce to 0 modulo G.  The problem with this
        reasoning is that it ignores the degenerate case where g1 and g2 have
        the same simplification.  In that case the s-polynomial of g1' and g2'
        is zero, and the zero polynomial does not have t-representation by
        Definition 2.7. So Theorem 2.4 does not apply.  What happens in this
        case is that the matrix will have two identical rows, which appeared in
        a previous matrix, and the new head term which should have appeared in
        the s-polynomial of g1 and g2 will be missing from the final basis.  We
        observed this happening, and producing a non-Groebner basis, with the
        Cyclic-4 example when using the identity selection process that selects
        all pairs, as well as in some of the Buchberger selection processes that
        select one pair.

        We mention that Faugère's paper omits the proof of correctness of the F4
        algorithm.  Instead it (mis)states a Theorem from the book by T. Becker
        and V. Weispfenning along with Theorem 2.4, and says that these two results
        could be used in a proof of correctness, without actually giving the proof.

        This is the Symbolic Preprocessing subalgorithm of Faugère's F4 algorithm,
        modified to correct the error discussed above.

        INPUT:
          * L, a list of unevaluated left and right products from a set of pairs
          * G, a list of polynomials (to eventually be extended to a Groebner basis).

        OUTPUT:
          * F, a list of polynomials whose echelon form will contain the simplfication
            of the reduction modulo G of the s-polynomial of each pair, or the
            reduction of the equivalent s-polynomial of the simplified pair.
        """
        cdef Term t, t1, t2
        cdef Term t_over_ghead = Term(ring=self.ring)
        cdef Polynomial g
        cdef list S
        cdef tuple tails
        cdef Polynomial reducer, f1, f2
        cdef Term_t* g_head
        cdef int rank = self.ring.rank
        cdef int G_size = len(G)
        cdef int i, errors
        S = []
        errors = 0
        for p1, p2 in zip(*L):
             t1, f1 = self.simplify(p1[0], p1[1])
             t2, f2 = self.simplify(p2[0], p2[1])
             if (not Term_equals(t1.c_term, t2.c_term) or
                 not Term_equals(f1.c_poly.terms, f2.c_poly.terms)):
                 S.extend((self.mult((t1, f1)), self.mult((t2, f2))))
             else:
                 S.extend((self.mult(p1), self.mult(p2)))
        tails = self.tails(S)
        for t in tails:
            for g in G:
                g_head = g.c_poly.terms
                # if HT(g) divides t we add Simplify(t//g_head, g) to S
                if Term_divide(t.c_term, g_head, t_over_ghead.c_term):
                    t_over_ghead.total_degree = Term_total_degree(t_over_ghead.c_term, rank)
                    reducer = self.mult(self.simplify(t_over_ghead, g))
                    S.append(reducer)
                    break
#        self.history[-1].S = list(S)
        return S

    def update(self, G, P, h):
        """
        Use the Gebauer-Möller criteria to update a partial basis and a list of
        critical pairs to accommodate a new Polynomial h.  All pairs involving h
        are added to P, then useless pairs are removed.  Elements of G whose
        head term is divisible by h are removed, and finally h is added to G.

        Pseudocode for this subalgorithm is given on page 230 of "Gröbner Bases"
        by T. Becker and V. Weispfenning.
        """
        cdef Term p_lcm, q_lcm
        cdef Term h_head = h.head_term
        C = sorted((Pair(h, g) for g in G), key=lambda p: p.lcm.total_degree, reverse=True)
        D, E, P_new = [], [], []
        # Discard (h, g) if its lcm is divisible by the lcm of a saved or unseen
        # pair (h, f), provided that the head terms of f and g are not disjoint;
        # disjoint pairs are removed later in the hope that a disjoint pair will
        # be used to eliminate other pairs here before it gets deleted for being
        # disjoint below.
        while C:
            p = C.pop()
            if p.is_disjoint:
                D.append(p)
            else:
                p_lcm = p.lcm
                useless = False
                for q in chain(D, C):
                    q_lcm = q.lcm
                    if Term_divides(q_lcm.c_term, p_lcm.c_term):
                        useless = True
                        break
                    if p_lcm.total_degree < q_lcm.total_degree:
                        break
                if not useless:
                    D.append(p)
        # Now discard pairs with disjoint head terms.
        E = [p for p in D if not p.is_disjoint]
        Eg = {p.right_poly for p in E}
        # Add pairs (f, g) in P to P_new if they do not satiisfy:
        #   * HT(h) divides lcm(f, g) and (h, f) and (h, g) are both in E.
        P_new = [p for p in P if (not h_head.divides(p.lcm)) or
                 (p.left_poly not in Eg) or (p.right_poly not in Eg)]
        # Add the pairs in E to P_new.
        P_new.extend(E)
        # Remove redundant elements of G.
        G_new = [g for g in G if not h_head.divides(g.head_term)]
        G_new.append(h)
        return G_new, set(P_new)

    def simplify(self, Term t, Polynomial q):
        r"""
        The Simplify subalgorithm of F4.

        INPUT:
          * t, a term
          * q, a polynomial

        OUTPUT:
          * (u, g), where :math:`HT(u\star g) = HT(t\star q)` and (u,g) is simpler
            in the sense that the term u divides t, usually properly, and g is a row
            of a previously computed echelon form. The term u will be as small as
            possible.  This means that the new row :math:`u\star g` will be as close
            as possible to a previously computed row.

        NOTE: There are three serious typos in the description of this algorithm
        in Faugère's paper.  First he says:

        "there exists a (unique) :math:`p` in :math:`\tilde F^+_j` such that
        :math:`HT(P) = HT(u\star f)`".

        That is absurd since :math:`HT(u\star q)` is in :math:`HT(F_j)` and
        :math:`\tilde F_j^+` is constructed by removing all elements of
        :math:`\tilde F_j` having a head term
        in :math:`HT(F_j)`.  He means that :math:`p` is in :math:`\tilde F_j`.

        Second, he says: ":math:`u\star f` is in :math:`F_j`" when he means
        ":math:`HT(u\star f)` is in :math:`HT(F_j)`".

        The first typo is fixed in his 2013 slide presentation.  While the
        second typo is not fixed in the slides, his example shows that he is
        only looking at head terms (i.e. top reducibility).  As he explains in
        the slides, the goal of the simplification is to:

        "replace any product :math:`m\star f` by a product :math:`u\star t\star f'`
        where :math:`t\star f'` is a previously computed row and :math:`u\star t`
        divides the monomial :math:`m`"

        By "previously computed row" he means a row of of an echelon form
        :math:`\tilde F_j`.  The example suggests that he means that (t, f) should be
        replaced by (u, g) where g is a previously computed row and HT(t*f) =
        HT((u*g) with u < t.  The pseudocode in the paper generates (u, g) such
        that u|t and the proof of Lemma 2.3 assumes that u|t, but the statement
        does not assert it.  It is unclear whether he only needs the statement
        of the Lemma, of if he actually uses the proof at some point.

        Third, the algorithm, as presented, enters an infinite loop if passed a
        pair (t,f) with f in F and f in F_ech but t != 1.  In that case u = 1,
        since 1*f is in F, so u != t and f = p and the algorithm says to
        recursively call simplify(t/u, f) with t/u = t.  We need to handle this
        by stopping the recursion in this case.
        """
        cdef Term q_head = q.head_term
        cdef Term t_over_u = Term(ring=self.ring)
        cdef Polynomial f
        cdef PolyMatrix F_ech
        cdef list echelons = self.echelons
        cdef tuple rows
        cdef Term_t u
        cdef Term_t uq_head
        cdef Term_t **c_heads
        cdef int i
        cdef int rank = self.ring.rank
        result = (t, q) # The default result, if there is nothing better.
        if t.total_degree == 0:
            for F_ech in echelons:
                rows = F_ech.rows
                c_heads = F_ech.c_heads
                for i in range(F_ech.num_rows):
                    if Term_equals(q_head.c_term, c_heads[i]):
                        f = rows[i]
                        return(t, f)
        else:
            for F_ech in echelons:
                c_heads = F_ech.c_heads
                rows = F_ech.rows
                for i in range(F_ech.num_rows):
                    # Search for u, a divisor of t, such that u*HT(q) = HT(f)
                    if (Term_divide(c_heads[i], q_head.c_term, &u) and
                        Term_divide(t.c_term, &u, t_over_u.c_term)):
                        # Make t_over_u a valid Term by adding its total_degree attribute.
                        t_over_u.total_degree = Term_total_degree(t_over_u.c_term, rank)
                        f = rows[i]
                        if Term_equals(t.c_term, &u):
                            # t // u == 1:  No further reduction is possible.
                            return (t_over_u, f)
                        elif Term_equals(t.c_term, t_over_u.c_term):
                            # t // u == t:  Avoid infinite recursion, but continue to
                            # look for something better.
                            result = (t, f)
                            continue
                        else:
                            # Recursively search for a simpler pair.
                            return self.simplify(t_over_u, f)
        return result

    def normal_form(self, f):
        return self._normalize(f, self.groebner_basis())

    def reduced_groebner_basis(self):
        """
        Return the canonical reduced Gröbner basis of this ideal.
        >>> R3 = PolyRing('x', 'y', 'z')
        >>> A = [
        ... R3.Polynomial({(2,2,0):8, (1,3,0):5, (3,0,1):3, (2,1,1):1}),
        ... R3.Polynomial({(5,0,0):1, (0,3,2):2, (0,2,3):13, (0,1,4):5}),
        ... R3.Polynomial({(3,0,0):8, (0,3,0):12, (1,0,2):2, (0,0,0):3}),
        ... R3.Polynomial({(2,4,0):7, (1,3,2):18, (0,3,3):1})
        ... ]
        >>> A
        [8*x^2*y^2 + 5*x*y^3 + 3*x^3*z + x^2*y*z, x^5 + 2*y^3*z^2 + 13*y^2*z^3 + 5*y*z^4, 8*x^3 + 12*y^3 + 2*x*z^2 + 3, 7*x^2*y^4 + 18*x*y^3*z^2 + y^3*z^3]
        >>> I = R3.Ideal(A)
        >>> I.reduced_groebner_basis()
        [y^3 + 536870912, z^2, x]
        >>> R3.reduce_coeff(4*536870912)
        1
        """
        if self._reduced_groebner_basis is not None:
            return self._reduced_groebner_basis
        G = set(self.groebner_basis())
        H = set()
        while G:
            g = G.pop()
            keeper = True
            for f in chain(G, H):
                if f.head_term.divides(g.head_term):
                    keeper = False
                    break
            if keeper:
                H.add(g)
        R = sorted(H, key=lambda f: f.head_term, reverse=True)
        return [self._normalize(r, R[:n] + R[n+1:]) for n, r in enumerate(R)]
