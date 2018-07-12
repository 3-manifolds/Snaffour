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
#cython: profile=False

from __future__ import print_function
from collections import Iterable
from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t

cdef extern from "F4.h":
    ctypedef int bool
    cdef struct Term_s:
        pass
    ctypedef Term_s Term_t
    cdef bool Term_equals(Term_t *t, Term_t *s)
    cdef int  Term_hash(Term_t *t)
    cdef int  Term_total_degree(Term_t *t, int rank)
    cdef bool Term_divides(Term_t *t, Term_t *s)
    cdef bool Term_divide(Term_t *t, Term_t *s, Term_t *answer)
    cdef int  Term_multiply(Term_t *t, Term_t *s, Term_t *answer)
    cdef int  Term_gcd(Term_t *t, Term_t *s, Term_t *answer)
    cdef int  Term_lcm(Term_t *t, Term_t *s, Term_t *answer)
    cdef int  Term_revlex_diff(Term_t *t, Term_t *s, int rank)
    cdef void Term_print(Term_t *t)

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

cdef extern from "Python.h":
    pass

from collections import Mapping
from itertools import chain
import re

cdef class PolyRing(object):
    """
    A polynomial ring over a field with a fixed set of named (string)
    variables. Currently the coefficient field is the finite field of order
    P = 2^31 - 1.

    A PolyRing is required when instantiating a Term, a Monomial, a Polynomial
    or an Ideal.  It specifies the number and names of variables and provides
    methods for normalizing the coefficients.  The best way to do this is to use
    the Term, Monomial, Polynomial or Ideal methods of the PolyRing.
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

    def Monomial(self, *args):
        return Monomial(*args, ring=self)

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

    A Term is modeled as a vector of up to 32 signed chars, stored in an array of
    two 16-byte vector types suitable for use in a 128-bit SSE2 MMX register.

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
    cdef public _hash
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
        self._hash = None
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
        # This could be much more efficient!
        if self._hash is None:
            self._hash = hash(self.degree)
        return self._hash

    def __eq__(self, Term other):
        if isinstance(self, Term) and isinstance(other, Term):
            return Term_equals(self.c_term, other.c_term) == 1
        return False

    def __lt__(self, Term other):
        if self.ring is other.ring:
            if self.total_degree < other.total_degree:
                return True
            if self.total_degree == other.total_degree:
                return Term_revlex_diff(
                    self.c_term, other.c_term, self.rank) > 0
            else:
                return False
        else:
            raise ValueError('Terms over different rings are incomparable.')

    def __floordiv__(Term self, Term other):
        cdef Term answer = Term(ring=self.ring)
        if Term_divide(self.c_term, other.c_term, answer.c_term) != 1:
            raise RuntimeError("Term division failed for %s // %s!"%(self, other))
        answer.total_degree = Term_total_degree(answer.c_term, self.ring.rank)
        return answer

    def __mul__(self, other):
        cdef Term left, right, answer
        if isinstance(self, Term) and isinstance(other, Term) and self.ring == other.ring:
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
        elif isinstance(self, Term) and isinstance(other, Monomial) and self.ring == other.ring:
            return Monomial(*((self*other.term).degree),
                            coefficient=other.coefficient,
                            ring=self.ring)
        elif isinstance(self, Monomial) and isinstance(other, Term) and self.ring == other.ring:
            return Monomial(*((self.term*other).degree),
                            coefficient=self.coefficient,
                            ring=other.ring)
        elif isinstance(self, Monomial) and isinstance(other, Monomial) and self.ring == other.ring:
            return Monomial(*((self.term*other.term).degree),
                            coefficient=self.ring.reduce_coeff(self.coefficient*other.coefficient),
                            ring=self.ring)
        else:
            return NotImplemented

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
    Base class for polynomials over a PolyRing.  A polynomial is modeled as an array of
    terms, sorted by decreasing degree, and an array of coefficients.

    Instantiate a polynomial with a mapping from degrees to coefficients or with
    any number of monomial arguments.

    Note that the F4 algorithm only requires adding polynomials and multiplying
    polynomials by monomials or terms.

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
    cdef _hash
    cdef public PolyRing ring
    cdef public is_nonzero
    cdef public Monomial head_monomial
    cdef public Term head_term

    def __cinit__(self):
        self.c_poly.num_terms = 0
        self.c_poly.terms = NULL
        self.c_poly.coefficients = NULL

    def __dealloc__(self):
        #print('Polynomial dealloc')
        Poly_free(&self.c_poly)

    def __init__(self, *args, PolyRing ring=no_ring):
        cdef int i, size
        self._hash = None
        self.ring = ring
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
        # This could surely be much more efficient!
        if self._hash is None:
            self._hash = 0
            for m in self.monomials:
                self._hash ^= hash(m)
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

    def terms(Polynomial self):
        cdef int i
        cdef Term a_Term
        result = []
        for i in range(self.c_poly.num_terms):
            a_Term = Term(ring=self.ring)
            a_Term.c_term[0] = self.c_poly.terms[i]
            a_Term.total_degree = Term_total_degree(a_Term.c_term, self.ring.rank)
            result.append(a_Term)
        return result

    @property
    def monomials(self):
        return [self.monomial(n) for n in range(self.c_poly.num_terms)]

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
        return '<%s, %s>'%(self.left_poly, self.right_poly)

    def spoly(self):
        ring = self.left_head.ring
        a = ring.invert_coeff(self.left_poly.head_monomial.coefficient)
        b = ring.invert_coeff(self.right_poly.head_monomial.coefficient)
        return (a*(self.lcm // self.left_head)*self.left_poly -
                b*(self.lcm // self.right_head)*self.right_poly)

    def is_reducing(self):
        f = self.spoly()
        if f.head_term < self.left_head or f.head_term < self.right_head:
           return True

cdef class Ideal(object):
    cdef public generators
    cdef public monic_generators
    cdef public ring
    cdef public _groebner_basis
    cdef public _reduced_groebner_basis
    cdef public echelons
    cdef public matrices
    cdef public select
    cdef public history

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
        self.history = []
        self._groebner_basis = None
        self._reduced_groebner_basis = None
        self.select = self.normal_select

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

    def echelon_form(self, poly_list):
        """
        Return a list of Polynomials in echelon form computed from a list of
        Polynomials.  The ideal generated by the echelon form is equal to that
        generated by the input polynomials.  To be in echelon form means that
        each polynomial is monic and that no other polynomial in the list has
        the same head term.

        >>> R = PolyRing('a', 'b', 'c', 'd')
        >>> f1 = Polynomial({(1,0,0,1): 1, (0,1,0,1): 1, (0,0,1,1): 1, (0,0,0,2): 1}, ring=R)
        >>> f2 = Polynomial({(1,1,0,0): 1, (0,1,1,0): 1, (1,0,0,1): 1, (0,0,1,1): 1}, ring=R)
        >>> f3 = Polynomial({(1,1,0,0): 1, (0,2,0,0): 1, (0,1,1,0): 1, (0,1,0,1): 1}, ring=R)
        >>> I = Ideal(ring=R)
        >>> I.echelon_form([f1, f2, f3])
        [a*d + b*d + c*d + d^2, a*b + b*c - b*d - d^2, b^2 + 2*b*d + d^2]
        >>> f3 = Polynomial({(1,1,0,0): 2, (0,2,0,0): 4, (0,1,1,0): 4, (0,1,0,1): 1}, ring=R)
        >>> f3
        2*a*b + 4*b^2 + 4*b*c + b*d
        >>> I.echelon_form([f1, f2, f3])
        [a*d + b*d + c*d + d^2, a*b + b*c - b*d - d^2, b^2 - 1073741823*b*c - 536870911*b*d - 1073741823*d^2]
        >>> P = 2**31 - 1
        >>> -2*1073741823 % P
        1
        >>> -4*536870911 % P
        3
        """
        cdef Polynomial_t** polys
        cdef Term_t** terms
        cdef Polynomial_t* answer
        cdef int N = len(poly_list)
        cdef int n
        cdef Polynomial p
        polys = <Polynomial_t **>malloc(N*sizeof(Polynomial_t*))
        answer = <Polynomial_t *>malloc(N*sizeof(Polynomial_t))
        for n, p in enumerate(poly_list):
            polys[n] = &p.c_poly
        if not Poly_echelon(polys, answer, self.ring.BIG_PRIME, self.ring.rank, N):
            raise RuntimeError('Out of memory')
        free(polys)
        result = []
        for n in range(N):
            f = Polynomial(ring=self.ring)
            f.c_poly = answer[n]
            f.decorate()
            if f.is_nonzero:
                result.append(f)
        free(answer)
        return result

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
        Use the F4 algorithm to find a grevlex Gröbner basis of the ideal
        generated by the polynomials in F.  NOTE: This does not compute
        a reduced groebner basis, and hence does not produce a canonical
        result.
        """
        if self._groebner_basis is not None:
            return self._groebner_basis
        self.make_monic_generators()
        G, P = [], set()
        for f in self.monic_generators:
            G, P = self.update(G, P, f)
        while P:
            self.history.append((list(G), list(P)))
            P_new = self.select(P)
            L = [(p.lcm // p.left_head, p.left_poly) for p in P_new]
            L += [(p.lcm // p.right_head, p.right_poly) for p in P_new]
            tilde_F_plus = self.reduce(L, G)
            P = P - P_new
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
            return {p}
    
    def normal_select(self, pairs):
        """
        The normal selector.
        """
        d = min(p.lcm.total_degree for p in pairs)
        return {p for p in pairs if p.lcm.total_degree == d}

    def reduce(self, L, G):
        r"""
        Start with a list L containing left and right elements from a selected
        set of critical pairs.  Use preprocessing and row echelon form to compute
        a set of polynomials to be added to the partial grobner basis G.

        This is the Reduction subalgorithm of Faugère's F4.

        INPUT:
          * L, a set of pairs (t, f), t a term and f a polynomial
          * G, a set of polynomials

        SIDE EFFECTS:
          self.rows is appended with the rows of the echelon form computed from the
          preprocessed list of polynomials.

        OUTPUT:
          * Returns :math:`\tilde F^plus` after saving the result F of preprocessing
            and :math:`\tilde F`.  (A key point of the algorithm is to use *all* rows
            of the echelon form.)

        """
        F = self.preprocess(L, G)
        # Simplify just iterates through the rows of F, so a list is fine.
        self.matrices.append(F)
        F_ech = self.echelon_form(F)
        # Sinplify looks up echelon rows by their unique head term.  So use a dict.
        self.echelons.append({f.head_term: f for f in F_ech})
        # Return the set that Faugère calls $\tilde F_d^+$.  It contains only
        # those rows of the echelon form whose head term is new, i.e. not the
        # head term of a polynomial in F.  These get added to G.
        heads = {f.head_term for f in F}
        return {f for f in F_ech if f.head_term not in heads}

    cdef heads_and_tails(self, F):
        """
        Return a list of all terms which appear in one of the Polynomials in F, but
        do not appear as a head term of any of those Polynomials.
        """
        heads = {f.head_term for f in F}
        terms = set.union(*(set(f.terms()) for f in F))
        return heads, terms - heads

    def preprocess(self, L, G):
        """
        First simplify the unevaluated products in L to get S.  Next, for every
        non-head term t of a polynomal in S, if t has a reduction modulo G then
        add a (simplified) reducer u*g to S, where HT(u*g) = t and g is a row
        from one of the echelon forms.  Return S.

        Adding these reducers means that the polynomials in S will get reduced modulo
        G during the computation of the echelon form.

        This is the Symbolic Preprocessing subalgorithm of Faugère's F4 algorithm.

        INPUT:
          * L, a set of unevaluated products (t, f), t a term and f a polynomial
          * G, a set of polynomials (to be improved to a Groebner basis.)

        OUTPUT:
          * F, a set of polynomials such that:
            o F contains a simplification of t*f for every (t, f) in L
            o For any non-head term T of a polynomial in F, if T is reducible
              modulo G then F contains a simplified reducing polynomial for T.
        """
        cdef Term s, g_head
        cdef Term s_over_ghead = Term(ring=self.ring)
        cdef reducer, heads, tails
        cdef int rank = self.ring.rank
        # Start by simplifying the unevaluated products in L.
        S = {self.mult(self.simplify(t, f)) for t, f in L}
        # Add (simplified) u*g, g in G, which reduce non-head terms of these products.
        heads, tails = self.heads_and_tails(S)
        while tails:
            s = tails.pop()
            for g in G:
                g_head = g.head_term
                # if HT(g) divides s
                if Term_divide(s.c_term, g_head.c_term, s_over_ghead.c_term):
                    # Make s_over_ghead into a valid Term by adding its total_degree attribute
                    s_over_ghead.total_degree = Term_total_degree(s_over_ghead.c_term, rank)
                    # Add the simplified reducer.
                    reducer = self.mult(self.simplify(s_over_ghead, g))
                    S.add(reducer)
                    heads.add(reducer.head_term)
                    tails.union({t for t in reducer.terms()[1:] if t not in heads})
                    break
        return S

    def update(self, G, P, h):
        """
        Use the two "Buchberger Criteria" to update a partial basis and a list of
        critical pairs to accommodate a new Polynomial h.  All pairs involving h
        are added to P, then useless pairs are removed.  Elements of G whose
        head term is divisible by h are removed, and finally h is added to G.

        See page 230 of "Gröbner Bases" by T. Becker and V. Weispfenning.
        """
        cdef Term p_lcm, q_lcm
        cdef Term h_head = h.head_term
        C, D, E, P_new = [Pair(g, h) for g in G if g != h], [], [], []
        # Discard (p, h) if its lcm is divisible by the lcm of another pair (q, h)
        # that we have already saved (unless the head terms of p and h are disjoint;
        # disjoint pairs are dealt with later). The 2nd Criterion implies
        # that (p, h) is useless in this case.
        while C:
            p, useless = C.pop(), False
            p_lcm = p.lcm
            if not p.is_disjoint:
                for q in chain(C, D):
                    q_lcm = q.lcm
                    if Term_divides(q_lcm.c_term, p_lcm.c_term):
                        useless = True
                        break
            if not useless:
                D.append(p)
        # Finally, discard pairs with disjoint head terms.  (1st Criterion)
        E = [p for p in D if not p.is_disjoint]
        # We don't need (f, g) if we have both (f, h) and (g, h).  (2nd Criterion)
        for p in P:
            p_lcm = p.lcm
            # If both lcm(f,h) and lcm(g,h) properly divide lcm(f,g) then they
            # are both in E, so we don't need to put (f,g) in P_new.
            if not (Term_divides(h_head.c_term, p_lcm.c_term) and
                    # This could be done in C
                    (h_head | p.left_head != p_lcm and
                     h_head | p.right_head != p_lcm)):
                P_new.append(p)
        # Add the pairs in E.
        P_new += E
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
        cdef Term f_head
        cdef Term t_over_u = Term(ring=self.ring)
        cdef Term_t u
        cdef Term_t uq_head
        cdef int td
        cdef int rank = self.ring.rank
        result = (t, q) # The default, if there is nothing better.
        if t.total_degree == 0:
            return result
        for F_ech, F in zip(self.echelons, self.matrices):
             for f in F:
                 f_head = f.head_term
                 # We want u, a divisor of t, such that u*HT(q) = HT(f)
                 if (Term_divide(f_head.c_term, q_head.c_term, &u) and
                     Term_divide(t.c_term, &u, t_over_u.c_term)):
                     # Make t_over_u a valid Term by adding its total_degree attribute.
                     t_over_u.total_degree = td = Term_total_degree(t_over_u.c_term, rank)
                     # We are guaranteed a unique p in F_ech with HT(f) = HT(p)
                     p = F_ech[f_head]
                     if Term_equals(t.c_term, &u):
                         # t // u == 1.  No further reduction is possible.
                         return (t_over_u, p)
                     elif Term_equals(t.c_term, t_over_u.c_term):
                         # t // u == t.  Avoid infinite recursion, but continue with the
                         # loop in case there is something better later on.
                         # Note that we cannot use t_over_u since it changes in the loop.
                         result = (t, p)
                         continue
                     else:
                         # Recursively search for a simpler pair.
                         return self.simplify(t_over_u, p)
        # Nothing more left to do.
        return result

    def normal_form(self, f):
        return self._normalize(f, self.groebner_basis())

    def _normalize(self, f, G):
        """
        If a term of f is divisible by a head term of g in G, then kill it by
        subtracting a multiple of g.  Return f when no further simplification is
        possible.  We assume that the elements of G are monic, since in practice
        G will usually be part of a Gröbner basis.
        """
        progress = True
        while progress:
            progress = False
            # One day we might want to do this in C
            for m in f.monomials:
                for g in G:
                    if g.head_term.divides(m.term):
                        f = f - m.coefficient*(m.term // g.head_term)*g
                        progress = True
                        break
        return f

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
        G = self.groebner_basis()
        result = []
        for n, g in enumerate(G):
            h = self._normalize(g, G[:n] + G[n+1:])
            if h.is_nonzero:
                result.append(h)
        result.sort(key=lambda f: f.head_term, reverse=True)
        self._reduced_groebner_basis = result
        return result
