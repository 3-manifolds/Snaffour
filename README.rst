Snaffour
========

Description
-----------

Snaffour is an implementation of Faugère's F4 algorithm for computing
Gröbner bases of ideals in polynomial rings.  This implementation is
for polynomials over a finite field of prime order, where currently
the order is always 2^31 - 1.

Snaffour is a Python extension module, written in Cython, with low
level computations implemented directly in C.

This is an alpha release and is not ready to be used.

To build the extension module::

  python3 setup.py build

To run a quick test::

  python3 setup.py test

To build the documentation::

  python3 setup.py document

To clean up the build directory and start over::

  python3 setup.py clean
