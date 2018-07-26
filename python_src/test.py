from __future__ import print_function
from doctest import testmod
from . import *
results = testmod(F4base)
print('Doctest results:', results)

from .examples import *
print('\nComputing some Groebner bases over Fp, p = 2^31 - 1\n')
print("Elizabeth Arnold's example:")
print(R3.Ideal(A).reduced_groebner_basis())
print()
print('Cyclic-4 with normal selector:')
print(R4.Ideal(C4).reduced_groebner_basis())
print()
print('Cyclic-4 with identity selector:')
I = R4.Ideal(C4)
I.set_select('id')
GI = I.reduced_groebner_basis()
print(GI)
print()
print('Cyclic-4 with Buchberger selector:')
J = R4.Ideal(C4)
J.set_select('buchberger')
GJ = J.reduced_groebner_basis()
print(GJ)
print()
print('The size of the Cyclic-5 basis:')
G = R5.Ideal(C5).reduced_groebner_basis()
print(len(G))
