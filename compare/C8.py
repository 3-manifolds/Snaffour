from snaffour.examples import *
I = R8.Ideal(C8)
G = I.reduced_groebner_basis()
print(len(G))
