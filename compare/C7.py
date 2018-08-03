from snaffour.examples import *
I = R7.Ideal(C7)
G = I.reduced_groebner_basis()
print(len(G))
