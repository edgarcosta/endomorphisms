"""
 *  Some examples of bidimension and bidegree
"""

from endomorphisms import EndomorphismData, mHyperellipticCurve

F = QQ
R.<x> = PolynomialRing(F)
f = -x^5
h = x^3 + x + 1
#f = x^4 + x^3 + 3*x^2 + x + 2
#h = x^3 + x^2 + x
#f = x^5 - x
#h = x^3 + x + 1
#f = x^3 + x^2 - 2*x
#h = x^3 + x + 1
#f = x^3 - x
#h = x^3 + x^2 + 1

X = mHyperellipticCurve(f, h)
print ""
print "Curve:"
print X
Endo = EndomorphismData(X, prec = 100, have_oldenburg = True)

Dec = Endo.decomposition()
facs = Dec.factors()
idems = Dec.idempotents()
test = Dec.verify()

Dec = Endo.decomposition()
print ""
print "Decomposition:"
Dec = Endo.decomposition()
print "Factors:"
facs = Dec.factors()
print facs
print "Verify:"
test = Dec.verify()
print test
print Dec.correspondences()

print ""
print "Endomorphisms:"
print Endo.lattice().pretty_print()
print "Verify:"
test = Endo.verify()
print test

g = Endo.g
X = Endo.X
for dikt in Endo._geo_rep_dict_:
    print dikt['tangent']
    R = dikt['homology']
    S = dikt['corresp']
    if not magma.IsScalar(R):
        A = magma.Submatrix(R, 1,1,     g,g)
        B = magma.Submatrix(R, 1,g+1,   g,g)
        C = magma.Submatrix(R, g+1,1,   g,g)
        D = magma.Submatrix(R, g+1,g+1, g,g)
        print Endo.degree_estimate(dikt['tangent'])
        print magma.Trace(magma.Transpose(A)*D - magma.Transpose(C)*B)
        print S
        print magma.BiDimDeg(X, X, S)
        for I in magma.IrreducibleComponents(S):
            print magma.Dimension(I)
            print magma.BiDimDeg(X, X, I)
