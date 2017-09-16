"""
 *  Some examples of endomorphism rings
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

load("../Initialize.sage")

F = QQ
R.<x> = PolynomialRing(F)
f, h = x^6 + 3*x^5 + 2*x^4 + 7*x^3 + 11*x^2 + 14, x^2 + x
f, h = -2*x^7 - 4*x^6 + 3*x^4 + x^3 - 2*x^2 - x, x^2 + x + 1
#f, h = x^7 - x^6 + 2*x^4 - 3*x^3 + 2*x^2 - x, x^4 + x^2 + 1
X = mHyperellipticCurve(f, h)

#R.<x,y,z> = PolynomialRing(F)
#f = x^3*z + x^2*y*z + x^2*z^2 + x*y^3 - 3*x*y^2*z - 4*x*z^3 - y^4 + 2*y^3*z + 2*z^4
#X = mPlaneCurve(f)

print ""
print "Curve:"
print X
Endo = EndomorphismData(X, prec = 300, have_oldenburg = True)

print ""
print "Lattice:"
print Endo.lattice().pretty_print()

Dec = Endo.decomposition()
print ""
print "Decomposition:"

# Genus 2 test:
while True:
    facs = Dec.factors()
    R.<x> = QQ[]
    d = 6
    pol = R(algdep(facs[0][1], d))
    if pol.factor()[0][0].degree() <= 3:
        print pol
#print Dec.verify()
