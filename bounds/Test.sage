"""
 *  Bound functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

"""
load("constants.sage")
load("DiscriminantBound.sage")
load("EndomorphismRankBound.sage")
attach("Genus2Factors.sage")
load("GeometricallyIrreducible.sage")
#load("MagmaInterface.m")
load("NonCM.sage")
load("NonIsogenous.sage")
load("NonQM.sage")
load("PointCounting.sage")
load("TwistPolynomials.sage")
"""


def testHigherGenus() :
    R.<u> = PolynomialRing(QQ)

    # Reducible example
    f = u^6-7*u^4+14*u^2-7
    h = 0

    # Trivial endomorphisms in genus 3
    f = u^7-2*u^6+5*u^5-2*u^4+u^3+4*u^2+2*u
    h = u^2+1

    # RM by Q(sqrt(5))
    f = R([0, -2, 1])
    h = R([1, 1, 0, 1])


    # RM by Q(\sqrt{17})
    f = R([-3, 8, 5, 7, 2, 1])
    h = 0

    # Stupid full CM
    f = u^5 + 1
    f = u^7 + 1
    h = 0

    C = HyperellipticCurve(f,h)
    LPolys = ComputeLPolys(C)

    type, bd = DiscriminantBound(LPolys)

    print "Type", type
    print "Bound", bd.factor()

    print "Bound on the Z-rank of the endomorphism ring", EndomorphismRankBound(LPolys, C.genus())

def TestEC() :
    LPolys1 = [ 0 for i in range(0,maxP) ]
    LPolys2 = [ 0 for i in range(0,maxP) ]
    # p1 = [0, -1, 1, 0, 0]
    # p2 = [0, -1, 1, -10, -20]

    p1 = [0, -1, 0, 1, 0]
    p2 = [0, 1, 0, 1, 0]

    E1 = EllipticCurve(p1)
    E2 = EllipticCurve(p2)

    LPolys1 = ComputeLPolys(E1)
    LPolys2 = ComputeLPolys(E2)


    #d = E1.discriminant() * E2.discriminant()
    #for p in range(2,maxP) :
    #    if is_prime(p) and d%p != 0 :
    #        E1p = E1.base_extend(FiniteField(p))
    #        E2p = E2.base_extend(FiniteField(p))
    #        LPolys1[p] = E1p.frobenius_polynomial()
    #        LPolys2[p] = E2p.frobenius_polynomial()
    #print "Finished computing L-functions"
    print "Can prove there is no isogeny over the ground field?", CertifyNonIsogenous(LPolys1, LPolys2, false)
    print "Can prove there is no isogeny geometrically?", CertifyNonIsogenous(LPolys1, LPolys2)


def TestPicard() :
    R.<x> = PolynomialRing(QQ)
    LPolys = magma.computeLPolys2(100)
    print "Finished computing L-polynomials"
    LPolysInternal = [R(l) for l in LPolys]
    LPolysInternal.insert(0,0)
    print "Finished converting L-polys to Sage"
    type, bound = DiscriminantBound(LPolysInternal)

    print "Type", type
    print "Discriminant bound", bound

def TestDetectGenus2Factor() :
    R.<x> = PolynomialRing(QQ)

    f = x^7+x^6+x^5+x^3+x^2+x
    h = x^4+x^2+1
    # d = 3993
    # Cannot prove anything, which is correct because the quotient abelian surface is QM


    f = x^5-x^4+x^3
    h = x^4+1
    # d = 7744
    # OK, confirms that the quotient abelian surface has no extra endomorphisms (it used to: why did this stop working?)

    # f = x^8+x^2+1
    # h = 0


    C = HyperellipticCurve(f,h)
    LPolys = ComputeLPolys(C)


    Genus2LPolys = Genus2FactorTwistedLPolys(LPolys)

    type, bound = DiscriminantBound(Genus2LPolys, true)
    print "Type", type
    print "Discriminant bound", bound
