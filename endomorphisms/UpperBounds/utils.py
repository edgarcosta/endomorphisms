
#  Copyright (C) 2018, Edgar Costa
#  See LICENSE file for license details.


from sage.all import GCD, GF, HyperellipticCurve, NumberField, PolynomialRing, QQ
from sage.all import next_prime, pari

# Counting roots of polynomials

def count_cyclotomic_roots(f):
    r"""
    INPUT:

    -   ``f`` -- a polynomial with its roots in the unit circle

    OUTPUT: The number of roots that are cyclotomic
    """

    total = 0
    for factor, power in f.change_ring(QQ).factor():
        if factor.is_cyclotomic():
            total += power *factor.degree()
    return total

# get the polredabs polynomial
def polredabs(f):
    return f.parent()(pari(f).polredabs().list())


# field intersection

def field_intersection(L, K, Lsubfields = None):
    if L.is_isomorphic(K):
        return L
    Ksubfields = K.subfields();
    if Lsubfields is None:
        Lsubfields = L.subfields()

    # discarding QQ form the list
    for sL, _, _ in reversed(Lsubfields[1:]):
        for sK, _, _ in reversed(Ksubfields[1:]):
            if sL.is_isomorphic(sK):
                # sL and sK is the largest common subfield
                return NumberField(polredabs(sL.polynomial()),'a')
    else:
        # The intersection is QQ
        return NumberField(PolynomialRing(QQ,'x').gen(),'a')

    return L

def field_intersection_list(defining_polynomials):
    # try to just use discriminants of the first 1000 polynomials
    D = 0;
    for f in defining_polynomials[:1000]:
        D = GCD(D, f.discriminant())
        if D == 1:
            return NumberField(PolynomialRing(QQ,'x').gen(),'a')

    L = NumberField(defining_polynomials[0], 'a');
    for f in defining_polynomials:
        K = NumberField(f, 'b');
        # the 1st argument should be the 'smaller' field
        L = field_intersection(L, K)

        if L.degree() == 1:
            return NumberField(PolynomialRing(QQ,'x').gen(),'a');
    else:
        return NumberField(polredabs(L.absolute_polynomial()),'a')


# figuring out RR representation of the endomorphism algebra
def RR_representation(g, K, d):
    r"""
    g = dim of abelian variety (simple or isogenous to a power)
    K = center
    d = dimension of End over center
    """
    if K.is_CM():
        n = K.degree()/2
        KRR = 'CC'
    else:
        KRR = 'RR'
        n = K.degree();
        if d % 2 == 0 and g % 2 == 0 and g > 3:
            # Can't have type III for g = 2
            #We can deduce if it is type II or III
            KRR = 'M_2(RR) or HH'
            d = d / 2;
    if d > 1:
        out = 'M_%d(%s)' % (d, KRR,)
    else:
        out = KRR;

    return [out]*int(n)

# flatten a list
def flatten_list(l, level = None):
    if level == 0:
        return l
    if level is not None:
        level -= 1
    listtype = type([]);
    if type(l) == listtype:
        out = [];
        for elt in l:
            if type(elt) == listtype:
                out += flatten_list(elt, level)
            else:
                out += [elt]
        return out
    else:
        return l;



def get_frob_list_HyperellipticCurve(f, h, bound = 100):
    C = HyperellipticCurve(f,h)
    frob_list = []
    p = 2;
    while p < bound:
        try:
            f = C.change_ring(GF(p)).frobenius_polynomial()
            if f.degree() == 2*C.genus():
                frob_list.append([p, f.reverse()])
        except:
            pass
        p = next_prime(p)
    return frob_list



