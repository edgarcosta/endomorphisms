
#  Copyright (C) 2018, Edgar Costa
#  See LICENSE file for license details.


from sage.all import GCD, GF, HyperellipticCurve, NumberField, QQ, Set
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
        return L.absolute_polynomial();
    Ksubfields = K.subfields()[1:];
    if Lsubfields is None:
        Lsubfields = L.subfields()

    Lsubfields = Lsubfields[1:];
    Kupper = len(Ksubfields);
    # discarding QQ form the list
    for sL, _, _ in reversed(Lsubfields):
        for j in reversed(range(Kupper)):
            sK = Ksubfields[j][0];
            if sK.absolute_degree() == sL.absolute_degree():
                if sL.is_isomorphic(sK):
                    # sL and sK is the largest common subfield
                    return polredabs(sL.absolute_polynomial())
            elif sK.absolute_degree() <  sL.absolute_degree():
                break;
            else:
                Kupper -= 1;
    else:
        # The intersection is QQ
        return sL.absolute_polynomial().parent().gen()

    return L

def field_intersection_list(defining_polynomials):
    # try to just use discriminants of the first 1000 polynomials
    D = 0;
    for f in defining_polynomials[:1000]:
        D = GCD(D, f.discriminant())
        if D == 1:
            return f.parent().gen()

    L = NumberField(defining_polynomials[0], 'a');
    for f in defining_polynomials:
        K = NumberField(f, 'b');
        # the 1st argument should be the 'smaller' field
        L = NumberField(field_intersection(L, K), 'a');

        if L.degree() == 1:
            return f.parent().gen();
    else:
        return polredabs(L.absolute_polynomial())

def field_intersection_matrix(matrix_defining_polynomials):
    r"""
    INPUT:

        - a matrix of polynomials defining numberfields (f_ij)

    OUTPUT:

        - returns a pair per column entry (Ak, Bk)
            after picking a row with a small intersection, WLOG i = 0
            Let S_ij denote all the polynomials defining subfields of QQ[x]/fij

            Then Bk = S_0k \cap_{i > 0) \cup_j S_ij
            If the list Bk can be obtained as list of subfields of one subfield
            then such subfield, otherwise Ak is None


    """
    if len(matrix_defining_polynomials[0]) == 1:
        f = field_intersection_list( [ elt[0] for elt in matrix_defining_polynomials ] );
        return [[f, subfields_polynomials(f)]]


    # we would like to pick the row which has the smallest intersection
    # an easy to get something close to that,is by picking the one with 
    # lowest GCD of the discriminants
    Di = None;
    Dmin = None;
    for i, row in enumerate(matrix_defining_polynomials):
        D = GCD([elt.discriminant() for elt in row])
        if Dmin is None or D < Dmin:
            Dmin = D;
            Di = i;
        if Dmin == 1:
            break;

    working_row = matrix_defining_polynomials[Di];

    subfields_matrix = [None] * len(matrix_defining_polynomials)
    output = [None] * len(working_row);
    for k, f in enumerate(working_row):
        Lsubfields = Set(subfields_polynomials(f));
        for i, row in enumerate(matrix_defining_polynomials):
            if i != Di:
                if subfields_matrix[i] is None:
                    U = Set([]);
                    for g in row:
                        Sg = Set(subfields_polynomials(g));
                        U = U.union(Sg);
                    subfields_matrix[i] = U;
                # try to discard some subfields
                Lsubfields = Lsubfields.intersection(subfields_matrix[i]);
                if Lsubfields.cardinality() == 1:
                    # QQ is the only possible subfield
                    break;

        Lsubfields = list(Lsubfields);
        Lsubfields.sort();
        output[k] = [None, Lsubfields];
        # try to produce a single polynomial, if Lsubfields are all the subfields of a common subfield
        if len(Lsubfields) == 1:
            output[k][0] = f.parent().gen();
        else:
            if Lsubfields[-1].degree() != Lsubfields[-2].degree():
                # there is one field of maximal degree, that might contain all the others
                K = NumberField(Lsubfields[-1],'a');
                maxsubfields = [polredabs(sK.absolute_polynomial()) for sK, _, _ in K.subfields()];
                if maxsubfields == Lsubfields:
                    output[k][0] = Lsubfields[-1]

    return output;


def subfields_polynomials(f):
    return [polredabs(sL.absolute_polynomial()) for sL, _, _ in NumberField(f,'a').subfields()]




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



