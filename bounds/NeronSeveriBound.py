#  Copyright (C) 2016-2017 Edgar Costa
#  See LICENSE file for license details.

from sage.all import FiniteField, Infinity, Matrix, NumberField, PolynomialRing, QQ, cyclotomic_polynomial, sage_eval, sqrt
from sage.all import HyperellipticCurve, ZZ, GF, next_prime

def is_ordinary(a):
    assert len(a) == 5
    p = sqrt(a[-1]);
    return not ((p**2).divides(a[2]) or p.divides(a[1]))

def is_totally_split(L, p):
    if L.discriminant() % p == 0:
        return False;
    return len(PolynomialRing(FiniteField(p),'z')(L.absolute_polynomial().list()).factor()) == L.absolute_polynomial().degree()

def RM_and_notCM(ef, RM_coeff):
    QQt = PolynomialRing(QQ, 't')
    RM_field = NumberField(QQt(RM_coeff),'r')
    possible_CM = None;
    for p1 in ef:
        if len(p1) == 5:
            p = sqrt(p1[-1])
            if  QQt(p1).is_irreducible() and is_ordinary(p1) and is_totally_split(RM_field, p):
                N = NumberField(QQt(p1), "a")
                if possible_CM is None:
                    possible_CM = N;
                elif not N.is_isomorphic(possible_CM):
                    return True;
                else:
                    pass;
    return False

def cyclotomic_polynomial_bound_4():
    return [cyclotomic_polynomial(k) for k in [1, 2, 3, 4, 5, 6, 8, 10, 12] ]

cyclo_bound_4 = map( list, cyclotomic_polynomial_bound_4())


def rank_from_endo(string_list):
    """
    string_list given by the following process:
    X = mHyperellipticCurve(f,g)
    Endo = EndomorphismData(X, prec = 300, have_oldenburg = False)
    overK = Endo.geometric()
    print overK._desc_[2]
    """
    rank = 0
    for s in string_list:
        # in the LMFDB the description is stored without a space between M_2 and (
        if s in ['M_2 (CC)','M_2(CC)']:
            rank += 4
        elif s in ['M_2 (RR)', 'M_2(RR)']:
            rank += 3
        elif s == 'RR' or s =='CC':
            rank += 1
    return rank

def P2_4factor(a):
    # ref: Edgar's thesis section 2.3
    assert len(a) == 5
    q = sqrt(a[-1])
    b = [0]*5;
    b[0] = 1
    b[1] = 2*q - a[2]
    b[2] = (2*q ** 2 + a[1]**2 * q - 2*a[2]*q)
    b[3] = q**2*b[1]
    b[4] = q**4
    return q, b

def rank_disc(P1):
    R= PolynomialRing(QQ, 'T')
    T = R.gen();
    if len(P1) != 5:
        return +Infinity, None
        
    p, p2 =  P2_4factor(P1)
    
    cp = R(p2)(T/p)
    rank = 2
    
    rank = 2;
    t_cp = cp
    n_cp = R(1)
    for z in cyclo_bound_4:
        q, r = t_cp.quo_rem( R(z) )
        while r == 0:
            
            n_cp *= R(z);
            t_cp = q
            rank += len(z) -1 
            q, r = t_cp.quo_rem(R(z))
            
    assert n_cp * t_cp == cp
    e = 1
    assert rank == n_cp.degree() + 2
    while (T-1)**(rank - 2) != characteristic_polynomial_extension(n_cp, e):
        e += 1
    
    disc = (-1) * characteristic_polynomial_extension(t_cp, e)(1) * p**e
    assert rank % 2 == 0
    return p, rank, disc.squarefree_part()


def characteristic_polynomial_extension(cp, r):
    cplist = list(cp)
    assert cplist[-1] == 1
    
    n = len(cplist) -1
    
    A = Matrix(cp.base_ring(), n, n)
    for i in range(n-1):
        A[i+1,i] = 1
    for i in range(n):
        A[i, n-1] = -cplist[i]
        
    A = A**r
    cpnewlist = list(A.charpoly())
    cpnew =  cp.parent()(cpnewlist)    
    return cpnew


def upperrank(ef, verbose = False):
    rd = [];
    rank = +Infinity
    disc = None
    for p1 in ef:
        if len(p1) == 5:
            p, r, d = rank_disc(p1)
            # print p, r, d
            rd.append( [p, r, d] )
            if r < rank:
                rank = r;
                disc = d;
            elif r == rank:
                if disc != d:
                    rank -= 1
    if verbose:
        print rd
    return rank, disc

def verify_curve(g, factorsRR_geom, bound = 300, RM_coeff = None):
    D = ZZ(g.discriminant());
    C = HyperellipticCurve(g)
    p = 3;
    ef = [];
    while p < bound:
        if D.mod(p) != ZZ(0):
            ef +=  [C.change_ring(GF(p)).frobenius_polynomial().list()]
            ef[-1].reverse()
        p = next_prime(p)
    rank_euler, _ =  upperrank(ef)
    rank_endo = rank_from_endo(factorsRR_geom)
    if factorsRR_geom ==  ['RR', 'RR'] and RM_coeff is not None:
        checkRM =  RM_and_notCM(ef, RM_coeff)
        return rank_endo == rank_euler and checkRM, rank_endo, rank_euler
    else:
        return rank_endo == rank_euler, rank_endo, rank_euler

    


def verify_curve_lmfdb(label, Lhash = None):
    from lmfdb import getDBconnection
    import pymongo
    C = getDBconnection()
    curve = C.genus2_curves.curves.find_one({"label" : label})
    endo_query = C.genus2_curves.endomorphisms.find_one({'label': label})
    rank_endo = rank_from_endo(endo_query['factorsRR_geom'])
    
    if not endo_query['is_simple_geom']:
        assert rank_endo > 1

    
    if Lhash is None:
        Lhash = curve['Lhash'];

    Lfunction_query = C. Lfunctions.Lfunctions.find_one({'Lhash': Lhash})
    ef = sage_eval(str(Lfunction_query['euler_factors']))
    rank_euler, _ =  upperrank(ef)

    if endo_query['factorsRR_geom'] ==  [u'RR', u'RR'] and curve['is_simple_geom']:
        RM_coeff = endo_query['factorsQQ_geom'][0][1]
        checkRM =  RM_and_notCM(ef, RM_coeff)
        return rank_endo == rank_euler and checkRM, rank_endo, rank_euler
    
    return rank_endo == rank_euler, rank_endo, rank_euler

def verify_batch(start = 0, limitN = 7000, verbose = False):
    from lmfdb import getDBconnection
    C = getDBconnection()
    i = 0 
    bound = 0;
    label = None
    for curve in C.genus2_curves.curves.find().sort([("cond", pymongo.ASCENDING), ("label", pymongo.ASCENDING)]).limit(limitN).skip(start):
        label = curve['label']
        Lhash = curve['Lhash'];
        q, rendo, reuler = verify_curve_lmfd_lmfdb(label, Lhash)
        if not q:
            print "FAILED at label = %s" % label 
        if verbose:
            print label, q
        i+=1
        if int(100.0*i/limitN) >= bound:
            print "%s%%\t %s / 66158\t at label = %s" %(int(100.0*i/limitN), start + i, label)
            bound+=1
    print "Done from %s to %s / 66158\t at label = %s" %(start + 1, start + i, label)



