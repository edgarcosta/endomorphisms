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

def DiscriminantBound(LPolys, alreadyTwisted = false ):
    # returns a pair (type, bound), where type is one of "Z", "RM",
    # "QuadraticCM", "FullCM" and bound is a bound on the discriminant of the
    # field of endomorphisms
    discBound=0;

    CMField = Rationals()
    RMField = Rationals()
    FullCM = true
    RM = true

    for p in range (2,maxP):
      if LPolys[p] <> 0 :
        q = LPolys[p]
        g = q.degree() / 2
        extDegree = extensionBounds[g]
        if alreadyTwisted :
            extDegree = 1
        break

    for p in range (2,maxP):
      if is_prime(p) and LPolys[p] <> 0 :
        q = LPolys[p]
        # print q
        # print "Testing p =", p

        # print "Twisting the L-poly"
        q = twistPolynomial(q, extDegree)
        # print "Prime", p, "Polynomial", q;
        # print "Twisting ends, now test for irreducibility"
        if(q.is_irreducible()) :
            # print "Irreducible polynomial"
            if(q.coefficients(sparse=false)[g] % p<>0) :
                # print "Ordinary"
                # print "Define K1"
                K1=NumberField(q,'b')
                # print K1
                # print K1.discriminant()
                # print "Define K1Real"
                K1Real = K1.maximal_totally_real_subfield()[0];

                # print "Comparisons"
                if (CMField == Rationals()) :
                    CMField = K1
                if (RMField == Rationals()) :
                    RMField = K1Real
                if (FullCM and not CMField.is_isomorphic(K1) ) :
                    FullCM = false ;
                if ( RM and RMField.degree() == 1 ) :
                    RM = false ;
                if (RM and not RMField.is_isomorphic(K1Real) ) :
                    RM = false ;

                FieldDiscriminant=K1.discriminant()
                discBound=GCD(discBound,FieldDiscriminant)
                # print "Current discriminant bound", discBound
                if discBound<=4 :
                    return "Z", 1
    if FullCM :
        return "FullCM", discBound
    if RM :
        return "RM", discBound^(1/g)

    return "Quadratic", discBound^(1/g)
