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

def CurveRankBound ( C ) :
    End = EndomorphismData(C, 100, have_oldenburg = true)
    genus = End.g

    if genus > 2 :
        raise ValueError("The genus is too large")

    LPolys = ComputeLPolys(C)
    Dec = End.decomposition()

    Irreducible = IsGeometricallyIrreducible(LPolys)
    ProvenReducible = len(Dec.idempotents()) > 1

    if ProvenReducible and Irreducible:
	raise ValueError("The Jacobian is geometrically irreducible, but a projection to an elliptic curve has been found (?)")

    if Irreducible :						# this covers ZZ, RM, CM
        Type, DiscBound = DiscriminantBound( LPolys )
        if Type == "Z" :
            return 1;
        if Type == "Quadratic" :
            return 2;
        if Type == "RM" :
            return genus;
        if Type == "FullCM" :
            return 2*genus;

    if ProvenReducible :
	return RankBoundProductEC(LPolys);

    # if we got here, the lower bound part of the computation *should* have proven that we have QM.
    # We just need to (1) check that (how?) and
    # (2) prove that we are not the square of an elliptic curve with CM.

    if IsGeometricallyNotSquareOfCM( LPolys ) :
	    return 4;
    else :
	    return 8;


    raise NotImplementedError()

