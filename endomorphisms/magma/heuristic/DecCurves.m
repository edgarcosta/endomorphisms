/***
 *  Find curves in Jacobian
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


function DecompositionFactorsG1(P, idem, F)

/* Projection map */
Q, proj := ProjectionFromIdempotent(P, idem);
A, R := Explode(proj);
return [* ReconstructCurveG1(Q, F) *];

end function;


function DecompositionFactorsG2(P, idem, F)

/* Projection map */
Q, proj := ProjectionFromIdempotent(P, idem);
A, R := Explode(proj);

/* Induced polarization and overlattices needed to make them principal */
EQ := InducedPolarization(StandardSymplecticMatrix(3), R);
Us := IsogenousPPLatticesG2(EQ);

Qsnew := [ ];
for U in Us do
    Qnew := Q*ChangeRing(U^(-1), BaseRing(Q));
    assert IsBigPeriodMatrix(Qnew);
    Append(~Qsnew, Qnew);
end for;
return [* ReconstructCurveG2(Q, F) : Q in Qsnew *];

end function;


intrinsic DecompositionFactors(P::ModMatFldElt, idem::List, F::Fld) -> .
{Finds curves corresponding to given Jacobian factor.}

g := Rank(idem[2]) div 2;
if g eq 1 then
    return DecompositionFactorsG1(P, idem, F);
elif g eq 2 then
    return DecompositionFactorsG2(P, idem, F);
else
    error "Finding factors not implemented for given genus";
end if;

end intrinsic;
