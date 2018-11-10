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

Q, proj := ProjectionFromIdempotent(P, idem);
A, R := Explode(proj);
return [* ReconstructCurveG1(Q, F) *];

end function;


function DecompositionFactorsG2(P, idem, F)

Q, proj := ProjectionFromIdempotent(P, idem);
A, R := Explode(proj);

/* Induced polarization and overlattices needed to make it principal */
/* TODO: Take actual kernel and induced polarization on that, then same idea */
EQ := InducedPolarization(StandardSymplecticMatrix(3), R);
Us := IsogenousPPLatticesG2(EQ);

Qsnew := [ ];
for U in Us do
    Qnew := Q*ChangeRing(U^(-1), BaseRing(Q));
    assert IsBigPeriodMatrix(Qnew);
    Append(~Qsnew, Qnew);
end for;
recs := [* *];
for Q in Qsnew do
    Y := ReconstructCurveG2(Q, F);
    Append(~recs, Y);
end for;
return recs;

end function;


intrinsic DecompositionFactors(P::ModMatFldElt, idem::List, F::Fld) -> .
{Finds curves corresponding to Prym variety of given idempotent.}
/* TODO: So no image, but a kernel */

g := Rank(idem[2]) div 2;
if g eq 1 then
    return DecompositionFactorsG1(P, idem, F);
elif g eq 2 then
    return DecompositionFactorsG2(P, idem, F);
else
    error "Finding factors not implemented for genus larger than 2";
end if;

end intrinsic;


intrinsic ReconstructCurveFromRoot(root::.) -> .
{Curve reconstruction with extension if needed.}

Qroot, hcomp := Explode(root); K := BaseRing(hcomp[1]);
g := #Rows(Qroot);
if g eq 1 then
    return ReconstructCurveG1(Qroot, K);
elif g eq 2 then
    return ReconstructCurveG2(Qroot, K);
end if;
error "Reconstruction for genus larger than 2 not yet implemented";

end intrinsic;
