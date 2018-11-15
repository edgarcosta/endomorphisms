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


function DecompositionFactorsG1(P, idem, K : ProjOrInc := "Proj")

Q, mor := ComponentFromIdempotent(P, idem : ProjOrInc := ProjOrInc);
E, h := ReconstructCurveG1(Q, K);
return [ [* E, h, Q, mor *] ];

end function;


function DecompositionFactorsG2(P, idem, K : ProjOrInc := "Proj")

g := #Rows(P);
Q, mor := ComponentFromIdempotent(P, idem : ProjOrInc := ProjOrInc);
A, R := Explode(mor);

EQ := InducedPolarization(StandardSymplecticMatrix(g), R : ProjOrInc := ProjOrInc);
Us := IsogenousPPLatticesG2(EQ);

facs := [ ];
for U in Us do
    if ProjOrInc eq "Proj" then
        Qnew := Q*ChangeRing(U^(-1), BaseRing(Q));
        Rnew := T*R;
    else
        Qnew := Q*ChangeRing(Transpose(U), BaseRing(Q));
        Rnew := Transpose(T)*R;
    end if;
    assert IsBigPeriodMatrix(Qnew);
    /* Reconstruct curves (some of them may give rise to an extension, but the
     * tangent representation is always the identity) */
    Y, h := ReconstructCurveG2(Q, K);
    Append(~facs, [* E, h, Q, [* A, Rnew *] *]);
end for;
return facs;

end function;


intrinsic DecompositionFactors(P::ModMatFldElt, idem::List, K::Fld : ProjOrInc := "Proj") -> .
{Finds curves corresponding to Prym variety of given idempotent.}

g := Rank(idem[2]) div 2;
if g eq 1 then
    return DecompositionFactorsG1(P, idem, K : ProjOrInc := ProjOrInc);
elif g eq 2 then
    return DecompositionFactorsG2(P, idem, K : ProjOrInc := ProjOrInc);
else
    error "Finding factors not yet implemented for genus larger than 2";
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
