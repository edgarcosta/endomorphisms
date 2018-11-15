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


function ReconstructionsFromComponentG1(P, Q, mor : ProjOrInc := "Proj")

A, R := Explode(mor);
E, h := ReconstructCurveG1(Q, K);
Anew := ConjugateMatrix(h, A);
return [ [* E, Q, [* Anew, R *] *] ];

end function;


function ReconstructionsFromComponentG2(P, Q, mor : ProjOrInc := "Proj")

gP := #Rows(P);
A, R := Explode(mor); K := BaseRing(A);
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
    Anew := ConjugateMatrix(h, A);
    Append(~facs, [* Y, Q, [* Anew, Rnew *] *]);
end for;
return facs;

end function;


intrinsic ReconstructionsFromComponent(P::., Q::., mor::. : ProjOrInc := "Proj") -> .
{Given a factor (Q, mor) of the Jacobian, finds corresponding curves along with morphisms to them.}

gQ := #Rows(Q);
if gQ eq 1 then
    return ReconstructionsFromComponentG1(P, Q, mor : ProjOrInc := ProjOrInc);
elif gQ eq 2 then
    return ReconstructionsFromComponentG2(P, Q, mor : ProjOrInc := ProjOrInc);
else
    error "Finding factors not yet implemented for genus larger than 2";
end if;

end intrinsic;
