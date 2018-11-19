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


/* TODO: Make geometric versions for use with period matrices only */

function ReconstructionsFromComponentG1(P, Q, mor : ProjOrInc := "Proj")

A, R := Explode(mor); K := BaseRing(A);
Rnew := R;
if Im(Q[1,2]/Q[1,1]) lt 0 then
    Q := Matrix([ [ Q[1,2], Q[1,1] ] ]);
    T := Matrix(Rationals(), [[0, 1],[1, 0]]);
    if ProjOrInc eq "Proj" then
        Rnew := T*R;
    else
        Rnew := R*T;
    end if;
end if;
assert IsBigPeriodMatrix(Q);
E, h := ReconstructCurve(Q, K : Base := true);
Anew := ConjugateMatrix(h, A);
return [ [* E, [* Anew, Rnew *] *] ];

end function;


function ReconstructionsFromComponentG2(P, Q, mor : ProjOrInc := "Proj")

gP := #Rows(P);
A, R := Explode(mor); K := BaseRing(A);
EQ := InducedPolarization(StandardSymplecticMatrix(gP), R : ProjOrInc := ProjOrInc);
E0, _ := FrobeniusFormAlternatingAlt(EQ);
vprint CurveRec: "Frobenius form:";
vprint CurveRec: E0;
Ts := IsogenousPPLattices(EQ);

facs := [ ];
for T in Ts do
    if ProjOrInc eq "Proj" then
        Qnew := Q*ChangeRing(T^(-1), BaseRing(Q));
        Rnew := T*R;
    else
        Qnew := Q*ChangeRing(Transpose(T), BaseRing(Q));
        Rnew := R*Transpose(T);
    end if;
    assert IsBigPeriodMatrix(Qnew);
    /* Reconstruct curves (some of them may give rise to an extension, but the
     * tangent representation is always the identity) */
    Y, h := ReconstructCurve(Qnew, K);
    vprint CurveRec: "";
    vprint CurveRec: "Reconstructed curve:";
    vprint CurveRec: Y;
    vprint CurveRec: "";
    Anew := ConjugateMatrix(h, A);
    Append(~facs, [* Y, [* Anew, Rnew *] *]);
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
