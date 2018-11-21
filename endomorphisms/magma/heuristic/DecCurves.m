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

function ReconstructionsFromComponentG1(P, Q, mor : ProjToIdem := true, ProjToPP := true)

A, R := Explode(mor); K := BaseRing(A);
Rnew := R;
if Im(Q[1,2]/Q[1,1]) lt 0 then
    Q := Matrix([ [ Q[1,2], Q[1,1] ] ]);
    T := Matrix(Rationals(), [[0, 1],[1, 0]]);
    CorrectMap := ProjToIdem eq ProjToPP;
    Rnew := R;
    if CorrectMap then
        if ProjToPP then
            Rnew := T*R;
        else
            Rnew := R*T;
        end if;
    end if;
end if;
assert IsBigPeriodMatrix(Q);
E, h := ReconstructCurve(Q, K : Base := true);
Anew := ConjugateMatrix(h, A);
return [ [* E, [* Anew, Rnew *] *] ];

end function;


function ReconstructionsFromComponentG2(P, Q, mor : ProjToPP := true, ProjToIdem := true)

gP := #Rows(P);
A, R := Explode(mor); K := BaseRing(A);
EQ := InducedPolarization(StandardSymplecticMatrix(gP), R : ProjToIdem := ProjToIdem);
E0, _ := FrobeniusFormAlternatingAlt(EQ);

vprint CurveRec: "Frobenius form of induced polarization:";
vprint CurveRec: E0;

Ts := IsogenousPPLattices(EQ);
facs := [ ];
for T in Ts do
    CorrectMap := ProjToIdem eq ProjToPP;
    Qnew := Q*ChangeRing(T^(-1), BaseRing(Q));
    Rnew := R;
    if CorrectMap then
        if ProjToPP then
            Rnew := T*R;
        else
            Rnew := R*T^(-1);
        end if;
    end if;
    assert IsBigPeriodMatrix(Qnew);
    Y, h := ReconstructCurve(Qnew, K);

    vprint CurveRec: "";
    vprint CurveRec: "Reconstructed curve found!";
    vprint CurveRec: Y;
    vprint CurveRec: "";

    Anew := ConjugateMatrix(h, A);
    Append(~facs, [* Y, [* Anew, Rnew *] *]);
end for;
return facs;

end function;


intrinsic ReconstructionsFromComponent(P::., Q::., mor::. : ProjToPP := true, ProjToIdem := true) -> .
{Given a factor (Q, mor) of the Jacobian, finds corresponding curves along with morphisms to them.}

gQ := #Rows(Q);
vprint CurveRec : "";
vprint CurveRec : "Reconstructing curves...";
if gQ eq 1 then
    recs := ReconstructionsFromComponentG1(P, Q, mor : ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
    vprint CurveRec : "done reconstructing curves.";
    return recs;
elif gQ eq 2 then
    recs := ReconstructionsFromComponentG2(P, Q, mor : ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
    vprint CurveRec : "done reconstructing curves.";
    return recs;
else
    error "Finding factors not yet implemented for genus larger than 2";
end if;

end intrinsic;
