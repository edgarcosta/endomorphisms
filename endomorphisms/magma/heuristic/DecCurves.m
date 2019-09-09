/***
 *  Find curves in Jacobian
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


function ReconstructionFromComponentG1(P, Q, mor : AllPPs := false, ProjToIdem := true, ProjToPP := true)

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
E, h := ReconstructCurve(Q, K);
E`period_matrix := Q;

Anew := ConjugateMatrix(h, A);
if not AllPPs then
    return [* E, [* Anew, Rnew *] *];
end if;
return [ [* E, [* Anew, Rnew *] *] ];

end function;


function ReconstructionFromComponentG2(P, Q, mor : AllPPs := false, ProjToPP := true, ProjToIdem := true)

gP := #Rows(P);
A, R := Explode(mor); K := BaseRing(A);
EQ := InducedPolarization(StandardSymplecticMatrix(gP), R : ProjToIdem := ProjToIdem);
E0, _ := FrobeniusFormAlternatingAlt(EQ);

vprint CurveRec: "";
vprint CurveRec: "Frobenius form of induced polarization:";
vprint CurveRec: E0;

Ts := IsogenousPPLattices(EQ);
if not AllPPs then
    N := 1;
else
    N := #Ts;
end if;

facs := [ ];
for T in Ts[1..N] do
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
    Y`period_matrix := Qnew;

    vprint CurveRec: "";
    vprint CurveRec: "Reconstructed curve found!";
    vprint CurveRec: Y;

    Anew := ConjugateMatrix(h, A);
    Append(~facs, [* Y, [* Anew, Rnew *] *]);
end for;

if not AllPPs then
    return facs[1];
end if;
return facs;

end function;


intrinsic ReconstructionFromComponent(P::., Q::., mor::. : AllPPs := false, ProjToPP := true, ProjToIdem := true) -> .
{Given a factor (Q, mor) of the Jacobian, finds corresponding curves along with morphisms to them.}

gQ := #Rows(Q);
if gQ eq 1 then
    recs := ReconstructionFromComponentG1(P, Q, mor : AllPPs := AllPPs, ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
    return recs;
elif gQ eq 2 then
    recs := ReconstructionFromComponentG2(P, Q, mor : AllPPs := AllPPs, ProjToIdem := ProjToIdem, ProjToPP := ProjToPP);
    return recs;
else
    error "Finding factors not yet implemented for genus larger than 2";
end if;

end intrinsic;
