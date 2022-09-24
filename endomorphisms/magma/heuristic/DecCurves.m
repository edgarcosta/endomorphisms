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


function ReconstructionFromFactorG1(P, Q, mor : Base := true)

A, R := Explode(mor); K := BaseRing(A);
if not IsBigPeriodMatrix(Q) then
    Q := Matrix([ [ Q[1,2], Q[1,1] ] ]);
    T := Matrix(Rationals(), [[0, 1], [1, 0]]);
end if;
assert IsBigPeriodMatrix(Q);
E, h, test := ReconstructCurve(Q, K : Base := Base);
if not test then
    error "Failed to algebraize one of g4 and g6. Try increasing the precision.";
end if;
vprint CurveRec: "";
vprint CurveRec: "Reconstructed curve found!";
vprint CurveRec: E;
E`period_matrix := Q;

Anew := ConjugateMatrix(h, A); Rnew := R;
return [* E, [* Anew, Rnew *] *];

end function;


function ReconstructionFromFactorG2(P, Q, mor : AllPPs := false, Base := true)

gP := #Rows(P);
A, R := Explode(mor); K := BaseRing(A);
EQ := InducedPolarization(StandardSymplecticMatrix(gP), R);
assert IsPolarization(EQ, Q);

E0, _ := FrobeniusFormAlternatingAlt(EQ);
vprint CurveRec: "";
vprint CurveRec: "Frobenius form of induced polarization:";
vprint CurveRec: E0;

Ts := IsogenousPPLattices(EQ);
facs := [ ];
for T in Ts do
    Qnew := Q*ChangeRing(Transpose(T), BaseRing(Q)); Rnew := Transpose(T)^(-1)*R;
    assert IsBigPeriodMatrix(Qnew);

    Y, h, test := ReconstructCurve(Qnew, K : Base := Base);
    if test then
        vprint CurveRec: "";
        vprint CurveRec: "Reconstructed curve found!";
        vprint CurveRec: Y;

        Y`period_matrix := Qnew;
        Anew := ConjugateMatrix(h, A);
        if not AllPPs then
            return [* Y, [* Anew, Rnew *] *];
        end if;
        Append(~facs, [* Y, [* Anew, Rnew *] *]);
    end if;
end for;
if not AllPPs then return [* *]; else return [ ]; end if;

end function;


intrinsic ReconstructionFromFactor(P::., Q::., mor::. : AllPPs := false, Base := true) -> .
{Given a factor (Q, mor) of the Jacobian, finds a corresponding curve along with a morphism to it.}

gQ := #Rows(Q);
if gQ eq 1 then
    return ReconstructionFromFactorG1(P, Q, mor : Base := Base);
elif gQ eq 2 then
    return ReconstructionFromFactorG2(P, Q, mor : AllPPs := AllPPs, Base := Base);
else
    error "Finding factors not yet implemented for genus larger than 2";
end if;

end intrinsic;
