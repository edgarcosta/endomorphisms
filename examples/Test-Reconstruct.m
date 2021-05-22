SetSeed(1);
SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 2);

prec := 500;
R<x> := PolynomialRing(RationalsExtra(prec));

f := x^6 + x + 1;
X := HyperellipticCurve(f);

F := BaseRing(X); P := PeriodMatrix(X); CC := BaseRing(P);
R := Matrix([[2,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]);
RCC := ChangeRing(R, CC);

Q := P*RCC^(-1);
EP := StandardSymplecticMatrix(2);
assert IsPolarization(EP, P);
EQ := InducedPolarization(EP, R);
assert IsPolarization(EQ, Q);

Q := P*RCC;
EQ := -PolarizationBasis(Q)[1];
assert IsPolarization(EQ, Q);

Ts := IsogenousPPLattices(EQ); Ys := [* *];
for T in Ts do
    Qnew := Q*ChangeRing(Transpose(T), BaseRing(Q));
    assert IsBigPeriodMatrix(Qnew);
    Y, h, test := ReconstructCurve(Qnew, F : Base := true);
    Append(~Ys, Y);
end for;
print Ys;


