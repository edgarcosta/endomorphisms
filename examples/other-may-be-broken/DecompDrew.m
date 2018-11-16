/*
  An example in Magma.
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

SetVerbose("EndoFind", 2);
SetVerbose("CurveRec", 2);

prec := 1000;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := 2*x^10 + 6*x^9 + 6*x^8 + 12*x^7 + 7*x^6 + 7*x^4 - 12*x^3 + 6*x^2 - 6*x + 2; h := R ! 0;
f := 10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10; h := R ! 0;

X := HyperellipticCurve(f, h);
X := ReducedMinimalWeierstrassModel(X);
print "Curve:";
print X;

time P := PeriodMatrix(X);
time EndoRep := GeometricEndomorphismRepresentation(P, F);
print EndomorphismStructureFromRepresentation(EndoRep);

idems := IdempotentsFromRepresentation(EndoRep);
Q := Ker0(idems[1], P, P);
E := PolarizationBasis(Q)[1];
print FrobeniusFormAlternatingAlt(E);

Us := IsogenousPPLatticesG2(E);
Ys := [* *];
for U in Us do
    Qnew := Q*ChangeRing(U^(-1), BaseRing(Q));
    //Qnew := Q*ChangeRing(Us[1]^(-1), BaseRing(Q));
    assert IsBigPeriodMatrix(Qnew);
    Y := ReconstructCurveG2(Qnew, F);
    Append(~Ys, Y);
end for;
print [* G2Invariants(Y) : Y in Ys *];
print [* ReducedMinimalWeierstrassModel(Y) : Y in Ys *];

