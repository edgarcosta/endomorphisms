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
f := 10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10; h := 0;
g := 4*x^6 + 8*x^5 + 11*x^4 + 7*x^3 - 7*x^2 - 23*x;

X := HyperellipticCurve(f);
Y := HyperellipticCurve(g);

time P := PeriodMatrix(X);
time Q := PeriodMatrix(Y);
time HomRepCC := GeometricHomomorphismRepresentationCC(P, Q);
print HomRepCC;
