SetSeed(1);
SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 2000;
R<x> := PolynomialRing(RationalsExtra(prec));

f := 10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10;
f := x^5 - x;
f := x^8 + x^6 + 5*x^4 - 3*x^2 + 17;
f := x^9 + 2*x;
X := HyperellipticCurve(f);

time dec := HeuristicDecomposition(X : SmallestField := false);

/*
time HeuristicDecomposition(X);

time P := PeriodMatrix(X);
time EndoRep := HeuristicEndomorphismRepresentation(X);
time GeoEndoRep := GeometricEndomorphismRepresentation(X);

time facs := IsotypicalFactorsOverBase(P, EndoRep);
time K := IsotypicalField(X);

time idems := SplittingIdempotentsOverBase(P, EndoRep);
time idems := SplittingIdempotentsOverClosure(P, GeoEndoRep);

time facs := SplittingFactorsOverBase(P, EndoRep);
time facs := SplittingFactorsOverClosure(P, GeoEndoRep);

time decbasedesc := DecompositionOverBase(P, EndoRep);
time decgeodesc := DecompositionOverClosure(P, GeoEndoRep);

time HeuristicDecomposition(X);
*/
