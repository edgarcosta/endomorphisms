import "../../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve, OurAffinePatch;

SetVerbose("EndoCheck", 4);

prec := 1000;
F := RationalsExtra(prec);

R<x> := PolynomialRing(F);
f := x*(x + 1)*(x + 7);
X := HyperellipticCurve(f);
P0 := X ! [1, 4, 1];
P0 := X ! [-1, 0, 1];
P0 := X ! [1, 0, 0];

time InitializeCurve(X, P0 : NonWP := false);
