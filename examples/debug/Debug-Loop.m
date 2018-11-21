SetVerbose("EndoFind", 2);
SetVerbose("CurveRec", 2);

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := x^8 + x^6 + 1;
X := HyperellipticCurve(f);
X := ReducedMinimalWeierstrassModel(X);

facs := HeuristicJacobianFactors(X);
print facs;
