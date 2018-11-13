//SetVerbose("EndoFind", 1);
//SetVerbose("CurveRec", 1);

prec := 600;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;
R<x> := PolynomialRing(F);

f := x^8 + x^2 + 1; h := R ! 0;
X := HyperellipticCurve(f, h);
P := PeriodMatrix(X);

f := x*(x^4 + x + 1); h := R ! 0;
Y := HyperellipticCurve(f, h);
Q := PeriodMatrix(Y);

print MorphismOfSmallDegree(P, Q, F);
print MorphismOfSmallDegree(X, Y);
