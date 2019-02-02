SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 1);

prec := 1000;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := (-7 + x)*(-5 + x)*(4 + x)*(8 + x)*(17 + x)*(19 + x)*(20 + x); h := R ! 0;
//f := 8 * x - 15 * x^2 + 6 * x^3 - 63 * x^4 + 48 * x^5 + 15 * x^6 + x^7;

//f := 2*x^10 + 6*x^9 + 6*x^8 + 12*x^7 + 7*x^6 + 7*x^4 - 12*x^3 + 6*x^2 - 6*x + 2; h := R ! 0;
//f := 10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10;

X := HyperellipticCurve(f, h);
X := ReducedMinimalWeierstrassModel(X);
print "";
print "Curve:";
print X;

P := PeriodMatrix(X);
EndoRep := GeometricEndomorphismRepresentation(P, F);

Y := HyperellipticCurve(x^3 - x^2 - 14916*x - 205884);
tup, d0 := MorphismOfSmallDegree(X, Y);
print "";
print "Degree:";
print d0;
print "";
print "Map:";
print tup;
