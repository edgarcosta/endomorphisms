//SetVerbose("EndoFind", 1);
//SetVerbose("CurveRec", 1);

prec := 600;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := (-7 + x)*(-5 + x)*(4 + x)*(8 + x)*(17 + x)*(19 + x)*(20 + x); h := R ! 0;
//f := 2*x^10 + 6*x^9 + 6*x^8 + 12*x^7 + 7*x^6 + 7*x^4 - 12*x^3 + 6*x^2 - 6*x + 2; h := R ! 0;

X := HyperellipticCurve(f, h);
X := ReducedMinimalWeierstrassModel(X);
print "Curve:";
print X;

P := PeriodMatrix(X);
EndoRep := GeometricEndomorphismRepresentation(P, F);

comps := SplitComponents(P, EndoRep);
Q, mor := Explode(comps[1]);

//Y, h := 
