//SetVerbose("EndoFind", 1);
//SetVerbose("CurveRec", 1);

prec := 600;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := (-7 + x)*(-5 + x)*(4 + x)*(8 + x)*(17 + x)*(19 + x)*(20 + x); h := R ! 0;
f := 2*x^10 + 6*x^9 + 6*x^8 + 12*x^7 + 7*x^6 + 7*x^4 - 12*x^3 + 6*x^2 - 6*x + 2; h := R ! 0;

X := HyperellipticCurve(f, h);
X := ReducedMinimalWeierstrassModel(X);
print "Curve:";
print X;

P := PeriodMatrix(X);
EndoRep, L := GeometricEndomorphismRepresentation(P, F);

idems := IsotypicalIdempotents(P, EndoRep);
print idems;

//comps := IsotypicalComponentsWithInclusions(P, EndoRep);
comps := IsotypicalComponentsWithProjections(P, EndoRep);
comp := comps[1];
print comp[3];

Q, h, incdata := Explode(comp);
print SplittingIdempotents(Q, h, incdata); 

roots := RootsOfIsotypicalComponentWithProjections(Q, h, incdata); 
root := roots[1];
print root[3];

//comps := SplitComponentsWithInjections(P, EndoRep);
comps := SplitComponentsWithProjections(P, EndoRep);
print comps;
