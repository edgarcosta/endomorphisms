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
print [ tup[2] : tup in EndoRep ];

idems := IsotypicalIdempotents(P, EndoRep);
print idems;

/*
comps_proj := IsotypicalComponentsWithProjections(P, EndoRep);
print comps_proj;

comp := comps_proj[1];
Q, h := Explode(comp);
print IsotypicalComponentFoD(Q, h); 

comps := RootsOfIsotypicalComponentWithProjections(Q, h); 
print comps;
*/

comps_inc := IsotypicalComponentsWithInclusions(P, EndoRep);
print comps_inc;

comp := comps_inc[1];
Q, h := Explode(comp);
print IsotypicalComponentFoD(Q, h); 

roots := RootsOfIsotypicalComponentWithInclusions(Q, h); 
root := roots[1];
print root;

print ReconstructCurveFromRoot(root);
