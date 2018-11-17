SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 1);

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, 0, 0, 0, 1, 1], R![1, 1, 0, 1]);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![-2, -8, -10, -5, 0, 1], R![0, 0, 0, 1]);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![-2], R![0, 0, 0, 1]);

//R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, -1, 0, 0, 0, 1], R![]);

print "Curve:";
print C;

P := PeriodMatrix(C);
GeoEndoRep := GeometricEndomorphismRepresentation(P, F);

comps := SplitComponents(P, GeoEndoRep);
for comp in comps do
    Q, mor := Explode(comps[1]);
    recs := ReconstructionsFromComponent(P, Q, mor);
    print recs[1];
end for;
