SetVerbose("EndoFind", 1);
SetVerbose("CurveRec", 1);

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

/* All decomposing Sato--Tate groups are represented */
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![-7, 0, 0, 0, 1], R![0, 1, 0, 1]);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, 0, 1, 0, 1], R![1, 0, 0, 1]);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![1, 1, 2, 1, 1], R![1, 1, 1, 1]);

print "";
print "Curve:";
print X;

facs := HeuristicJacobianFactors(X : AllIdems := false, AllPPs := false);
print "";
print "Heuristic Jacobian factors:";
print facs;

for fac in facs do
    Y, mor := Explode(fac);
    test, fs := Correspondence(X, Y, mor);
    R<x,y> := Parent(fs[1]);
    print fs;
end for;
