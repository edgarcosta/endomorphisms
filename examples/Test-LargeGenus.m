SetVerbose("EndoFind", 3);
SetVerbose("CurveRec", 0);

prec := 50;

// Takes about one minute
F := RationalsExtra(prec);
R<x,y> := PolynomialRing(F, 2);
f := 1 + 7*x*y + 21*x^2*y^2 + 35*x^3*y^3 + 28*x^4*y^4 + 2*x^7 + 2*y^7;
X := PlaneCurve(f);

print "";
print "Curve:";
print X;

Lat := HeuristicEndomorphismLattice(X);
print "";
print "Heuristic endomorphism lattice:";
print Lat;
