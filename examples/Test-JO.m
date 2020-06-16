SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 200;

// Big Sato-Tate group, this calculation takes about 20 minutes and needs
// precision at least 200, unless setting UseQQ to true.
// This flag is hidden in Recognition.m and is especially useful when minimal
// polynomials over QQ are small, as they typically are in small genus.
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
f := x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1;
X := HyperellipticCurve(f);

print "";
print "Curve:";
print X;

time P := PeriodMatrix(X);
time desc := HeuristicEndomorphismAlgebra(X : CC := true);
time desc := HeuristicEndomorphismAlgebra(X : Geometric := true);
//time rep := HeuristicEndomorphismRepresentation(X);
//time lat := HeuristicEndomorphismLattice(X);

