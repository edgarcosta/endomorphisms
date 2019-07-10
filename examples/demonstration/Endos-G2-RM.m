//SetVerbose("EndoFind", 1);
//SetVerbose("CurveRec", 2);

prec := 500;

F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
f := -x^5;
h := x^3 + x + 1;
X := HyperellipticCurve(f, h);

print "";
print "Curve:";
print X;

L := HeuristicEndomorphismFieldOfDefinition(X);
print "";
print "Heuristic field of definition of the endomorphisms:";
print L;

Lat := HeuristicEndomorphismLattice(X);
print "";
print "Heuristic endomorphism lattice:";
print Lat;
