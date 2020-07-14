//SetVerbose("EndoFind", 1);
//SetVerbose("CurveRec", 2);

prec := 500;

R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - t + 1, prec);
R<x> := PolynomialRing(F);
f := R ! [ -30*r + 42, -156*r + 312, -66*r + 186, -1456*r + 1040, -90*r + 126, 156*r - 312, -22*r + 62 ];
X := HyperellipticCurve(f);

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
