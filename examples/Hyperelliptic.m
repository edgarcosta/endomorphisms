/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

AttachSpec("../endomorphisms/magma/spec");
SetVerbose("EndoFind", 2);

prec := 100;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

// Big Sato-Tate group, this calculation takes about half an hour:
f := x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1; h := R ! 0;
// CM:
f := x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8; h := 0;

R<t> := PolynomialRing(Rationals());
F<r> := NumberFieldExtra(t^2 - 5);
R<x> := PolynomialRing(F);
f := x^5 + x + 1; h := R ! 0;
f := x^5 + r*x^3 + x; h := R ! 0;

R<t> := PolynomialRing(Rationals());
F<r> := NumberFieldExtra(t^2 - t + 1);
R<x> := PolynomialRing(F);
f := x^6 + r; h := R ! 0;
f := R ! [ -30*r + 42, -156*r + 312, -66*r + 186, -1456*r + 1040, -90*r + 126, 156*r - 312, -22*r + 62 ]; h := R ! 0;

X := HyperellipticCurve(f, h);
print "Curve:";
print X;

P := PeriodMatrix(X);
print "";
print "Period matrix:";
print ChangeRing(P, CCSmall);

GeoEndoRep := GeometricEndomorphismRepresentation(P, F);
L<s> := BaseRing(GeoEndoRep[1][1]);

/* Entries can be made relative by using RelativeField if so desired */
print "";
print "Geometric endomorphism representations:";
print GeoEndoRep;

time lat, sthash := EndomorphismLattice(GeoEndoRep);
print "";
print "Endomorphism lattice:";
print lat;

print "";
print "Sato-Tate hash:";
print sthash;
