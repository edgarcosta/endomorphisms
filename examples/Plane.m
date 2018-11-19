/*
  An example in Magma.
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

//SetVerbose("EndoFind", 1);

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
R<x,y> := PolynomialRing(F, 2);
z := 1;

f := y^3*z - x^4 - z^4;
f := 1 + 7*x*y + 21*x^2*y^2 + 35*x^3*y^3 + 28*x^4*y^4 + 2*x^7 + 2*y^7;

X := PlaneCurve(f);
print "Curve:";
print X;

P := PeriodMatrix(X);
print "";
print "Period matrix:";
print ChangeRing(P, CCSmall);

GeoEndoRep := GeometricEndomorphismRepresentation(P, F);
L<s> := BaseRing(GeoEndoRep[1][1]);

print "";
print "Endomorphism representations:";
print GeoEndoRep;

lat, sthash := EndomorphismLattice(GeoEndoRep);
print "";
print "Endomorphism lattice:";
print lat;

print "";
print "Sato-Tate hash:";
print sthash;
