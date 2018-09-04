/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

AttachSpec("../endomorphisms/magma/spec");
SetVerbose("EndoFind", 2);

F := RationalsExtra();
R<x> := PolynomialRing(F);
// Big ST, this calculation takes about:
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

prec := 500;
CCSmall := ComplexField(5);

print "Curve:";
print X;

P := PeriodMatrix(X : prec := prec);

print "";
print "Period matrix:";
print ChangeRing(P, CCSmall);

GeoEndoRep := GeometricEndomorphismRepresentation(P, F);
L<s> := BaseRing(GeoEndoRep[1][1]);
F := BaseRing(L);

print "";
print "Geometric endomorphism representations:";
print GeoEndoRep;

R<x> := PolynomialRing(F);
K<s> := NumberField(x^2 + 7);
EndoRep := EndomorphismRepresentation(GeoEndoRep, K);
if not IsQQ(BaseRing(EndoRep[1][1])) then
    K<s> := BaseRing(EndoRep[1][1]);
end if;
print "";
print "Endomorphism representation over", K, ":";
print EndoRep;

lat, sthash := EndomorphismLattice(GeoEndoRep);
print "";
print "Endomorphism lattice:";
print lat;

print "";
print "Sato-Tate hash:";
print sthash;
