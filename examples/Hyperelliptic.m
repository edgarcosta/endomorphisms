SetVerbose("EndoFind", 0);

prec := 300;
CCSmall := ComplexField(5);

F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

// Big Sato-Tate group, this calculation takes about 20 minutes:
f := x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1; h := R ! 0;
// CM:
f := x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8; h := 0;

R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - 5, prec);
R<x> := PolynomialRing(F);
f := x^5 + x + 1; h := R ! 0;
f := x^5 + r*x^3 + x; h := R ! 0;

R<t> := PolynomialRing(Rationals());
F<r> := BaseNumberFieldExtra(t^2 - t + 1, prec);
//F<r> := BaseNumberFieldExtra(t^2 - 5, prec);
R<x> := PolynomialRing(F);
f := R ! [ -30*r + 42, -156*r + 312, -66*r + 186, -1456*r + 1040, -90*r + 126, 156*r - 312, -22*r + 62 ]; h := R ! 0;
f := x^6 + r; h := R ! 0;

X := HyperellipticCurve(f, h);
print "Curve:";
print X;

P := PeriodMatrix(X);
print "";
print "Period matrix:";
print ChangeRing(P, CCSmall);

time GeoEndoRep := GeometricEndomorphismRepresentation(P, F);

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
