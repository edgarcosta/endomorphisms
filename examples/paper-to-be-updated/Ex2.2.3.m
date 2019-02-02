SetVerbose("EndoFind", 0);

prec := 300;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

f := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1;
X := HyperellipticCurve(f);

print "Curve:";
print X;

P := PeriodMatrix(X);
time GeoEndoRep := GeometricEndomorphismRepresentation(P, F);

/* Entries can be made relative by using RelativeField if so desired */
print "";
print "Geometric endomorphism representations:";
print GeoEndoRep;

time lat, sthash := EndomorphismLattice(GeoEndoRep);
print "";
print "Endomorphism lattice:";
print lat;
