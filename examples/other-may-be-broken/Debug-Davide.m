SetVerbose("EndoFind", 2);
SetVerbose("CurveRec", 2);

prec := 1000;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := (-7 + x)*(-5 + x)*(4 + x)*(8 + x)*(17 + x)*(19 + x)*(20 + x);

X := HyperellipticCurve(f);
X := ReducedMinimalWeierstrassModel(X);

time P := PeriodMatrix(X);
time GeoEndoRep := GeometricEndomorphismRepresentation(P, F);
GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep);

print "";
print "Geometric endomorphism algebra:";
print GeoEndoAlg;

GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];
comps := IsotypicalComponents(P, GeoEndoRep : ProjOrInc := "Inc");

print "";
print "Isotypical components:";
print [ [* comp[2], comp[3] *] : comp in comps ];

Q, mor, incdata := Explode(comps[2]);
A, R := Explode(mor);
E := InducedPolarization(StandardSymplecticMatrix(3), R : ProjOrInc := "Inc");

E0 := FrobeniusFormAlternatingAlt(E);
print "";
print "Frobenius form of polarization:";
print E0;

E := E^(-1);
E := -Denominator(E) * E;

E0 := FrobeniusFormAlternatingAlt(E);
print "";
print "Frobenius form of polarization:";
print E0;

Us := IsogenousPPLattices(E);
Ys := [* *];
for U in Us do
    Qnew := Q*ChangeRing(U^(-1), BaseRing(Q));
    assert IsBigPeriodMatrix(Qnew);
    Y := ReconstructCurve(Qnew, F);
    print "";
    print Y;
    Append(~Ys, Y);
end for;
