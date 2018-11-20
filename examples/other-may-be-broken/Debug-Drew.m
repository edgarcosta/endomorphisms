/*
  An example in Magma.
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

SetVerbose("EndoFind", 2);
SetVerbose("CurveRec", 2);

prec := 1000;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := 2*x^10 + 6*x^9 + 6*x^8 + 12*x^7 + 7*x^6 + 7*x^4 - 12*x^3 + 6*x^2 - 6*x + 2;
f := 10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10;
g := 4*x^6 + 8*x^5 + 11*x^4 + 7*x^3 - 7*x^2 - 23*x;

X := HyperellipticCurve(f);
X := ReducedMinimalWeierstrassModel(X);
Y := HyperellipticCurve(g);
Y := ReducedMinimalWeierstrassModel(Y);

time P := PeriodMatrix(X);
time Q := PeriodMatrix(Y);
/*
print GeometricHomomorphismRepresentation(P, Q, F);
*/

time GeoEndoRep := GeometricEndomorphismRepresentation(P, F);
GeoEndoAlg, GeoEndoDesc := EndomorphismStructure(GeoEndoRep);

print "";
print "Geometric endomorphism algebra:";
print GeoEndoAlg;

GeoEndoData := [* GeoEndoRep, GeoEndoAlg, GeoEndoDesc *];
comps := IsotypicalComponents(P, GeoEndoRep);

print "";
print "Isotypical components:";
print [ [* comp[2], comp[3] *] : comp in comps ];

P, mor, incdata := Explode(comps[1]);
idems := SplittingIdempotents(P, mor, incdata);

print "";
print "Splitting idempotents:";
print idems;

Q := Ker0([* ChangeRing(idems[1][1], CC), idems[1][2] *], P, P);
E := PolarizationBasis(Q)[1];

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
