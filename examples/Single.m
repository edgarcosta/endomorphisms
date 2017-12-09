/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

F := QQ;
R<x> := PolynomialRing(F);
f := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1; h := 0;
f := x^6 + x^2 + 1; h := 0;
//f := 15*x^5 + 50*x^4 + 55*x^3 + 22*x^2 + 3*x; h := x;
X := HyperellipticCurve(f, h);
prec := 100;
CCSmall := ComplexField(5);

print "Curve:";
print X;

eqsCC := EmbedCurveEquations(X, prec);
eqsF := DefiningEquations(X);
P := PeriodMatrix(eqsCC, eqsF : HaveOldenburg := false);

print "";
print "Period matrix:";
print ChangeRing(P, CCSmall);

GeoEndoRepPartial := GeometricEndomorphismRepresentationPartial(P);
fs := RelativeMinimalPolynomialsPartial(GeoEndoRepPartial, F);
K := RelativeSplittingFieldExtra(fs);
GeoEndoRep := GeometricEndomorphismRepresentationRecognition(GeoEndoRepPartial, K);

print "";
print "Endomorphism representations:";
print GeoEndoRep;

print "";
print "More tests:";
A := GeoEndoRep[2][1];
print "Rosati involution:";
print RosatiInvolution(GeoEndoRep, A);
print "Degree estimate:";
print DegreeEstimate(GeoEndoRep, A);
print "Verifying saturatedness:";
print VerifySaturated(GeoEndoRep, P);

print "Endomorphisms over QQ:";
struct := EndomorphismStructure(GeoEndoRep, K, F);
print struct;

print "Endomorphisms lattice:";
lat := EndomorphismLattice(GeoEndoRep, F);
print lat;

idems, L := IdempotentsFromLattice(lat);
print "Idempotents and their field of definition:";
print idems;
print L;

print "";
print "Decomposition:";
facs_an := [ ];
projs_an := [ ];
for idem in idems do
    fac_an, proj_an := ProjectionFromIdempotentNew(P, idem);
    Append(~facs_an, fac_an); Append(~projs_an, proj_an);
end for;
facs_alg := [ FactorReconstruct(P, facs_an[i], projs_an[i][3], projs_an[i][2], L) : i in [1..#facs_an] ];
print "Analytic factors:";
print facs_an;
print "Analytic projections:";
print projs_an;
print "Algebraic factors:";
print facs_alg;

exit;
