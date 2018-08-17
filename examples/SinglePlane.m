/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

AttachSpec("../endomorphisms/magma/spec");
SetVerbose("EndoFind", 0);

F := QQ;
R<x,y> := PolynomialRing(F, 2);
z := 1;

f := x^3*z + x^2*y^2 - 4*x^2*z^2 + 6*x*z^3 + 7*y^4 + y^2*z^2 - 3*z^4;
f := 3*x^3*z + x^2*y^2 + 3*x^2*y*z - 2*x^2*z^2 + 2*x*y^3 - 3*x*y^2*z - 5*x*y*z^2 - 4*x*z^3 + y^4 - 4*y^3*z + 3*y^2*z^2 + 2*y*z^3 + 3*z^4;
f := 2*x^3*z + 2*x^2*y*z + x*y^3 - 2*x*y^2*z + 4*x*y*z^2 - 5*x*z^3 + 2*y^4 + 2*y^3*z - 2*y^2*z^2 + 2*y*z^3 - 3*z^4;
f := x^4 + 3*x^3*y + 6*x^3*z + 2*x^2*y^2 + x^2*y*z + 4*x^2*z^2 - 6*x*y^3 + 6*x*y^2*z + 6*x*y*z^2 + 2*y^4 - 4*y^3*z - y^2*z^2 - 2*z^4;
f := x^4 + 3*x^2*y^2 + 3*x^2*y*z + 3*x^2*z^2 + y^4 + 2*y^3*z + 4*y^2*z^2 + 3*y*z^3 + 2*z^4;

X := PlaneCurve(f);

prec := 200;
CCSmall := ComplexField(5);

print "Curve:";
print X;

eqsCC := EmbedCurveEquations(X, prec);
eqsF := DefiningEquations(X);
P := PeriodMatrix(eqsCC, eqsF : MolinNeurohr := true);

print "";
print "Period matrix:";
print ChangeRing(P, CCSmall);

GeoEndoRepPartial := GeometricEndomorphismRepresentationPartial(P);
fs := RelativeMinimalPolynomials(GeoEndoRepPartial, F);
K := RelativeSplittingFieldExtra(fs);
print "Endomorphism field:";
print K;
GeoEndoRep := GeometricEndomorphismRepresentationRecognition(GeoEndoRepPartial, K);

print "";
print "Endomorphism representations:";
print GeoEndoRep;

lat, sthash := EndomorphismLattice(GeoEndoRep, F);
print "";
print "Endomorphism lattice:";
print lat;

print "";
print "Sato-Tate hash:";
print sthash;
print CanonizeSatoTateHashes(sthash);
