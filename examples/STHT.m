/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

AttachSpec("../endomorphisms/magma/spec");
SetVerbose("EndoFind", 0);

F := QQ;
R<x> := PolynomialRing(F);
f1 := x^6 + 10*x^3 - 2; h1 := 0;
f2 := x^6 + 6*x^5 - 30*x^4 - 40*x^3 + 60*x^2 + 24*x - 8; h2 := 0;

X1 := HyperellipticCurve(f1, h1);
X2 := HyperellipticCurve(f2, h2);

prec := 300;
CCSmall := ComplexField(5);

/*
eqsCC := EmbedCurveEquations(X1, prec);
eqsF := DefiningEquations(X1);
P1 := PeriodMatrix(eqsCC, eqsF : MolinNeurohr := true);

eqsCC := EmbedCurveEquations(X2, prec);
eqsF := DefiningEquations(X2);
P2 := PeriodMatrix(eqsCC, eqsF : MolinNeurohr := true);

GeoEndoRepPartial := GeometricEndomorphismRepresentationPartial(P1);
fs := RelativeMinimalPolynomials(GeoEndoRepPartial, F);
K := RelativeSplittingFieldExtra(fs);
GeoEndoRep1 := GeometricEndomorphismRepresentationRecognition(GeoEndoRepPartial, K);

GeoEndoRepPartial := GeometricEndomorphismRepresentationPartial(P2);
fs := RelativeMinimalPolynomials(GeoEndoRepPartial, F);
K := RelativeSplittingFieldExtra(fs);
GeoEndoRep2 := GeometricEndomorphismRepresentationRecognition(GeoEndoRepPartial, K);

lat1, sthash1 := EndomorphismLattice(GeoEndoRep1, F);
lat2, sthash2 := EndomorphismLattice(GeoEndoRep2, F);

sthash1 := CanonizeSatoTateHashes(sthash1);
sthash2 := CanonizeSatoTateHashes(sthash2);
*/

R<x1,x2> := PolynomialRing(F, 2);
x3 := 1;
f := x1^3*x2 + x1^3*x3 + x1^2*x2^2 + 3*x1^2*x2*x3 + x1^2*x3^2 - 4*x1*x2^3 - 3*x1*x2^2*x3 - 3*x1*x2*x3^2 - 4*x1*x3^3 + 2*x2^4 + 3*x2^2*x3^2 + 2*x3^4;
X := PlaneCurve(f);

eqsCC := EmbedCurveEquations(X, prec);
eqsF := DefiningEquations(X);
P := PeriodMatrix(eqsCC, eqsF : MolinNeurohr := true);

GeoEndoRepPartial := GeometricEndomorphismRepresentationPartial(P);
fs := RelativeMinimalPolynomials(GeoEndoRepPartial, F);
K := RelativeSplittingFieldExtra(fs);
GeoEndoRep := GeometricEndomorphismRepresentationRecognition(GeoEndoRepPartial, K);

lat, sthash := EndomorphismLattice(GeoEndoRep, F);
sthash := CanonizeSatoTateHashes(sthash);
print sthash;

