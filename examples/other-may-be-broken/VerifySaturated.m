AttachSpec("../../endomorphisms/magma/spec");

F := QQ;
R<x> := PolynomialRing(F);
f := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1; h := 0;
f := x^6 + x^2 + 1; h := 0;
f := 15*x^5 + 50*x^4 + 55*x^3 + 22*x^2 + 3*x; h := x;
X := HyperellipticCurve(f, h);

prec := 300;
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
GeoEndoRep := GeometricEndomorphismRepresentationRecognition(GeoEndoRepPartial, K);

print VerifySaturated(GeoEndoRep, P);

exit;
