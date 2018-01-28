AttachSpec("../../endomorphisms/magma/spec");

F := QQ;
R<x> := PolynomialRing(F);
f := x^6 + x^2 + 1; h := 0;
f := x^6 + 3*x^5 + 6*x^4 + 7*x^3 + 6*x^2 + 3*x + 1; h := x^2 + x;
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

print "Verifying saturatedness:";
print VerifySaturated(GeoEndoRep, P);

print "Endomorphisms lattice:";
lat := EndomorphismLattice(GeoEndoRep, F);
print lat;

idems, L := IdempotentsFromLattice(lat);
print "Idempotents and their field of definition:";
print idems;
print L;

print ProjectionFromIdempotent(P, idems[1]);

print "";
print "Decomposition:";
facs_an := [ ];
projs_an := [ ];
for idem in idems do
    fac_an, proj_an := ProjectionFromIdempotent(P, idem);
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
