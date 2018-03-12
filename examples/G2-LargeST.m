/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

AttachSpec("../endomorphisms/magma/spec");
SetVerbose("EndoFind", 1);

F := QQ;
R<x> := PolynomialRing(F);

// TODO: Large ST (integrate with Sage, but also try without)
f := x^6 - 5*x^4 + 10*x^3 - 5*x^2 + 2*x - 1; h := R ! 0;
X := HyperellipticCurve(f, h);

prec := 1000;
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
/*
fs := RelativeMinimalPolynomials(GeoEndoRepPartial, F);
K := RelativeSplittingFieldExtra(fs);
print "Endomorphism field:";
print K;
*/

K<r> := NumberFieldExtra(x^48 - 8*x^46 + 38*x^44 - 258*x^42 + 1481*x^40 - 5354*x^38 + 13470*x^36 - 27892*x^34 + 52210*x^32 - 85828*x^30 + 120366*x^28 - 147390*x^26 + 158497*x^24 - 147390*x^22 + 120366*x^20 - 85828*x^18 + 52210*x^16 - 27892*x^14 + 13470*x^12 - 5354*x^10 + 1481*x^8 - 258*x^6 + 38*x^4 - 8*x^2 + 1);
GeoEndoRep := GeometricEndomorphismRepresentationRecognition(GeoEndoRepPartial, K);

exit;

print "";
print "Endomorphism representations:";
print GeoEndoRep;

print "";
print "More tests:";
for tup in GeoEndoRep do
    A := tup[1];
    print "Endomorphism:";
    print A;
    print "Rosati involution:";
    print RosatiInvolution(GeoEndoRep, A);
    print "Degree estimate:";
    print DegreeEstimate(GeoEndoRep, A);
    print "Minimal polynomial:";
    print MinimalPolynomial(A);
    print "Base ring:";
    print BaseRing(A);
end for;

/*
Endomorphism:
[1/1639416*(-5601*K.1^5 + 7589*K.1^4 + 9618*K.1^3 - 1491087*K.1^2 - 7236999*K.1 +
    8622692) 1/3278832*(5601*K.1^5 - 7589*K.1^4 - 9618*K.1^3 + 1491087*K.1^2 +
    8876415*K.1 - 8622692) 0]
[0 K.1 0]
[0 0 1/819708*(-287*K.1^5 - 5892*K.1^4 + 20616*K.1^3 - 143433*K.1^2 - 1716338*K.1
    - 3134772)]
Rosati involution:
[1/1639416*(4329*K.1^5 + 4379*K.1^4 + 7494*K.1^3 + 1198119*K.1^2 + 7501575*K.1 +
    5997476) 1/3278832*(-6009*K.1^5 + 1117*K.1^4 - 46758*K.1^3 - 1677855*K.1^2 -
    8791551*K.1 - 5474084) 0]
[0 1/68309*(-70*K.1^5 + 229*K.1^4 - 1636*K.1^3 - 19989*K.1^2 - 53749*K.1 + 21808)
    0]
[0 0 1/819708*(1763*K.1^5 - 2840*K.1^4 - 9540*K.1^3 + 529785*K.1^2 + 1409330*K.1 -
    3617300)]
Degree estimate:
168
Minimal polynomial:
$.1^3 + 1/1639416*(6175*K.1^5 + 4195*K.1^4 - 50850*K.1^3 + 1777953*K.1^2 +
    9030259*K.1 - 2353148)*$.1^2 + 1/546472*(2585*K.1^5 + 6181*K.1^4 - 56686*K.1^3
    + 528359*K.1^2 + 4995349*K.1 - 15599116)*$.1 + 1/409854*(6633*K.1^5 +
    745*K.1^4 + 84762*K.1^3 + 1481319*K.1^2 + 11530737*K.1 + 9932500)
Base ring:
Number Field with defining polynomial x^6 - x^5 + 2*x^4 + 255*x^3 + 1291*x^2 -
    784*x + 2192 over the Rational Field
*/

print "Verifying saturatedness:";
print VerifySaturated(GeoEndoRep, P);
print "Endomorphisms over the algebraic closure:";
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
