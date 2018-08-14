/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

AttachSpec("../endomorphisms/magma/spec");
SetVerbose("EndoFind", 0);

F := QQ;
R<x> := PolynomialRing(F);
f := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1; h := 0;
f := x^6 + x^2 + 1; h := 0;
f := 15*x^5 + 50*x^4 + 55*x^3 + 22*x^2 + 3*x; h := x;
f := -x^5; h := x^3 + x + 1;
f := x^6 + x^2 + 1; h := R ! 0;
f := 15*x^5 + 50*x^4 + 55*x^3 + 22*x^2 + 3*x; h := x;
f := x^4 + x^3 + 2*x^2 + x + 1; h := x^3 + x^2 + x + 1;
f := x^6 - 8*x^4 - 8*x^3 + 8*x^2 + 12*x - 8; h := 0;
f := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1; h := R ! 0;
f := x^5 + x^4 + 2*x^3 + x^2 + x; h := x^2 + x;
f := 3*x^3 - 2*x^2 + 6*x + 2; h := x^3 + x;
f := x^7 - 14*x^6 + 210*x^5 - 658*x^4 + 245*x^3 + 588*x^2 + 637*x - 686; h := 0;
f := 21*x^7 + 37506*x^5 + 933261*x^3 + 5841759*x; h := 0;
f := x^7 + 6*x^5 + 9*x^3 + x; h := 0;
f := 16*x^7 + 357*x^5 - 819*x^3 + 448*x; h := 0;
f := -4*x^8 + 105*x^6 - 945*x^4 + 2100*x^2 - 5895*x + 420; h := x^4;
f := 2*x^10 + 6*x^9 + 6*x^8 + 12*x^7 + 7*x^6 + 7*x^4 - 12*x^3 + 6*x^2 - 6*x + 2; h := 0;
f := x^4 + x^2; h := x^3 + 1;
f := 10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10; h := 0;
f := 11*x^6 + 11*x^3 - 4; h := 0;
f := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1; h := 0;

/*
R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 - t + 1);
R<x> := PolynomialRing(F);
f := R ! [ -30*r + 42, -156*r + 312, -66*r + 186, -1456*r + 1040, -90*r + 126, 156*r - 312, -22*r + 62 ]; h := R ! 0;
f := x^6 + r; h := R ! 0;

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 - 5);
R<x> := PolynomialRing(F);
f := x^5 + r*x^3 + x; h := R ! 0;
*/

X := HyperellipticCurve(f, h);

R<x1,x2> := PolynomialRing(F, 2);
x3 := 1;
f := x1^4 - x1^3*x2 + 2*x1^3*x3 + 2*x1^2*x2*x3 + 2*x1^2*x3^2 - 2*x1*x2^2*x3 +
4*x1*x2*x3^2 - x2^3*x3 + 3*x2^2*x3^2 + 2*x2*x3^3 + x3^4;
f := 3278898472*x1^4 + 35774613556*x1^3*x2 - 172165788624*x1^3*x3 -
42633841878*x1^2*x2^2 + 224611458828*x1^2*x2*x3 + 362086824567*x1^2*x3^2 +
6739276447*x1*x2^3 + 195387780024*x1*x2^2*x3 + 1153791743988*x1*x2*x3^2 -
3461357269578*x1*x3^3 - 18110161476*x2^4 - 549025255626*x2^3*x3 -
482663555556*x2^2*x3^2 + 15534718882176*x2*x3^3 - 61875497274721*x3^4;

// TODO: This takes very long (problem for MN, along with spec bugging out when gc should be gj)
f := x1^3*x2 + x1^3*x3 + x1^2*x2^2 + 3*x1^2*x2*x3 + x1^2*x3^2 - 4*x1*x2^3 -
3*x1*x2^2*x3 - 3*x1*x2*x3^2 - 4*x1*x3^3 + 2*x2^4 + 3*x2^2*x3^2 + 2*x3^4;
// TODO: This gives an error (problem for MN)
//f := x1^4 - x1^3*x3 + 2*x1^3*x2 + 2*x1^2*x3*x2 + 2*x1^2*x2^2 - 2*x1*x3^2*x2 +
//4*x1*x3*x2^2 - x3^3*x2 + 3*x3^2*x2^2 + 2*x3*x2^3 + x2^4;
// TODO: This gives an error (problem for MN)
f := x1^4 - x1^3*x3 + 2*x1^3*x2 + 2*x1^2*x3*x2 + 2*x1^2*x2^2 - 2*x1*x3^2*x2 +
4*x1*x3*x2^2 - x3^3*x2 + 3*x3^2*x2^2 + 2*x3*x2^3 + x2^4;
f := -x1^3*x2 + x1^2*x2^2 + 5*x1^2*x2*x3 - x1^2*x3^2 + 4*x1*x2^3 + 5*x1*x2^2*x3
+ x1*x2*x3^2 + 4*x1*x3^3 + 2*x2^4 + x2^2*x3^2 + 3*x3^4;
f := x1^3*x3 - x1*x3^3 + x2^4;
f := x1^3*x3 + x1*x3^3 + x2^4;

X := PlaneCurve(f);

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

rep_part := GeometricEndomorphismRepresentationPartial(P);
fs := RelativeMinimalPolynomials(rep_part, F);
K := RelativeSplittingFieldExtra(fs);
print "Endomorphism field:";
print K;

rep := GeometricEndomorphismRepresentationRecognition(rep_part, K);
galrep := GaloisRepresentationOfConjugacyClasses(rep);
G, trs := TracesOfConjugacyClasses(rep);
print G, trs;

hash1 := CanonizeTraces(G, trs);
print hash1;

hash2 := CanonizeRepresentation(rep);
print hash2;
