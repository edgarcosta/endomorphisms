AttachSpec("../../endomorphisms/magma/spec");
import "../../endomorphisms/magma/puiseux/Branches.m": DevelopPoint;
import "../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;
import "../../endomorphisms/magma/puiseux/RiemannRoch.m": RRGenerators, RRBasis, RREvaluations, ProductBasis, ProductEvaluations, GlobalGenerators, GlobalBasis, GlobalProductBasis;
SetVerbose("EndoCheck", 4);

F := Rationals();
R<x> := PolynomialRing(F);
f := -x^5;
h := x^3 + x + 1;
X := HyperellipticCurve(f, h);
P0 := X ! [0, 0, 1];

T := Matrix(F, [
[ -1, -1],
[ -1,  0]
]);

print "Curve:";
print X;
InitializeCurve(X, P0);

p := NextPrime(10^100);
FF := FiniteField(p);
K<x,y> := RationalFunctionField(FF, 2);
print "Riemann-Roch generators:";
print [ X`KA ! f : f in RRGenerators(X) ];

d := 12; prec := 20;
print "Riemann-Roch basis:";
print [ X`KA ! f : f in RRBasis(X, d) ];
P := DevelopPoint(X, P0, prec);
R<pi> := Parent(P[1]);
print "Riemann-Roch evaluations:";
print RREvaluations(X, d, P);

print "Product basis:";
print #ProductBasis(X, X, d);
print "Product evaluations:";
print #ProductEvaluations(X, X, d, P, [ P ])[1];

print "Global generators:";
print GlobalGenerators(X);
print "Global basis:";
print GlobalBasis(X, d);

print "Global product basis:";
print #GlobalProductBasis(X, X, d);

exit;
