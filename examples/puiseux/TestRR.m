AttachSpec("../../endomorphisms/magma/spec");
import "../../endomorphisms/magma/puiseux/Branches.m": PuiseuxRamificationIndex, InitializeLift, CreateLiftIterator;
import "../../endomorphisms/magma/puiseux/FractionalCRT.m": ReduceRationalFunctionSplit;
import "../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve, ChangeTangentAction;
import "../../endomorphisms/magma/puiseux/RiemannRoch.m": RRGenerators, RRBasis, RREvaluate;
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
print "Riemann-Roch bases:";
for d in [1..20] do
    print d;
    print RRBasis(X, d);
    print [* ReduceRationalFunctionSplit(f, p, 1) : f in RRBasis(X, d) *];
    //print [ K ! ReduceRationalFunctionSplit(f, p, 1) : f in RRBasis(X, d) ];
end for;

exit;
