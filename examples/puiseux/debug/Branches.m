import "../../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;
import "../../../endomorphisms/magma/puiseux/Branches.m": DevelopPoint;

F := RationalsExtra();
R<x> := PolynomialRing(F);
PR<t> := PuiseuxSeriesRing(F, 100);

X := HyperellipticCurve(x^6 + 7*x + 1);
P0 := X ! [ 0, 1, 1 ];
//P0 := X ! [ 1, 3, 1 ];
InitializeCurve(X, P0);

print "Testing that result is correct to given precision and strictly increases...";
N := 10;
Pold := DevelopPoint(X, [t^(2/3), 1], 2^N - 4);
for prec in [2^N - 3..2^N + 3] do
    P := DevelopPoint(X, [t^(2/3), 1], prec);
    print P[2] - Pold[2];
    f := X`DEs[1];
    print Evaluate(f, P);
    Pold := P;
end for;
