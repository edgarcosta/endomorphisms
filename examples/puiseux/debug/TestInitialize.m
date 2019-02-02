import "../../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;
import "../../../endomorphisms/magma/puiseux/Branches.m": DevelopPoint, InitializeLift, CreateLiftIterator;

R<t> := PolynomialRing(Rationals());
prec := 1000;
F<r> := BaseNumberFieldExtra(t^2 + t + 2, prec);
D := [-5..5];

R<x> := PolynomialRing(F);
f := x^3 - x^2 - 2680*x + 66322; h := x + 1;
//f := x^3 - 3472875*x + 3063075750; h := 0;
X := HyperellipticCurve(f, h);
P0 := X ! [1, 0, 0];
time InitializeCurve(X, P0 : NonWP := true);

M := Matrix(F, 1, 1, [[r]]);
print "Initializing...";
Iterator := InitializedIterator(X, X, M, 2*X`g + 2);
print "done.";

/*
while true do
    IteratorNew := IterateIterator(Iterator);
    QsNew := IteratorNew[2]; Qs := Iterator[2];
    print QsNew[1][2] - Qs[1][2];
    Iterator := IteratorNew;
    print QsNew[1];
end while;
*/

time test, fs := CantorFromMatrixAmbientSplit(X, P0, X, P0, M);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

print CorrespondenceVerifyG1(X, X, M, fs : CheckDegree := true);
