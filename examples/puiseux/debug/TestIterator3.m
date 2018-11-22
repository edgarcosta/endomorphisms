import "../../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;
import "../../../endomorphisms/magma/puiseux/Branches.m": DevelopPoint, InitializeLift, CreateLiftIterator;

prec := 1000;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
D := [-5..5];

f := x^8 + x^6 + 1; h := R ! 0;
X := HyperellipticCurve(f, h);
P0 := X ! [0, 1, 1];
time InitializeCurve(X, P0 : NonWP := true);

f := x^3 - 64*x + 64; h := R ! 0;
Y := HyperellipticCurve(f, h);
Q0 := Y ! [1, 1, 1];
//Q0 := Y ! [1, 0, 0];
time InitializeCurve(Y, Q0 : NonWP := false);

M := Matrix(F, 1, 3, [ Random(D) : i in [1..3] ]);
print M;
print "Initializing...";
Iterator := InitializedIterator(X, Y, M, 2*Y`g + 2);
print "done.";

/* Further iteration seems perfect now, no exceptions found yet */
while true do
    IteratorNew := IterateIterator(Iterator);
    QsNew := IteratorNew[2]; Qs := Iterator[2];
    print QsNew[1][2] - Qs[1][2];
    Iterator := IteratorNew;
end while;
