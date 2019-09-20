import "../../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;
import "../../../endomorphisms/magma/puiseux/Branches.m": DevelopPoint, InitializeLift, CreateLiftIterator;

prec := 1000;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
D := [-5..5];
repeat
    f := x^8 + &+[ Random(D)*x^i : i in [0..7] ];
    f := f - Evaluate(f, 0) + 1;
    X := HyperellipticCurve(f, 0);
    P0 := X ! [0, 1];
until not IsWeierstrassPlace(Place(P0));
print f;
time InitializeCurve(X, P0 : AssertNonWP := true);

prec := 1000;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
D := [-5..5];
repeat
    f := x^8 + &+[ Random(D)*x^i : i in [0..7] ];
    f := f - Evaluate(f, 0) + 1;
    X := HyperellipticCurve(f, 0);
    P0 := X ! [0, 1];
until not IsWeierstrassPlace(Place(P0));
print f;
time InitializeCurve(X, P0 : AssertNonWP := true);

prec := 100;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
F := NumberFieldExtra(x^6 + 42*x^4 + 441*x^2 + 784);
R<x> := PolynomialRing(F);
f := 16*x^7 + 357*x^5 - 819*x^3 + 448*x;
f /:= 2; h := R ! 0;
X := HyperellipticCurve(f, h);
// Try to make one the point at infinity!
P0 := X ! [1, 1, 1];
time InitializeCurve(X, P0 : AssertNonWP := true);

/* This is fast and does not give trouble so far */
f := X`DEs[1];
P := DevelopPoint(X, P0, 98);
for n in [98..100] do
    Pnew := DevelopPoint(X, P0, n);
    print Pnew[2] - P[2];
    print Evaluate(f, Pnew);
    print "---";
    P := Pnew;
end for;
