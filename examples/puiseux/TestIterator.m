import "../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;
import "../../endomorphisms/magma/puiseux/Branches.m": DevelopPoint, InitializeLift, CreateLiftIterator;

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
time InitializeCurve(X, P0 : NonWP := true);

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
time InitializeCurve(X, P0 : NonWP := true);

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
time InitializeCurve(X, P0 : NonWP := true);

/*
prec := 1000;
F := RationalsExtra(prec);
R<x,y,z> := PolynomialRing(F, 3);
D := [-5..5];
repeat
    f := &+[ Random(D)*mon : mon in MonomialsOfDegree(R, 4) ];
    f := f - Evaluate(f, [1, 1, 1])*z^4;
    X := PlaneCurve(f);
    P0 := X ! [1, 1, 1];
until not IsWeierstrassPlace(Place(P0));
print f;
time InitializeCurve(X, P0 : NonWP := true);
*/

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

M := Matrix(F, 3, 3, [ Random(D) : i in [1..9] ]);
M := Matrix(F, [
[ 1/28*(F.1^3 + 21*F.1),                      0,                      0],
[                     0, 1/28*(-F.1^3 - 21*F.1),                      0],
[                     0,                      0,  1/28*(F.1^3 + 21*F.1)]
]);
M := Matrix(F, [
[1/168*(5*F.1^4 + 133*F.1^2 + 392), 0, 1/168*(-F.1^4 + 7*F.1^2 + 392)],
[0, 1/56*(-3*F.1^4 - 91*F.1^2 - 392), 0],
[1/84*(-F.1^4 + 7*F.1^2 + 392), 0, 1/42*(F.1^4 + 35*F.1^2 + 196)]
]);
Iterator := InitializedIterator(X, X, M, X`g + 7);

/* Further iteration seems perfect now, no exceptions found yet */
while true do
    IteratorNew := IterateIterator(Iterator);
    QsNew := IteratorNew[2]; Qs := Iterator[2];
    print QsNew[3][2] - Qs[3][2];
    Iterator := IteratorNew;
end while;
