import "../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;
import "../../endomorphisms/magma/puiseux/Branches.m": DevelopPoint, InitializeLift, CreateLiftIterator;

prec := 1000;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
D := [-5..5];
repeat
    f := x^8 + &+[ Random(D)*x^i : i in [0..7] ];
    f := f - Evaluate(f, 0) + 1;
    X := HyperellipticCurveExtra(f, 0, prec);
    P0 := X ! [0, 1];
until not IsWeierstrassPlace(Place(P0));
print f;
time InitializeCurve(X, P0);

prec := 1000;
F := RationalsExtra(prec);
R<x,y,z> := PolynomialRing(F, 3);
D := [-5..5];
repeat
    f := &+[ Random(D)*mon : mon in MonomialsOfDegree(R, 4) ];
    f := f - Evaluate(f, [1, 1, 1])*z^4;
    X := PlaneCurveExtra(f, prec);
    P0 := X ! [1, 1, 1];
until not IsWeierstrassPlace(Place(P0));
print f;
time InitializeCurve(X, P0);

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

M := Matrix(Rationals(), 3, 3, [ Random(D) : i in [1..9] ]);
P, Qs, IterateLift := InitializedIterator(X, X, M, X`g + 5);
print P;
print Qs;

/* Further iteration seems perfect now, no exceptions found yet */
while true do
    Pnew, Qsnew := IterateIterator(P, Qs, IterateLift);
    print Qsnew[3][2] - Qs[3][2];
    P := Pnew; Qs := Qsnew;
end while;
