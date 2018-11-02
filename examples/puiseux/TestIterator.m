import "../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;
import "../../endomorphisms/magma/puiseux/Branches.m": DevelopPoint, InitializeLift, CreateLiftIterator;

prec := 1000;
F := RationalsExtra(prec);
R<x,y,z> := PolynomialRing(F, 3);
D := [-5..5];
f := &+[ Random(D)*mon : mon in MonomialsOfDegree(R, 4) ];
f := f - Evaluate(f, [0, 0, 1])*z^4;
X := PlaneCurveExtra(f, prec);
P0 := X ! [0, 0, 1];
time InitializeCurve(X, P0);

/* This is fast and does not give trouble so far */
f := X`DEs[1];
P := DevelopPoint(X, P0, 0);
for n in [1..100] do
    Pnew := DevelopPoint(X, P0, n);
    //print Pnew[2] - P[2];
    //print Evaluate(f, Pnew);
    P := Pnew;
end for;

M := Matrix(Rationals(), 3, 3, [ Random(D) : i in [1..9 ] ]);
print InitializeLift(X, X, M);

P, Qs := InitializeLift(X, X, M);
IterateLift := CreateLiftIterator(X, X, M);
while true do
    Pnew, Qsnew := IterateLift(P, Qs, 30);
    if Pnew eq P and Qsnew eq Qs then
        break;
    end if;
    P := Pnew;
    Qs := Qsnew;
end while;

P, Qs:= IterateLift(P, Qs, 60);
Pnew, Qsnew:= IterateLift(P, Qs, 60);
print Qsnew[3][2] - Qs[3][2];
