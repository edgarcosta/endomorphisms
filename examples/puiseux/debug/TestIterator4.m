import "../../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;
import "../../../endomorphisms/magma/puiseux/Branches.m": DevelopPoint, InitializeLift, CreateLiftIterator;

X := HyperellipticCurve(x^5 + x + 1);
P0 := X ! [0,1];
M := Matrix(Rationals(), [[2,0],[0,2]]);
M := -M;

/*
time InitializeCurve(X, P0);
print M;
print "Initializing...";
Iterator := InitializedIterator(X, X, M, 2*X`g + 2);
print "done.";
print Iterator;

for i:=1 to 5 do
    IteratorNew := IterateIterator(Iterator);
    QsNew := IteratorNew[2]; Qs := Iterator[2];
    Iterator := IteratorNew;
    print QsNew[1][2] - Qs[1][2];
end for;
*/

time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 1);
print "";
print "Divisor:";
print D;

time test, fs := CantorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print "";
print "Cantor representation:";
print fs;

Y := X;
fs := [ X`KU ! f : f in fs ];
ceqs := Y`cantor_eqs;

print "";
print "Check 0:";
print [ Evaluate(ceq, fs) : ceq in ceqs ];

KU<u,v> := X`KU;
Y := BaseExtend(Y, X`KU);
R<x> := PolynomialRing(BaseRing(Y));
J := Jacobian(Y);

a := x^2 + fs[1]*x + fs[2];
b := fs[3]*x + fs[4];
div1 := J ! [a, b];

Q0 := Y ! [0, 1, 1];
Q0m := Y ! [0, -1, 1];
div0 := Q0 - Q0m;
rep := div1 - div0;
print "";
print "Improved Cantor representation:";
print rep;

a := x - u;
b := v/u*x;
div1 := J ! [a, b];
print "";
print "Even better:";
print 2*div1;
