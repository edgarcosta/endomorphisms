AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

R<x> := PolynomialRing(Rationals());
F<iF> := NumberField(x^2 + 1);
R<x> := PolynomialRing(F);
f := -x^6 + 5*x^5 - 8*x^4 + 4*x^3 - x^2 + x;
X := HyperellipticCurve(f);
P0 := X ! [0, 0, 1];
//P0 := X ! [1/2, 5/8, 1];
Y := HyperellipticCurve(f);
Q0 := Y ! [1/2, 5/8, 1];

T := Matrix(F, [
[ iF, -2*iF],
[ 2*iF, -iF]
]);

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, Y, Q0, T : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

Y := BaseExtend(Y, X`KU);
R<x> := PolynomialRing(BaseRing(Y));
J := Jacobian(Y);

fs := [ X`KU ! f : f in fs ];
a := x^2 + fs[1]*x + fs[2];
b := fs[3]*x + fs[4];
div1 := J ! [a, b];

Q0 := Y ! [1/2, 5/8, 1];
Q0m := Y ! [1/2, -5/8, 1];
div0 := Q0 - Q0m;
print "Improved Cantor representation:";
print div1 - div0;
