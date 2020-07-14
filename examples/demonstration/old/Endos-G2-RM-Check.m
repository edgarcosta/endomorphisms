AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 1);

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

print "";
print "Field:";
print F;
print "";
print "Curve:";
print X;
print "";
print "Tangent representation:";
print T;
print "";
print "Minimal polynomial:";
print MinimalPolynomial(T);

print "";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, X, P0, T : LowerBound := 1);

R<x,y> := Parent(fs[1]);
Y := X;
fs := [ X`KU ! f : f in fs ];
ceqs := Y`cantor_eqs;

Y := BaseExtend(Y, X`KU);
R<x> := PolynomialRing(BaseRing(Y));
J := Jacobian(Y);

a := x^2 + fs[1]*x + fs[2];
b := fs[3]*x + fs[4];
div1 := J ! [a, b];

f := -x^5;
h := x^3 + x + 1;

Q0 := Y ! [0, 0, 1];
Q0m := Y ! [0, -1, 1];
div0 := Q0 - Q0m;
print "";
print "Improved Cantor representation:";
print div1 - div0;
