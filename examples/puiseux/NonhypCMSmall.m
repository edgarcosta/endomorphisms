AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 2);

QQ := Rationals();
R<t> := PolynomialRing(QQ);
F<r> := NumberField(t^6 + t^5 + 2*t^4 - 8*t^3 - t^2 + 5*t + 7);

/*
P2<x,y,z> := ProjectiveSpace(F, 2);
fX := x^4 - x^3*y + 2*x^3*z + 2*x^2*y*z + 2*x^2*z^2 - 2*x*y^2*z + 4*x*y*z^2 - y^3*z + 3*y^2*z^2 + 2*y*z^3 + z^4;
X := Curve(P2, fX);
P0 := X ! [1, 1, 0];

M := Matrix(F, [
[1/77*(3*r^5 + 5*r^4 + 35*r^3 + 25*r^2 + 65*r - 70), 0, 0],
[0, 1/77*(9*r^5 + 15*r^4 + 28*r^3 - 79*r^2 - 36*r + 21), 0],
[0, 0, 1/77*(16*r^5 + r^4 + 7*r^3 - 149*r^2 + 90*r + 140)]
]);
*/

P2<x,y,z> := ProjectiveSpace(F, 2);
fX := y^4 - y^3*z + 2*y^3*x + 2*y^2*z*x + 2*y^2*x^2 - 2*y*z^2*x + 4*y*z*x^2 - z^3*x + 3*z^2*x^2 + 2*z*x^3 + x^4;
X := Curve(P2, fX);
P0 := X ! [0, 1, 1];

M := Matrix(F, [
[1/77*(16*r^5 + r^4 + 7*r^3 - 149*r^2 + 90*r + 140), 0, 0],
[0, 1/77*(3*r^5 + 5*r^4 + 35*r^3 + 25*r^2 + 65*r - 70), 0],
[0, 0, 1/77*(9*r^5 + 15*r^4 + 28*r^3 - 79*r^2 - 36*r + 21)]
]);

print "";
print "Field:";
print F;
print "Curve X:";
print X;
print "Point P0:";
print P0;
print "Tangent representation:";
print M;

print "Calculating divisor...";
time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 9, Margin := 2^8);
//time test, D := DivisorFromMatrixRRSplit(X, P0, X, P0, M : LowerBound := 1, Margin := 2^8);
print D;

exit;

print "Irreducible components:";
print IrreducibleComponents(D);

exit;
