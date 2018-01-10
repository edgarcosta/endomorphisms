AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

QQ := Rationals();
R<t> := PolynomialRing(QQ);
F<s> := NumberField(t^6 + t^5 + 2*t^4 - 8*t^3 - t^2 + 5*t + 7);

P2<x,y,z> := ProjectiveSpace(F, 2);
fX := x^4 - x^3*y + 2*x^3*z + 2*x^2*y*z + 2*x^2*z^2 - 2*x*y^2*z + 4*x*y*z^2 - y^3*z + 3*y^2*z^2 + 2*y*z^3 + z^4;
X := Curve(P2, fX);
P0 := X ! [1, 1, 0];

M := Matrix(F, [
[1/77*(3*s^5 + 5*s^4 + 35*s^3 + 25*s^2 + 65*s - 70), 0, 0],
[0, 1/77*(9*s^5 + 15*s^4 + 28*s^3 - 79*s^2 - 36*s + 21), 0],
[0, 0, 1/77*(16*s^5 + s^4 + 7*s^3 - 149*s^2 + 90*s + 140)]
]);

P2<x,y,z> := ProjectiveSpace(F, 2);
fX := x^4 - x^3*z + 2*x^3*y + 2*x^2*z*y + 2*x^2*y^2 - 2*x*z^2*y + 4*x*z*y^2 - z^3*y + 3*z^2*y^2 + 2*z*y^3 + y^4;
X := Curve(P2, fX);
P0 := X ! [1, 0, 1];

M := Matrix(F, [
[1/77*(3*s^5 + 5*s^4 + 35*s^3 + 25*s^2 + 65*s - 70), 0, 0],
[0, 1/77*(16*s^5 + s^4 + 7*s^3 - 149*s^2 + 90*s + 140), 0],
[0, 0, 1/77*(9*s^5 + 15*s^4 + 28*s^3 - 79*s^2 - 36*s + 21)]
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
time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 1);
print D;

exit;
