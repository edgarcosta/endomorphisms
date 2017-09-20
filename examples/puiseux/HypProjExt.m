AttachSpec("../../spec");
SetVerbose("EndoCheck", 3);

F := Rationals();
R<t> := PolynomialRing(F);
F<s> := NumberField(t^2 - 2);
F := SplittingField((t^2 - 2)*(t^2 - 10));
s2 := Roots(t^2 - 2, F)[1][1];
s10 := Roots(t^2 - 10, F)[1][1];

R<x> := PolynomialRing(F);
fX := 10*x^6 - 24*x^5 + 34*x^4 - 29*x^3 + 17*x^2 - 6*x + 1;
hX := R ! 1;
X := HyperellipticCurve(fX, hX);
P0 := X ! [1, s10, 0];

fY := x^3 + 1/3*(24*s2 - 52)*x + 1/27*(3312*s2 - 4448);
Y := HyperellipticCurve(fY);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [
[1/2, s2/2]
]);

print "Field:";
print F;
print "Curve X:";
print X;
print "Curve Y:";
print Y;
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print fs;

print "Verification:", CorrespondenceVerifyG1(X, P0, Y, Q0, M, fs);

exit;
