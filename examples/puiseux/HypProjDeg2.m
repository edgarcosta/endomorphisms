AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 0);

F := Rationals();
R<x> := PolynomialRing(F);
fX := x^6 + x^2 + 1;
hX := 0;
X := HyperellipticCurve(fX, hX);
P0 := X ! [1, 1, 0];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

fY := x^3 - 16/3*x + 1856/27;
fY := x^3 + 16*x + 64;
Y := HyperellipticCurve(fY);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [
[1, 0]
]);
M := Matrix(F, [
[0, 1]
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
time test, fs := CantorFromMatrixAmbientSplit(X, P0, Y, Q0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print fs;

print "Verification:", CorrespondenceVerifyG1(X, Y, M, fs);

exit;
