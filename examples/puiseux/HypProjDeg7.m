AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

F := Rationals();
R<x> := PolynomialRing(F);
fX := x^6 + 3*x^5 + 2*x^4 + 7*x^3 + 11*x^2 + 14;
hX := x^2 + x;
X := HyperellipticCurve(fX, hX);
P0 := X ! [1, 1, 0];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

fY := x^3 - 122451/16*x + 8353079/32;
//fY := x^3 + 16934253/16*x - 5499502505/32;
Y := HyperellipticCurve(fY);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [
[1, 1/7]
]);
//M := Matrix(F, [
//[0, -1/7]
//]);

print "Field:";
print F;
print "Curve X:";
print X;
print "Curve Y:";
print Y;
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, Y, Q0, M : LowerBound := 56);
print fs;

print "Verification:", CorrespondenceVerifyG1(X, Y, M, fs);

exit;
