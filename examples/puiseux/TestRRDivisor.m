AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 4);

F := Rationals();
R<x> := PolynomialRing(F);
f := -x^5;
h := x^3 + x + 1;
X := HyperellipticCurve(f, h);
P0 := X ! [0, 0, 1];
//P0 := X ! [1, -1, 0];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

T := Matrix(F, [
[ -1, -1],
[ -1,  0]
]);

print "Curve:";
print X;
print "Calculating divisor:";
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, T : LowerBound := 1);
print D;

exit;
