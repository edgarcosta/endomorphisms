AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

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
//time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, T : LowerBound := 6);
time test, D := DivisorFromMatrixRRSplit(X, P0, X, P0, T : LowerBound := 3);

R<v2,v1,u2,u1> := Parent(DefiningEquations(D)[1]);
print D;

exit;

print "Irreducible components:";
print IrreducibleComponents(D);

exit;
