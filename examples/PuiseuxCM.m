AttachSpec("../spec");
SetVerbose("EndoCheck", 1);
SetMemoryLimit(32*10^9);

R<t> := PolynomialRing(Rationals());
f := x^6 + 13*x^4 + 50*x^2 + 49;
F<r> := NumberField(f);
R<x> := PolynomialRing(F);
f := 21*x^7 + 37506*x^5 + 933261*x^3 + 5841759*x;
X := HyperellipticCurve(f);
P0 := X ! [3, -7203];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[1/7*(r^5 + 6*r^3 + r), 0, 0],
[0, 1/7*(-r^5 - 6*r^3 - r), 0],
[0, 0, 1/7*(r^5 + 6*r^3 + r)]
]);

//M := Matrix(F, [
//[1/49*(r^5 - r^3 - 13*r), 0, 1/49*(660*r^5 + 5808*r^3 + 10824*r)],
//[0, r, 0],
//1/539*(5*r^5 + 44*r^3 + 82*r), 0, 1/49*(6*r^5 + 43*r^3 + 69*r)]
//]);

//time test, fs := CantorMorphismFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
print D;

exit;
