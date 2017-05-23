AttachSpec("../spec");
SetVerbose("EndoCheck", 1);

F := Rationals();
R<x> := PolynomialRing(F);
f := x^8 - 12*x^7 + 50*x^6 - 108*x^5 + 131*x^4 - 76*x^3 - 10*x^2 + 44*x - 19;
X := HyperellipticCurve(f);
P0 := X ! [1, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[-1, -2, -4],
[ 2,  3,  4],
[-1, -1, -1]
]);
M := Transpose(M);

time test, fs := CantorMorphismFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
print D;

/*
DEs := DefiningEquations(D);
DEs := [ 2*30697399*DEs[1], 2*30697399*DEs[2] ];
Lat := Lattice(Matrix([ Coefficients(DE) : DE in DEs ]));
print LLL(Lat);
*/

exit;

// TODO: What follows worked much better with the old, less general, Puiseux
// folder in which DivisorFromMatrixSplit has three arguments. Find out why.

R<t> := PolynomialRing(Rationals());
f := t^3 - t^2 - 2*t + 1;
F<zeta7> := NumberField((t^7 - 1) div (t - 1));
r := Roots(f, F)[1][1];
P2<x0,x1,x2> := ProjectiveSpace(F, 2);
f := x0^4 + 8*x0^3*x2 + 2*x0^2*x1*x2 + 25*x0^2*x2^2 - x0*x1^3 + 2*x0*x1^2*x2 + 8*x0*x1*x2^2 + 36*x0*x2^3 + x1^4 - 2*x1^3*x2 + 5*x1^2*x2^2 + 9*x1*x2^3 + 20*x2^4;
X := Curve(P2, f);
//print PointSearch(X, 10);
P0 := X ! [-2, 0, 1];
//P0 := X ! [-3, -2, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[   -r^2 + r + 1,               0, 2*r^2 - 4*r - 2],
[              0,         r^2 - 2,               0],
[              0,               0,              -r]
]);
T := Matrix(F, [[0,0,1],[0,1,0],[1,0,0]]);
M := T * M * T;
M := Transpose(M);

time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
time test, fs := CantorMorphismFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);

exit;
