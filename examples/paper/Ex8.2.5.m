/***
 *  An example from the paper
 */

SetVerbose("EndoCheck", 3);

/*
NOTE: For some reason this code fails currently. We are fixing this.
For a working version that gives the corresponding divisor, use the git checkout
with sha e88bacc. The relevant file has the same name.
*/

R<t> := PolynomialRing(Rationals());
f := t^3 - t^2 - 2*t + 1;
F<zeta7> := NumberField((t^7 - 1) div (t - 1));
r := Roots(f, F)[1][1];
F<r> := NumberField(f);

P2<x0,x1,x2> := ProjectiveSpace(F, 2);
f := x0^4 + 8*x0^3*x2 + 2*x0^2*x1*x2 + 25*x0^2*x2^2 - x0*x1^3 + 2*x0*x1^2*x2 + 8*x0*x1*x2^2 + 36*x0*x2^3 + x1^4 - 2*x1^3*x2 + 5*x1^2*x2^2 + 9*x1*x2^3 + 20*x2^4;
X := Curve(P2, f);
P0 := X ! [-2, 0, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[          r^2 - 2,               0, 2*r^2 + 2*r - 4],
[                0,    -r^2 + r + 1,               0],
[                0,               0,              -r]
]);

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;

print "Calculating divisor:";
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

exit;
