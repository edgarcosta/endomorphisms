AttachSpec("../../spec");
SetVerbose("EndoCheck", 3);

F := Rationals();
R<t> := PolynomialRing(F);
F := SplittingField((t^2 - t + 2)*(t^2 - 21));
s := Roots(t^2 - t + 2, F)[1][1];
r := Roots(t^2 - 21, F)[1][1];

R<x> := PolynomialRing(F);
f := 5*x^6 - 4*x^5 - 5*x^4 + 14*x^3 - 5*x^2 - 4*x + 5;
h := x^3 + 1;
X := HyperellipticCurve(f, h);
P0 := X ! [ 1, (-r - 1)/2, 0 ];

/* This gives an obscene check: */
T := Matrix(F, [
[2*s - 1, 2*s - 1],
[2*s - 1, 2*s - 1]
]);

print "";
print "Curve:";
print X;
print "Base point:";
print P0;
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

print "Calculating divisor:";
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, T : LowerBound := 50);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

exit;
