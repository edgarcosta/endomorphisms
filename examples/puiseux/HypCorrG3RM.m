AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 1);
SetMemoryLimit(32*10^9);

R<t> := PolynomialRing(Rationals());
f := x^3 - x^2 - 2*x + 1;
F<r> := NumberField(f);
R<x> := PolynomialRing(F);
R<x> := PolynomialRing(Rationals());
f := x^7 - 14*x^6 + 210*x^5 - 658*x^4 + 245*x^3 + 588*x^2 + 637*x - 686;
X := HyperellipticCurve(f);
P0 := X ! [2, 18, 1];

M := Matrix(F, [
[-r^2 + 2*r + 1, 1/7*(-2*r^2 - 8*r + 6), 1/7*(2*r^2 - r - 3)],
[0, -3*r + 1, 0],
[28*r^2 - 14*r - 42, -10*r^2 - 4*r + 18, r^2 + r - 2]
]);

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;

print "Calculating divisor:";
time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 1);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

exit;
