load "../../puiseux/FractionalCRT.m";

R<x> := PolynomialRing(Rationals());
K<r> := NumberField(x^2 - 2);
R<x> := PolynomialRing(K);
L<s> := NumberField(x^3 + x + r);
R<x> := PolynomialRing(L);

print R;
print L;
print RandomSplitPrime(x^2 - 5, 10);

X := HyperellipticCurve([ x^5 + x + s, x ]);
print X;
print AbsoluteField(L);
print ChangeRing(X, AbsoluteField(L));

exit;
