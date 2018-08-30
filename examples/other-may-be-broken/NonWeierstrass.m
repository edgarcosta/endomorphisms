AttachSpec("../../endomorphisms/magma/spec");

F := QQ;
R<x> := PolynomialRing(F);
f := 15*x^5 + 50*x^4 + 55*x^3 + 22*x^2 + 3*x; h := x;
f := x^6 + x + 1; h := 0;
X := HyperellipticCurve(f, h);

print NonWeierstrassBasePoint(X, F);

exit;
