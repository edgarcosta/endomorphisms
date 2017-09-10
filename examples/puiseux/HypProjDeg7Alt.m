AttachSpec("../../spec");
SetVerbose("EndoCheck", 3);

F := Rationals();
R<x> := PolynomialRing(F);
pX := 56*x^6 + 45*x^4 + 30*x^3 + 9*x^2 + 12*x + 4;
X := HyperellipticCurve(pX);
P0 := X ! [0, 2, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

pY := 142*x^4 - 51*x^3 + x;
Y := HyperellipticCurve(pY);
Q0 := Y ! [0, 0, 1];

CC := ComplexFieldExtra(300);
PiX := Transpose(Matrix(CC, PeriodMatrix(pX : Prec := Precision(CC)))) / 2;

