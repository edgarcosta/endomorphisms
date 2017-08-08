AttachSpec("../../spec");
SetVerbose("EndoCheck", 0);

/* PLEASE NOTE: Give your elliptic curve as a HyperellipticCurve for now */

F := Rationals();
R<x> := PolynomialRing(F);
fX := x^6 - 8*x^4 + 2*x^3 + 16*x^2 - 36*x - 55;
X := HyperellipticCurve(fX);
P0 := X ! [1, 1, 0];
fY := x^3 + 215/3*x - 10582/27;
Y := HyperellipticCurve(fY);
Q0 := Y ! [1, 0, 0];
M := Matrix(F, [ [1, 1] ]);

time D := DivisorFromMatrix(X, P0, Y, Q0, M);
time D := DivisorFromMatrixSplit(X, P0, Y, Q0, M);
time fs := CantorFromMatrix(X, P0, Y, Q0, M);
time fs := CantorFromMatrixSplit(X, P0, Y, Q0, M);
print fs;

/* Check that the answer is a projection: */
R<x,y> := PolynomialRing(F, 2);
fX := R ! DefiningEquation(AffinePatch(X, 1)); fY := R ! DefiningEquation(AffinePatch(Y, 1));
IX := ideal<R | fX>;
print "Well-defined?", R ! Numerator(Evaluate(R ! fY, fs)) in IX;

/* Check that the action on differentials is correct: */
K := FieldOfFractions(R);
dx := K ! 1;
dy := K ! -Derivative(fX, 1)/Derivative(fX, 2);
ev := (Derivative(fs[1], 1)*dx + Derivative(fs[1], 2)*dy) / fs[2];
print "Correct pullback?", R ! Numerator(ev - (1 + x)/y) in IX;

AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));

exit;
