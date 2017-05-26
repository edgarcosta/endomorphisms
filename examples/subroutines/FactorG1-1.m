AttachSpec("../../spec");

F := Rationals();
R<x> := PolynomialRing(F);
K<r> := NumberField(x^2 - 7);
R<x> := PolynomialRing(K);
fX := x^5 + x^4 + 2*x^3 + x^2 + x;
hX := x^2 + x;
X := HyperellipticCurve(fX, hX);
P0 := X ! [1, -r - 1, 1];
fY := 1/4*x^3 - 35/4*x - 49/2;
Y := HyperellipticCurve(fY);
Q0 := Y ! [1, 0, 0];

/* Check that the answer is a projection: */
R<x,y> := PolynomialRing(K, 2);
fs := [
(7*x^2 + 14*x + 7)/(x^2 - 2*x + 1),
(-14*x^2 - 14*x - 28*y)/(x^3 - 3*x^2 + 3*x - 1)
];
fX := R ! DefiningEquation(AffinePatch(X, 1)); fY := R ! DefiningEquation(AffinePatch(Y, 1));
IX := ideal<R | fX>;
print "Well-defined?", R ! Numerator(Evaluate(R ! fY, fs)) in IX;

/* Check that the action on differentials is correct: */
K := FieldOfFractions(R);
dx := K ! 1;
dy := K ! -Derivative(fX, 1)/Derivative(fX, 2);
ev := (Derivative(fs[1], 1)*dx + Derivative(fs[1], 2)*dy) / fs[2];
h := x^2 + x;
c := 2;
print "Correct pullback?", R ! Numerator(ev - c*(1 + x)/(2*y + h)) in IX;

AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));

exit;
