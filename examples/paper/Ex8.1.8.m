/***
 *  An example from the paper
 */

SetVerbose("EndoCheck", 3);

F := Rationals();
R<x> := PolynomialRing(F);

p := x^6 - 8*x^4 + 2*x^3 + 16*x^2 - 36*x - 55;
X := HyperellipticCurve(p);
P0 := X ! [1, -1, 0];

q := x^3 + 3440/3*x - 677248/27;
Y := HyperellipticCurve(q);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [
[ -1/2, -1/2 ]
]);

print "Field:";
print F;
print "Curves:";
print X;
print Y;
print "Tangent representation:";
print M;
print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 16);
K<x,y> := Parent(fs[1]);
print fs;

/* Check that the answer is a projection: */
R<x,y> := PolynomialRing(F, 2);
K := FieldOfFractions(R);
fX := R ! DefiningEquation(AffinePatch(X, 1)); fY := R ! DefiningEquation(AffinePatch(Y, 1));
IX := ideal<R | fX>;
print "Well-defined?", R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX;

/* Check that the action on differentials is correct: */
fX := R ! DefiningEquation(AffinePatch(X, 1));
fY := R ! DefiningEquation(AffinePatch(Y, 1));
dx := K ! 1;
dy := K ! -Derivative(fX, 1)/Derivative(fX, 2);
ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! -Evaluate(Derivative(fY, 2)/2, fs));

print "Correct pullback?", R ! Numerator(K ! (ev - &+[ M[1,i]*x^(i - 1) : i in [1..Genus(X)] ]/y)) in IX;

AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));

F := Rationals();
R<x> := PolynomialRing(F);

p := x^6 - 8*x^4 + 2*x^3 + 16*x^2 - 36*x - 55;
X := HyperellipticCurve(p);
P0 := X ! [1, -1, 0];

q := x^3 + 752/3*x - 9088/27;
Y := HyperellipticCurve(q);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [
[ 3/2, 1/2 ]
]);

print "Field:";
print F;
print "Curves:";
print X;
print Y;
print "Tangent representation:";
print M;
print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 16);
K<x,y> := Parent(fs[1]);
print fs;

/* Check that the answer is a projection: */
R<x,y> := PolynomialRing(F, 2);
K := FieldOfFractions(R);
fX := R ! DefiningEquation(AffinePatch(X, 1)); fY := R ! DefiningEquation(AffinePatch(Y, 1));
IX := ideal<R | fX>;
print "Well-defined?", R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX;

/* Check that the action on differentials is correct: */
fX := R ! DefiningEquation(AffinePatch(X, 1));
fY := R ! DefiningEquation(AffinePatch(Y, 1));
dx := K ! 1;
dy := K ! -Derivative(fX, 1)/Derivative(fX, 2);
ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! -Evaluate(Derivative(fY, 2)/2, fs));

print "Correct pullback?", R ! Numerator(K ! (ev - &+[ M[1,i]*x^(i - 1) : i in [1..Genus(X)] ]/y)) in IX;

AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));

exit;
