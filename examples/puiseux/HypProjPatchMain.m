AttachSpec("../../spec");
SetVerbose("EndoCheck", 1);

F := Rationals();
R<x> := PolynomialRing(F);
fX := (x^2 + x)^3 + (x^2 + x) + 1;
X := HyperellipticCurve(fX);
P0 := X ! [0, 1, 1];

fY := x^3 + x + 1;
Y := HyperellipticCurve(fY);
Q0 := Y ! [0, 1, 1];

M := Matrix(F, [
[1, 2]
]);

print "Field:";
print F;
print "Curve X:";
print X;
print "Point P0:";
print P0;
print "Curve Y:";
print Y;
print "Point Q0:";
print Q0;
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(Q0));
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
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
ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! (2*fs[2]));
print "Correct pullback?", R ! Numerator(K ! (ev - &+[ M[1,i]*x^(i - 1) : i in [1..Genus(X)] ]/(2*y))) in IX;

AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));


F := Rationals();
R<x> := PolynomialRing(F);
fX := (x^2 + x)^3 + (x^2 + x) + 1;
X := HyperellipticCurve(fX);
P0 := X ! [0, 1, 1];

fY := x^4 + x^3 + x;
Y := HyperellipticCurve(fY);
Q0 := Y ! [0, 0, 1];

M := Matrix(F, [
[1, 2]
]);

print "";
print "Field:";
print F;
print "Curve X:";
print X;
print "Point P0:";
print P0;
print "Curve Y:";
print Y;
print "Point Q0:";
print Q0;
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(Q0));
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
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
ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! (2*fs[2]));
print "Correct pullback?", R ! Numerator(K ! (ev - &+[ M[1,i]*x^(i - 1) : i in [1..Genus(X)] ]/(2*y))) in IX;

AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));


F := Rationals();
R<x> := PolynomialRing(F);
fX := (x^2 + x)^3 + (x^2 + x) + 1;
fX := R ! (x^6*Evaluate(fX, 1/x));
X := HyperellipticCurve(fX);
P0 := X ! [0, 1, 1];

fY := x^3 + x + 1;
Y := HyperellipticCurve(fY);
Q0 := Y ! [0, 1, 1];

M := Matrix(F, [
[2, 1]
]);

print "";
print "Field:";
print F;
print "Curve X:";
print X;
print "Point P0:";
print P0;
print "Curve Y:";
print Y;
print "Point Q0:";
print Q0;
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(Q0));
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
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
ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! (2*fs[2]));
print "Correct pullback?", R ! Numerator(K ! (ev - &+[ M[1,i]*x^(i - 1) : i in [1..Genus(X)] ]/(2*y))) in IX;

AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));


F := Rationals();
R<x> := PolynomialRing(F);
fX := (x^2 + x)^3 + (x^2 + x) + 1;
fX := R ! (x^6*Evaluate(fX, 1/x));
X := HyperellipticCurve(fX);
P0 := X ! [0, 1, 1];

fY := x^4 + x^3 + x;
Y := HyperellipticCurve(fY);
Q0 := Y ! [0, 0, 1];

M := Matrix(F, [
[2, 1]
]);

print "";
print "Field:";
print F;
print "Curve X:";
print X;
print "Point P0:";
print P0;
print "Curve Y:";
print Y;
print "Point Q0:";
print Q0;
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(Q0));
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
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
ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! (2*fs[2]));
print "Correct pullback?", R ! Numerator(K ! (ev - &+[ M[1,i]*x^(i - 1) : i in [1..Genus(X)] ]/(2*y))) in IX;

AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));

exit;
