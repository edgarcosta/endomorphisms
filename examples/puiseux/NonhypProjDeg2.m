AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 0);

F := Rationals();
P2<x,y,z> := ProjectiveSpace(F, 2);

fX := x^4 + 2*x^3*y + 2*x^3*z - 4*x^2*y^2 + 2*x^2*y*z - 4*x^2*z^2 - x*y^3 - x*z^3 + 2*y^4 - 3*y^3*z + 5*y^2*z^2 - 3*y*z^3 + 2*z^4;
X := Curve(P2, fX);
P0 := X ! [1, 0, 1];

S<t> := PolynomialRing(F);
fY := t^3 - 241/3*t + 7378/27;
Y := HyperellipticCurve(fY);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [ [0, 1/2, -1/2 ] ]);

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
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print fs;

R<x,y> := PolynomialRing(F, 2);
K := FieldOfFractions(R);
print "Original equation:", fX;
print "First affine patch:", R ! DefiningEquation(AffinePatch(X, 1));

/* Check that the answer is a projection after transformation, automatic: */
fX := R ! DefiningEquation(AffinePatch(X, 1)); fY := R ! DefiningEquation(AffinePatch(Y, 1));
IX := ideal<R | fX>;
print "Well-defined?", R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX;

/*
AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
KX := FunctionField(AX); KY := FunctionField(AY);
m := map<AX -> AY | fs >;
print "Degree:", Degree(ProjectiveClosure(m));
*/

/* Check that the action on differentials is correct: */
fX := R ! DefiningEquation(AffinePatch(X, 1));
fY := R ! DefiningEquation(AffinePatch(Y, 1));
dx := K ! 1;
dy := K ! -Derivative(fX, 1)/Derivative(fX, 2);
ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! (2*fs[2]));
mults := [ x, y, 1 ];
print "Correct pullback?", R ! Numerator(K ! (ev - &+[ M[1,i]*mults[i] : i in [1..Genus(X)] ]/Derivative(fX, 2))) in IX;

exit;
