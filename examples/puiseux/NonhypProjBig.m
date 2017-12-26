AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

/* An absolutely obscene example */

F := Rationals();
R<t> := PolynomialRing(F);
F<r> := NumberField(t^2 - t + 1);
P2<x,y,z> := ProjectiveSpace(F, 2);

fX := x^3*z + x^2*z^2 + x*y^2*z + x*y*z^2 + x*z^3 - y^4 - 2*y^3*z - 2*y^2*z^2 - y*z^3;
X := Curve(P2, fX);
P0 := X ! [0, -r, 1];

S<t> := PolynomialRing(F);
fY := t^3 - 675*t - 27675/4;
fY := t^3 - 432/625*t + 556416/15625;
fY := t^3 - 432*t + 556416;
Y := HyperellipticCurve(fY);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [ [ 1/3, 0, 1/3 ] ]);
M := Matrix(F, [ [ 5/12, -5/12, -5/12 ] ]);
M := Matrix(F, [ [ 1/4, 5/12, 1/12 ] ]);

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
time test, fs := CantorFromMatrixAmbientSplit(X, P0, Y, Q0, M : LowerBound := 1);
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
