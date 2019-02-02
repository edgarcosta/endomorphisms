/***
 *  An example from the paper
 */

AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

F := Rationals();
R<x> := PolynomialRing(F);

p := x^8 - 12*x^7 + 50*x^6 - 108*x^5 + 131*x^4 - 76*x^3 - 10*x^2 + 44*x - 19;
X := HyperellipticCurve(p);
P0 := X ! [1, -1, 0];

q := x^3 + 6656/3*x - 185344/27;
Y := HyperellipticCurve(q);
Q0 := Y ! [1, 0, 0];

M := Matrix(F, [
[ 1, -1,  1/2 ]
]);

print "Field:";
print F;
print "Curves:";
print X;
print Y;
print "Tangent representation:";
print M;
print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, Y, Q0, M : LowerBound := 16);
K<x,y> := Parent(fs[1]);
print fs;

/* Check that the answer is a projection: */
R<x,y> := PolynomialRing(F, 2);
K := FieldOfFractions(R);
fs := [ (64*x^2 - 400*x + 272)/(3*(x^2 - x - 1)) , 224*y/(x^2 - x - 1)^2 ];
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
