SetVerbose("EndoCheck", 0);

R<t> := PolynomialRing(Rationals());
F<r> := NumberField(t^2 - 2);
R<x> := PolynomialRing(F);
f := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1;
f := -f;
X := HyperellipticCurve(f);
P0 := X ! [1, 0, 0];
P0 := X ! [0, 1, 1];
Q0 := X ! [0, 1, 1];

T := Matrix(F, [
[ 1, -r],
[ -r, 1]
]);
T := -T;

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print T;

print "Curve:";
print X;

print "";
print "Calculating divisor:";
time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, T : LowerBound := 1);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

print "";
print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, X, Q0, T : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

fs := [ X`KU ! f : f in fs ];
ceqs := X`cantor_eqs;

print "";
print "Check 0:";
print [ Evaluate(ceq, fs) : ceq in ceqs ];

X := BaseExtend(X, X`KU);
R<x> := PolynomialRing(BaseRing(X));
J := Jacobian(X);

a := x^2 + fs[1]*x + fs[2];
b := fs[3]*x + fs[4];
div1 := J ! [a, b];

P0 := X ! [0, 1, 1];
P0m := X ! [0, -1, 1];
div0 := P0 - P0m;
print "Improved Cantor representation:";
print div1 - div0;
