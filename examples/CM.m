SetVerbose("EndoCheck", 3);

prec := 100;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
F := NumberFieldExtra(x^6 + 42*x^4 + 441*x^2 + 784);
R<x> := PolynomialRing(F);

f := 16*x^7 + 357*x^5 - 819*x^3 + 448*x;
f /:= 2; h := R ! 0;
X := HyperellipticCurve(f, h);
// Try to make one the point at infinity!
P0 := X ! [1, 1, 1];

Y := HyperellipticCurve(f, h);
Q0 := Y ! [1, 1, 1];

PX := PeriodMatrix(X);
PY := PeriodMatrix(Y);

HomRep := GeometricHomomorphismRepresentation(PX, PY, F);
T := HomRep[3][1];
print T;

/*
This seems to work without needing a base field extension!
Also, we run into the copy problem
Finally, there is a problem in that the first iteration does not stabilize
*/

print "";
print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, Y, Q0, T);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

fs := [ X`KU ! f : f in fs ];
ceqs := Y`cantor_eqs;

print "";
print "Check 0:";
print [ Evaluate(ceq, fs) : ceq in ceqs ];

Y := BaseExtend(Y, X`KU);
R<x> := PolynomialRing(BaseRing(Y));
J := Jacobian(Y);

a := x^2 + fs[1]*x + fs[2];
b := fs[3]*x + fs[4];
div1 := J ! [a, b];

P := Y ! [1, 0, 1];
Q0 := Y ! [0, 3*r, 1];
div0 := P - Q0;
print "";
print "Improved Cantor representation, version 1:";
print div1 + 2*div0;

/* This always works, also without a rational Weierstrass point: */
Q0m := Y ! [0, -3*r, 1];
div0 := Q0 - Q0m;
print "Improved Cantor representation, version 2:";
print div1 - div0;

