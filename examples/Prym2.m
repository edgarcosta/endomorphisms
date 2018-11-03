SetVerbose("EndoCheck", 1);

prec := 100;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
F<r> := NumberFieldExtra(x^2 + 1);
R<x> := PolynomialRing(F);

f := x^8 + x^6 + 1; h := R ! 0;
X := HyperellipticCurve(f, h);
P0 := X ! [0, 1, 1];

f := 8*x^6 - 8*x^5 - 24*x^4 + 16*x^3 + 25*x^2 - 8*x - 9; h := R ! 0;
Y := HyperellipticCurve(f, h);
Q0 := Y ! [0, 3*r, 1];

PX := PeriodMatrix(X);
PY := PeriodMatrix(Y);

HomRep := GeometricHomomorphismRepresentation(PX, PY, F);
print HomRep;
T := HomRep[1][1];

print "";
print "Calculating divisor:";
time test, D := DivisorFromMatrixAmbientSplit(X, P0, Y, Q0, T);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print GroebnerBasis(Ideal(D));

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
P := Y ! [1, 0, 1];
Q0 := Y ! [0, 3*r, 1];
div0 := P - Q0;

a := x^2 + fs[1]*x + fs[2];
b := fs[3]*x + fs[4];
div1 := J ! [a, b];

print "";
print "Improved Cantor representation:";
print div1 + 2*div0;
