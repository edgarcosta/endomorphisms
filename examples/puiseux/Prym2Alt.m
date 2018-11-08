SetVerbose("EndoCheck", 3);

prec := 1000;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
f := x^8 + 7348293352896*x^6 + 12007076842932839406124542*x^4 + 2046917143744335830152529873759073348*x^2 + 4995663275630324728859423475246999393966892466625;
f := Polredbestabs(f);
g := x^4 - 2*x^3 - 120*x^2 - 8*x - 1196;
F<r> := NumberFieldExtra(g);
R<x> := PolynomialRing(F);
r := RootsPari(g, F)[1];

f := (-7 + x)*(-5 + x)*(4 + x)*(8 + x)*(17 + x)*(19 + x)*(20 + x); h := R ! 0;
X := HyperellipticCurve(f, h);
P0 := X ! [-4, 0, 1];
//P0 := X ! [0, 2^3*5*r, 1];

f := 1/37688194795*(-49049*r^3 + 1122882*r^2 + 395793*r + 51067324)*x^6 + 1/7537638959*(-5798660*r^3 - 29837904*r^2 + 113495376*r - 5768592683)*x^5 + 1/37688194795*(11317747532*r^3 - 4123812036*r^2 - 107423158224*r + 6335026039013)*x^4 + 1/443390527*(-9827357600*r^3 - 7770236640*r^2 - 300045996096*r - 8195152652288)*x^3 + 1/37688194795*(-43933484716272*r^3 + 249970367739456*r^2 + 7566324780761424*r + 43296433438988772)*x^2 + 1/7537638959*(980527002798528*r^3 - 3153589821468288*r^2 - 106968834815231232*r - 222415949032570944)*x + 1/37688194795*(-93107780429282496*r^3 + 50938205072753088*r^2 + 13766602640796312192*r - 4928358464914780224);
h := R ! 0;
Y := HyperellipticCurve(f, h);
Q0 := Y ! [2, Roots(x^2 - Evaluate(f, 2))[1][1], 1];

/*
PX := PeriodMatrix(X);
PY := PeriodMatrix(Y);

HomRep := GeometricHomomorphismRepresentation(PX, PY, F);
T := HomRep[1][1];
*/

T := Matrix(F, [[ 1,  0,  0], [94, 19,  1]]);
print T;

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

Q0m := Y ! [0, -Q0[2], 1];
div0 := Q0 - Q0m;
print "Improved Cantor representation, version 2:";
print div1 - div0;


