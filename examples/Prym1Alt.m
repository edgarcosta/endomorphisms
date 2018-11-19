SetVerbose("EndoCheck", 3);

prec := 700;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
F<r> := NumberFieldExtra(x^2 - 2*7*17*19);
R<x> := PolynomialRing(F);

f := (-7 + x)*(-5 + x)*(4 + x)*(8 + x)*(17 + x)*(19 + x)*(20 + x); h := R ! 0;
X := HyperellipticCurve(f, h);
P0 := X ! [-4, 0, 1];
//P0 := X ! [0, 2^3*5*r, 1];

f := x^3 - 19331568*x - 9837700992; h := R ! 0;
E := HyperellipticCurve(f, h);
O0 := E ! [1, 0, 0];

PX := PeriodMatrix(X);
PE := PeriodMatrix(E);

HomRep := GeometricHomomorphismRepresentation(PX, PE, F);
T := HomRep[1][1];

print "";
print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, E, O0, T);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;
print CorrespondenceVerifyG1(X, E, T, fs : CheckDegree := true);
