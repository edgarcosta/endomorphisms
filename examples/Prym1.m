SetVerbose("EndoCheck", 1);

prec := 100;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

f := x^8 + x^6 + 1; h := R ! 0;
X := HyperellipticCurve(f, h);
P0 := X ! [0, 1, 1];

f := x^3 - 64*x + 64; h := R ! 0;
E := HyperellipticCurve(f, h);
O0 := E ! [1, 0, 0];

PX := PeriodMatrix(X);
PE := PeriodMatrix(E);

HomRep := GeometricHomomorphismRepresentation(PX, PE, F);
T := HomRep[1][1];

print "";
print "Calculating divisor:";
time test, D := DivisorFromMatrixAmbientSplit(X, P0, E, O0, T);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

print "";
print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, E, O0, T);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;
/* TODO: Sign error */
print CorrespondenceVerifyG1(X, E, -T, fs : CheckDegree := true);
