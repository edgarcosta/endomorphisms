AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

F := Rationals();
P2<x,y,z> := ProjectiveSpace(F, 2);
/*
fX := -x^4 + x^3*z - 4*x^2*z^2 - x*y^2*z + 3*x*z^3 + y^4 - 3*z^4;
X := PlaneCurve(fX);
P0 := X ! [1, 1, 0];
*/
fX := -z^4 + z^3*y - 4*z^2*y^2 - z*x^2*y + 3*z*y^3 + x^4 - 3*y^4;
X := PlaneCurve(fX);
P0 := X ! [1, 0, 1];

S<t> := PolynomialRing(F);
fY := t^5 - 23/4*t^4 + t^3 + 33/2*t^2 + 12*t + 9/4;
Y := HyperellipticCurve(fY);
Q0 := Y ! [0, 3/2, 1];

/*
CC := ComplexFieldExtra(100);
P := PeriodMatrix([ ChangeRing(fX, CC) ], [ fX ] : HaveOldenburg := true);
Q := PeriodMatrix([ ChangeRing(fY, CC) ], [ fY ] : HaveOldenburg := true);
print GeometricIsogenyRepresentationPartial(P, Q);
print GeometricIsogenyRepresentationPartial(Q, P);
*/

MPQ := Matrix(F, [
[0,-2, 0],
[0, 0, 2]
]);
MQP := Matrix(F, [
[ 0, 0],
[ 1, 0],
[ 0,-1]
]);

/*
MPQ := Matrix(F, [
[0, 0,-2],
[2, 0, 0]
]);
MQP := Matrix(F, [
[ 0, 1],
[ 0, 0],
[-1, 0]
]);
*/

print "Calculating Cantor representation...";
time test, D := DivisorFromMatrixSplit(X, P0, Y, Q0, MPQ : LowerBound := 30);
print D;

exit;

print "Calculating Cantor representation...";
time test, D := DivisorFromMatrixSplit(Y, Q0, X, P0, MQP : LowerBound := 1);
print D;

exit;
