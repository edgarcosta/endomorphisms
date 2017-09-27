AttachSpec("../../spec");
SetVerbose("EndoCheck", 3);

F := Rationals();
P2<x,y,z> := ProjectiveSpace(F, 2);

fX := x^3*z + x^2*y*z + x^2*z^2 - x*y^3 + x*y^2*z + x*z^3 - y^2*z^2 + y*z^3;
X := Curve(P2, fX);
P0 := X ! [0, 1, 0];

M := Matrix(F, [
[ 0, 1,-1],
[ 0,-1, 0],
[-1,-1, 0]
]);

print "";
print "Field:";
print F;
print "Curve X:";
print X;
print "Point P0:";
print P0;
print "Tangent representation:";
print M;

print "Calculating Cantor representation...";
time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 20, Margin := 2^7);

eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

exit;
