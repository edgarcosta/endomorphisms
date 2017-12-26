AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

F := Rationals();
P2<x,y,z> := ProjectiveSpace(F, 2);

fX := x^3*y + x^3*z + x^2*y^2 + 3*x^2*y*z + x^2*z^2 - 4*x*y^3 - 3*x*y^2*z - 3*x*y*z^2 - 4*x*z^3 + 2*y^4 + 3*y^2*z^2 + 2*z^4;
X := Curve(P2, fX);
P0 := X ! [1, 1, 1];

M := Matrix(F, [
[-1, 0, 0],
[ 0, 0,-1],
[ 0,-1, 0]
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
time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, M : LowerBound := 1);
print D;

exit;
