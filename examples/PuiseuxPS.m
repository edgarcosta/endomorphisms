AttachSpec("../spec");
SetVerbose("EndoCheck", 1);
SetMemoryLimit(32*10^9);

F := Rationals();
P2<x0,x1,x2> := ProjectiveSpace(F, 2);
f := x1^3*x2 - (8*x0^4 + 8*x0^2*x2^2 + x2^4);
X := Curve(P2, f);
//print PointSearch(X, 10);
P0 := X ! [0, 1, 1];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [[0,2,0],[1,0,0],[0,0,0]]);
T := Matrix(F, [[0,0,1],[0,1,0],[1,0,0]]);
M := T * M * T;
M := Transpose(M);

time test, D := DivisorFromMatrixSplit(X, P0, X, P0, M : LowerBound := 1);
print D;
print IsIrreducible(D);

Is := IrreducibleComponents(D);
D := Is[1];

eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
I := EliminationIdeal(ideal<R | eqs>, { y1, y2 });
I := EliminationIdeal(ideal<R | eqs>, { x1, x2 });
test, gen := IsPrincipal(I);
print test;
print gen;

exit;
