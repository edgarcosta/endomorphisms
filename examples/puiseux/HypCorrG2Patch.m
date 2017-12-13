AttachSpec("../../endomorphisms/magma/spec");
SetVerbose("EndoCheck", 3);

F := Rationals();
R<x> := PolynomialRing(F);
f := -x^5;
h := x^3 + x + 1;
X := HyperellipticCurve(f, h);
P0s := [ X ! [0, 0, 1], X ! [1, 0, 0] ];

T := Matrix(F, [
[ -1, -1],
[ -1,  0]
]);

for P0 in P0s do
    print "";
    print "Curve:";
    print X;
    print "Base point:";
    print P0;
    print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

    print "Calculating divisor:";
    time test, D := DivisorFromMatrixSplit(X, P0, X, P0, T : LowerBound := 6);
    eqs := DefiningEquations(D);
    R<y2,y1,x2,x1> := Parent(eqs[1]);
    print "Divisor:";
    print D;

    SetVerbose("EndoCheck", 3);
    print "Calculating Cantor representation...";
    time test, fs := CantorFromMatrixSplit(X, P0, X, P0, T : LowerBound := 16);
    R<x,y> := Parent(fs[1]);
    print "Cantor representation:";
    print fs;
end for;

exit;
