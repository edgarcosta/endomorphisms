/***
 *  Iterator inspection
 */


AttachSpec("../../endomorphisms/magma/spec");
import "../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve, ChangeTangentAction;
import "../../endomorphisms/magma/puiseux/Branches.m": PuiseuxRamificationIndex, InitializeLift, CreateLiftIterator;
SetVerbose("EndoCheck", 4);


F := FiniteField(4001);
R<x> := PolynomialRing(F);
r := F ! 2845;
i := F ! 3102;

p := x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1;
X := HyperellipticCurve(p);
P0 := X ! [0, i];
print "Check that base point is not Weierstrass:", not IsWeierstrassPlace(Place(P0));

M := Matrix(F, [
[ 0, r ],
[ r, 0 ]
]);

print "Field:";
print F;
print "Curve:";
print X;
print "Tangent representation:";
print M;

InitializeCurve(X, P0);
M := ChangeTangentAction(X, X, M);
e := PuiseuxRamificationIndex(M);
P, Qs := InitializeLift(X, X, M);
IterateLift := CreateLiftIterator(X, X, M);

// TODO: Higher gives a denominator!
n := 2400;
n := 100;
for i:=1 to Ceiling(Log(2, n + e + 1)) - 1 do
    P_old := P; Qs_old := Qs;
    P, Qs := IterateLift(P_old, Qs_old, n + e + 1);
    print P[1] - P_old[1];
    print P[2] - P_old[2];
    for k in [1..#Qs] do
        print Qs[k][1] - Qs_old[k][1];
        print Qs[k][2] - Qs_old[k][2];
    end for;
    print Qs[1][1];
    print Qs[1][2];
end for;

exit;
