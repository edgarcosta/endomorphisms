SetVerbose("EndoFind", 3);
SetVerbose("CurveRec", 0);

prec := 100;

F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
f := -25*x^6 + 12*x^5 + 27*x^4 - 16*x^3 - 3*x^2 + 4*x + 1;
X := HyperellipticCurve(f);

print "";
print "Curve:";
print X;

desc := HeuristicEndomorphismAlgebra(X : CC := true);
rep := HeuristicEndomorphismRepresentation(X);
L := HeuristicEndomorphismFieldOfDefinition(X);
Lat := HeuristicEndomorphismLattice(X);
test_gl2 := HeuristicIsGL2(X);

print "";
print "Heuristic endomorphism lattice:";
print Lat;

/*
Should be
<[ <2, 1> ], <4, 2>, [
    <<1, 1>, [ 1, 0, 0, 0, 1 ], [ <2, 1> ], [
        <1, 4, [ -1, 1 ], 6, 2>
    ], <1, 1>, 3>,
    <<2, 1>, [ 1, 0, 1 ], [ <1, 1>, <1, 1> ], [
        <1, 2, [ -3, 0, 1 ], 1, 2>
    ], <1, -1>, 2>,
    <<2, 1>, [ 2, 0, 1 ], [ <1, 2> ], [
        <1, 2, [ 6, 0, 1 ], 1, 2>
    ], <1, -1>, 1>,
    <<2, 1>, [ -2, 0, 1 ], [ <1, 1>, <1, 1> ], [
        <1, 2, [ -2, 0, 1 ], 1, 2>
    ], <1, -1>, 2>,
    <<4, 2>, [ -1, 1 ], [ <1, 1> ], [
        <1, 1, [ -1, 1 ], 1, 2>
    ], <1, -1>, 1>
]>
*/

/*
facinfo := HeuristicJacobianFactors(X);
print "";
print "Heuristic Jacobian factors:";
print facinfo;

/*
exps, test, degs := IsogenyInformation(X : facinfo := facinfo);
print "";
print "Isogeny exponents:";
print exps;
print "";
print "Compatible with polarizations:";
print test;
print "";
print "Degrees (if applicable):";
print degs;
*/

