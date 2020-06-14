SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 50;

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

