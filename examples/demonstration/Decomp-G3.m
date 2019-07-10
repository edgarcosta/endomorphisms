SetSeed(1);
SetVerbose("EndoFind", 3);
SetVerbose("CurveRec", 2);

prec := 333;

F := RationalsExtra(prec);
R<x,y,z> := PolynomialRing(F, 3);

f := x^3*z + 2*x^2*y^2 + x^2*y*z + 2*x^2*z^2 - x*y^2*z + x*y*z^2 - x*z^3 + y^3*z - y^2*z^2 + y*z^3 - z^4;
X := PlaneCurve(f);

print "";
print "Curve:";
print X;

L := HeuristicEndomorphismFieldOfDefinition(X);
print "";
print "Heuristic field of definition of the endomorphisms:";
print L;

A, desc := HeuristicEndomorphismAlgebra(X : Geometric := true);
print "";
print "Heuristic geometric endomorphism algebra:";
print A;

facs := HeuristicJacobianFactors(X);
print "";
print "Heuristic Jacobian factors:";
print facs;
