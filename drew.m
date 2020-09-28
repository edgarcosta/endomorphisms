//SetVerbose("EndoFind", 3);

F := RationalsExtra(1000);
R<x> := PolynomialRing(F);

f := 5*x^6-15*x^5+25*x^4-25*x^2+15*x-5;
f := 4*x^6-9*x^5+20*x^4-20*x^2+9*x-4;
C := HyperellipticCurve(f);

print HeuristicEndomorphismAlgebra(C);
print HeuristicEndomorphismAlgebra(C : Geometric := true); 


