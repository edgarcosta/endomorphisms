/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

//SetVerbose("EndoFind", 1);
//SetVerbose("CurveRec", 1);

prec := 500;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
R<x,y> := PolynomialRing(F, 2);
z := 1;

f := x^3*z + 2*x^2*y^2 + x^2*y*z + 2*x^2*z^2 - x*y^2*z + x*y*z^2 - x*z^3 + y^3*z - y^2*z^2 + y*z^3 - z^4;
f := x^3*z + x^2*y^2 - 3*x*y^2*z - 4*x*z^3 - 2*y^4 + y^3*z - 4*y^2*z^2 - 3*z^4;
f := x^3*z + 2*x^2*y^2 + x^2*y*z + 3*x^2*z^2 - 4*x*y^3 - x*y^2*z + 5*x*y*z^2 + x*z^3 + 2*y^4 + 6*y^3*z + 6*y^2*z^2 + 2*y*z^3;
f := 2*x^4 + 3*x^3*y + 4*x^3*z + 6*x^2*y^2 + 4*x^2*y*z + 7*x^2*z^2 + 4*x*y^3 + 4*x*y^2*z + 7*x*y*z^2 + 4*x*z^3 + 3*y^4 + 2*y^3*z + 3*y^2*z^2 + 5*y*z^3 + 2*z^4;

X := PlaneCurve(f);

P := PeriodMatrix(X);
EndoRep := GeometricEndomorphismRepresentation(P, F);

for idem in IdempotentsFromRepresentation(EndoRep) do
    print idem[1];
    Ys := DecompositionFactors(P, idem, F);
    print Ys;

    g := Rank(idem[2]) div 2;
    if g eq 1 then
        print [ jInvariant(EllipticCurve(HyperellipticPolynomials(Y))) : Y in Ys ];
    elif g eq 2 then
        print [* WPSNormalize([2, 4, 6, 8, 10], IgusaInvariants(Y)) : Y in Ys *];
    end if;
end for;
