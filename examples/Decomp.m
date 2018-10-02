/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

AttachSpec("../endomorphisms/magma/spec");
SetVerbose("EndoFind", 0);

prec := 500;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := x^8 + x^6 + 1;
X := HyperellipticCurve(f);

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
        print [ WPSNormalize([2, 4, 6, 8, 10], IgusaInvariants(Y)) : Y in Ys ];
    end if;
end for;
