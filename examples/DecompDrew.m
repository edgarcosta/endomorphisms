/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

//SetVerbose("EndoFind", 1);
//SetVerbose("CurveRec", 1);

prec := 600;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := 2*x^10 + 6*x^9 + 6*x^8 + 12*x^7 + 7*x^6 + 7*x^4 - 12*x^3 + 6*x^2 - 6*x + 2; h := R ! 0;
//f := 10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10; h := R ! 0;

X := HyperellipticCurve(f, h);
X := ReducedMinimalWeierstrassModel(X);
print "Curve:";
print X;

time P := PeriodMatrix(X);
time GeoEndoRep := GeometricEndomorphismRepresentation(P, F);
print EndomorphismStructureFromRepresentation(GeoEndoRep);

/*
for idem in IdempotentsFromRepresentation(EndoRep) do
    print idem[2];
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
*/
