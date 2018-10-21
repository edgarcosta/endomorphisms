/*
  An example in Magma (to be made into an intuitve package there as well).
  Examples of verifications and projections can be found in the puiseux/
  directory; this file shows how to access the heuristic part.
*/

SetVerbose("EndoFind", 1);
SetVerbose("CurveRec", 1);

prec := 600;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

R<x> := PolynomialRing(F);
f := x^8 + x^6 + 2; h := R ! 0;
f := (-7 + x)*(-5 + x)*(4 + x)*(8 + x)*(17 + x)*(19 + x)*(20 + x); h := R ! 0;
/* This one gives error in reduction step, contradicting Dupont */
//f := x^7 + x^6 + 5*x^5 - 3*x^4 + 2*x^3 - 13*x^2 + 7*x - 1; h := x^3 + x;
//f := x^7 + x^6 + x^5 + x^3 + x^2 + x; h := x^4 + x^2 + 1;
//f := -2*x^7 - 4*x^6 + 3*x^4 + x^3 - 2*x^2 - x; h := x^2 + x + 1;
//f := x^7 - 2*x^5 - 4*x^4 - 2*x^3 + x; h := x^4 + x^2 + 1;
//f := x^7 - 4*x^6 - x^5 + 10*x^4 + 3*x^3 - 6*x^2 - x + 1; h := x^3 + x;
//f := x^7 + x^6 - 4*x^5 + 4*x^3 - 5*x^2 + 2*x - 1; h := x^4 + x^3 + x + 1; SetEpsComp(F`CC, 10^70*F`CC`epscomp);

X := HyperellipticCurve(f, h);
X := ReducedMinimalWeierstrassModel(X);
print "Curve:";
print X;

P := PeriodMatrix(X);
EndoRep := GeometricEndomorphismRepresentation(P, F);

for idem in IdempotentsFromRepresentation(EndoRep) do
    print idem[2];
    Ys := DecompositionFactors(P, idem, F);
    print Ys;

    g := Rank(idem[2]) div 2;
    if g eq 1 then
        print [ jInvariant(EllipticCurve(HyperellipticPolynomials(Y))) : Y in Ys ];
    elif g eq 2 then
        print [* WPSNormalize([2, 4, 6, 8, 10], IgusaInvariants(Y)) : Y in Ys *];
    end if;
end for;
