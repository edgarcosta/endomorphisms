SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 1);

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

/* All decomposing Sato--Tate groups are represented */
Cs := [ ];
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, -1, 0, 0, 0, 1], R![]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![-2, -8, -10, -5, 0, 1], R![0, 0, 0, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![-2], R![0, 0, 0, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, 0, 0, 0, 0, 0, 1], R![1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, 0, 0, 0, 0, 0, -1], R![1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![1, 3, 6, 7, 6, 3, 1], R![0, 1, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, -1, 0, -1], R![1, 1, 1, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, 0, 1, 2, 2, 1], R![1, 1, 0, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, -1, 1, 1, -3, 2], R![1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, 0, 0, 0, 1, 1], R![1, 1, 0, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![-1, 0, 2, 0, -2, 0, 1], R![]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, -4, 6, 0, -3, 1], R![]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![-1], R![1, 0, 0, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, -1, 1, 0, -1, 1], R![]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![1, 0, 0, 1], R![1, 0, 0, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![-7, 0, 0, 0, 1], R![0, 1, 0, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, 0, 1, 0, 1], R![1, 0, 0, 1]); Append(~Cs, C);
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![1, 1, 2, 1, 1], R![1, 1, 1, 1]); Append(~Cs, C);

for C in Cs do
    print "Curve:";
    print C;
    print "";

    P := PeriodMatrix(C);
    GeoEndoRep := GeometricEndomorphismRepresentation(P, F);

    comps := SplitComponents(P, GeoEndoRep);
    print "Components:";
    for comp in comps do
        Q, mor := Explode(comp);
        recs := ReconstructionsFromComponent(P, Q, mor);
        print recs[1];
        print "";
    end for;
end for;
