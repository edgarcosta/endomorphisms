SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 1);

prec := 300;
CCSmall := ComplexField(5);
F := RationalsExtra(prec);
CC := F`CC;

/* All decomposing Sato--Tate groups are represented */
Xs := [ ];
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, -1, 0, 0, 0, 1], R![]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![-2, -8, -10, -5, 0, 1], R![0, 0, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![-2], R![0, 0, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, 0, 0, 0, 0, 0, 1], R![1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, 0, 0, 0, 0, 0, -1], R![1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![1, 3, 6, 7, 6, 3, 1], R![0, 1, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, -1, 0, -1], R![1, 1, 1, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, 0, 1, 2, 2, 1], R![1, 1, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, -1, 1, 1, -3, 2], R![1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, 0, 0, 0, 1, 1], R![1, 1, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![-1, 0, 2, 0, -2, 0, 1], R![]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, -4, 6, 0, -3, 1], R![]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![-1], R![1, 0, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, -1, 1, 0, -1, 1], R![]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![1, 0, 0, 1], R![1, 0, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![-7, 0, 0, 0, 1], R![0, 1, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, 0, 1, 0, 1], R![1, 0, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![1, 1, 2, 1, 1], R![1, 1, 1, 1]); Append(~Xs, X);

for X in Xs do
    print "Curve:";
    print X;
    print "";

    facs := HeuristicJacobianFactors(X : AllMaps := false);
    print "";
    print "Heuristic Jacobian factors:";
    print facs;
end for;
