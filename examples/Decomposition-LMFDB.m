SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);

prec := 100;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

/* All decomposing Sato--Tate groups are represented */
Xs := [ ];
X := HyperellipticCurve(R![0, -1, 0, 0, 0, 1], R![]); Append(~Xs, X);
X := HyperellipticCurve(R![-2, -8, -10, -5, 0, 1], R![0, 0, 0, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![-2], R![0, 0, 0, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![0, 0, 0, 0, 0, 0, 1], R![1]); Append(~Xs, X);
X := HyperellipticCurve(R![0, 0, 0, 0, 0, 0, -1], R![1]); Append(~Xs, X);
X := HyperellipticCurve(R![1, 3, 6, 7, 6, 3, 1], R![0, 1, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![0, -1, 0, -1], R![1, 1, 1, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![0, 0, 1, 2, 2, 1], R![1, 1, 0, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![0, -1, 1, 1, -3, 2], R![1]); Append(~Xs, X);
X := HyperellipticCurve(R![0, 0, 0, 0, 1, 1], R![1, 1, 0, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![-1, 0, 2, 0, -2, 0, 1], R![]); Append(~Xs, X);
X := HyperellipticCurve(R![0, -4, 6, 0, -3, 1], R![]); Append(~Xs, X);
X := HyperellipticCurve(R![-1], R![1, 0, 0, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![0, -1, 1, 0, -1, 1], R![]); Append(~Xs, X);
X := HyperellipticCurve(R![1, 0, 0, 1], R![1, 0, 0, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![-7, 0, 0, 0, 1], R![0, 1, 0, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![0, 0, 1, 0, 1], R![1, 0, 0, 1]); Append(~Xs, X);
X := HyperellipticCurve(R![1, 1, 2, 1, 1], R![1, 1, 1, 1]); Append(~Xs, X);

/* A genus-3 case */
X := HyperellipticCurve(x^8 + x^6 + 5*x^4 - 3*x^2 + 17); Append(~Xs, X);

//Xs := [ Xs[18] ];

for i in [1..#Xs] do
    X := Xs[i];

    print "";
    print "";
    print "Curve number", i, ":";
    print X;

    print "";
    print "Lattice:";
    print HeuristicEndomorphismLattice(X);

    print "";
    print "Decomposition:";
    print HeuristicDecompositionDescription(X);
end for;

exit;
