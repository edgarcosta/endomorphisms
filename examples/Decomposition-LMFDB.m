SetVerbose("EndoFind", 0);
SetVerbose("CurveRec", 0);
SetVerbose("EndoCheck", 0);

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
/* The next line takes very long because our splitting function is too naive */
//R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![-1], R![1, 0, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, -1, 1, 0, -1, 1], R![]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![1, 0, 0, 1], R![1, 0, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![-7, 0, 0, 0, 1], R![0, 1, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![0, 0, 1, 0, 1], R![1, 0, 0, 1]); Append(~Xs, X);
R<x> := PolynomialRing(Rationals()); X := HyperellipticCurve(R![1, 1, 2, 1, 1], R![1, 1, 1, 1]); Append(~Xs, X);

Xs := [ Xs[1] ];
for X in Xs do
    print "";
    print "";
    print "Curve:";
    print X;

    facs := HeuristicJacobianFactors(X);
    for fac in facs do
        print "";
        print "Heuristic Jacobian factor:";
        Y, mor := Explode(fac);
        print Y;

        print "";
        print "Representation of morphism:";
        print mor;

        test, fs := Correspondence(X, Y, mor : CheckDegree := true);
        R<x,y> := Parent(fs[1]);
        print "";
        print "Equations of morphism:";
        print fs;
    end for;
end for;
