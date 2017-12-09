SetVerbose("EndoCheck", 0);

F := Rationals();
P2<x,y,z> := ProjectiveSpace(F, 2);

N := 3;
fXs := [* *]; Xs := [* *]; P0s := [* *]; Ms := [* *];

fX := (y - z)^4 - x^3*z - x*z^3 - z^4;
X := Curve(P2, fX);
P0 := X ! [0, 0, 1];
M := -Matrix(F, [ [0, 2, -2] ]);
Append(~fXs, fX); Append(~Xs, X); Append(~P0s, P0); Append(~Ms, M);

fX := (z - x)^4 - y^3*x - y*x^3 - x^4;
X := Curve(P2, fX);
P0 := X ! [1, 0, 0];
M := -Matrix(F, [ [-2, 0, 2] ]);
Append(~fXs, fX); Append(~Xs, X); Append(~P0s, P0); Append(~Ms, M);

fX := (x - y)^4 - z^3*y - z*y^3 - y^4;
X := Curve(P2, fX);
P0 := X ! [0, 1, 0];
M := Matrix(F, [ [2, -2, 0] ]);
Append(~fXs, fX); Append(~Xs, X); Append(~P0s, P0); Append(~Ms, M);

S<t> := PolynomialRing(F);
fY := t^3 + t + 1;
Y := HyperellipticCurve(fY);
Q0 := Y ! [0, 1, 1];
//Q0 := Y ! [1, 0, 0];

for i in [1..N] do
    fX := fXs[i]; X := Xs[i]; P0 := P0s[i]; M := Ms[i];

    print "";
    print "Field:";
    print F;
    print "Curve X:";
    print X;
    print "Point P0:";
    print P0;
    print "Curve Y:";
    print Y;
    print "Point Q0:";
    print Q0;
    print "Tangent representation:";
    print M;

    print "Calculating Cantor representation...";
    time test, fs := CantorFromMatrixSplit(X, P0, Y, Q0, M : LowerBound := 1);
    R<x,y> := Parent(fs[1]);
    print fs;

    R<x,y> := PolynomialRing(F, 2);
    K := FieldOfFractions(R);
    print "Original equation:", fX;
    print "First affine patch:", R ! DefiningEquation(AffinePatch(X, 1));

    /* Check that the answer is a projection after transformation, automatic: */
    fX := R ! DefiningEquation(AffinePatch(X, 1)); fY := R ! DefiningEquation(AffinePatch(Y, 1));
    IX := ideal<R | fX>;
    print "Well-defined?", R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX;

    /* Check that the action on differentials is correct: */
    fX := R ! DefiningEquation(AffinePatch(X, 1));
    fY := R ! DefiningEquation(AffinePatch(Y, 1));
    dx := K ! 1;
    dy := K ! -Derivative(fX, 1)/Derivative(fX, 2);
    ev := ((K ! Derivative(fs[1], 1))*dx + (K ! Derivative(fs[1], 2))*dy) / (K ! (2*fs[2]));
    mults := [ x, y, 1 ];
    print "Correct pullback?", R ! Numerator(K ! (ev - &+[ M[1,i]*mults[i] : i in [1..Genus(X)] ]/Derivative(fX, 2))) in IX;

    AX := AffinePatch(X, 1); AY := AffinePatch(Y, 1);
    KX := FunctionField(AX); KY := FunctionField(AY);
    m := map<AX -> AY | fs >;
    print "Degree:", Degree(ProjectiveClosure(m));
end for;

exit;
