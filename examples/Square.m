AttachSpec("../endomorphisms/magma/spec");

//load "SquareCRDataInput.m";
ListG2 := [ <-8,[50000,3750,-125]> ];

prec := 200;
CC := ComplexFieldExtra(prec);
CCSmall := ComplexFieldExtra(5);
F := RationalsExtra();

for tup in ListG2 do
    D := tup[1];
    X := HyperellipticCurveFromG2Invariants(tup[2]);
    f, h := HyperellipticPolynomials(X);
    g := (4*f + h^2)/4;
    x0 := -1;
    repeat
        x0 +:= 1;
        y0 := Evaluate(g, x0);
    until y0 ne 0;
    g *:= y0;
    X := HyperellipticCurve(g);
    PX := X ! [ x0, y0, 1 ];
    print "";
    print "Source curve:";
    print X;
    print "Point on source curve:";
    print PX;

    E := EllipticCMCurve(D);
    PE := E ! [1, 0, 0];
    print "";
    print "Target curve:";
    print E;
    print "Point on target curve:";
    print PE;

    eqsCC := EmbedCurveEquations(X, prec);
    eqsF := DefiningEquations(X);
    P := PeriodMatrix(eqsCC, eqsF : MolinNeurohr := true);

    eqsCC := EmbedCurveEquations(E, prec);
    eqsF := DefiningEquations(E);
    Q := PeriodMatrix(eqsCC, eqsF : MolinNeurohr := true);

    P := ChangeRing(P, CC); Q := ChangeRing(Q, CC);
    gen, deg := MorphismOfSmallDegreePartial(P, Q);
    print "";
    print "Degree of small morphism:", deg;
    print "Homology representation of morphism:";
    print gen[2];

    fs := RelativeMinimalPolynomials([ gen ], F);
    Kgen := RelativeSplittingField(fs);
    KX := BaseRing(X);
    KE := BaseRing(E);

    /* TODO: Write a function RelativeCompositumExtra that generalizes the
     * following step */
    L, phis := RelativeCompositum([* Kgen, KX, KE *]);
    phigen, phiX, phiE := Explode(phis);
    prec := 100;
    for iotaL in InfinitePlaces(L) do
        testX := Abs(Evaluate(phiX(KX.1), iotaL : Precision := prec) - Evaluate(KX.1, KX`iota : Precision := prec)) lt 10^(-prec + 10);
        testE := Abs(Evaluate(phiE(KE.1), iotaL : Precision := prec) - Evaluate(KE.1, KE`iota : Precision := prec)) lt 10^(-prec + 10);
        if testX and testE then
            SetInfinitePlace(L, iotaL);
        end if;
    end for;
    print "";
    print "Compositum over which verification will be run:";
    print L;
    print "Infinite place:";
    print L`iota;

    gen_alg := AlgebraizeMatrixInRelativeField(gen[1], L);
    gen := [* gen_alg, gen[2], gen[1] *];

    print "";
    print "Tangent representation of morphism:";
    print gen[1];

    L<s> := Compositum(Compositum(Kgen, KX), KE);
    print "Compositum over which verification will be run:";
    print L;

    XL := ChangeRingCurve(X, phiX);
    EL := ChangeRingCurve(E, phiE);
    PXL := [ phiX(c) : c in Eltseq(PX) ]; PXL := XL ! PXL;
    PEL := [ phiE(c) : c in Eltseq(PE) ]; PEL := EL ! PEL;
    AL := Matrix(L, [ [ phigen(c) : c in Eltseq(row) ] : row in Rows(gen[1]) ]);

    test, fs := CantorFromMatrixAmbientSplit(XL, PXL, EL, PEL, AL);
    print "";
    print "Algebraic equations for morphism:";
    print fs;
    print "Check for correctness:";
    print CorrespondenceVerifyG1(XL, PXL, EL, PEL, AL, fs);
end for;

exit;
