AttachSpec("../endomorphisms/magma/spec");

//load "SquareCRDataInput.m";
ListG2 := [ <-8,[50000,3750,-125]> ];

prec := 500;
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
    gen, deg := MorphismOfSmallDegree(P, Q, F);
    print "";
    print "Degree of small morphism:", deg;
    print "Tangent representation:";
    print gen;

    Kgen := BaseRing(gen[1]);
    KX := BaseRing(X);
    KE := BaseRing(E);
    /* This should be all right since all extensions involved are Galois */
    L<s> := Compositum(Compositum(Kgen, KX), KE);
    print "Compositum over which verification will be run:";
    print L;

    XL := ChangeRing(X, L);
    EL := ChangeRing(E, L);
    PXL := [ L ! c : c in Eltseq(PX) ]; PXL := XL ! PXL;
    PEL := [ L ! c : c in Eltseq(PE) ]; PEL := EL ! PEL;
    AL := Matrix(L, [ [ L ! c : c in Eltseq(row) ] : row in Rows(gen[1]) ]);

    test, fs := CantorFromMatrixAmbientSplit(XL, PXL, EL, PEL, AL);
    print "";
    print "Algebraic equations for morphism:";
    print fs;
    print "Check for correctness:";
    print CorrespondenceVerifyG1(XL, PXL, EL, PEL, AL, fs);
end for;

exit;
