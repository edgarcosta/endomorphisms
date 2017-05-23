/***
 *  Recognizing complex numbers as algebraic numbers in relative fields
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic RelativeMinimalPolynomial(a::FldComElt, F::Fld) -> RngUPolElt
{Determines a relative minimal polynomial of the element a with respect to the stored infinite place of F.}

// NOTE: This algorithm will always terminate, but the result may be nonsense
// if the element is rather transcendent. There should be a better way to
// terminate but for now I do not see how.

degF := Degree(F); R<x> := PolynomialRing(F);
CC := Parent(a); RR := RealField(CC); prec := Precision(CC);
RCC := PolynomialRing(CC);

// NOTE: Here height is a parameter to play with.
degf := 0; height := 1; height0 := 100;
gen := CC ! Evaluate(F.1, F`iota : Precision := prec);
powersgen := [ gen^i : i in [0..(degF - 1)] ];

// Create first entry corresponding to constant term
powera := CC ! 1;
MLine := [ ];
MLine cat:= [ powergen * powera : powergen in powersgen ];

// Successively adding other entries to find relations
while true do
    // Increase height and number of possible relations
    degf +:= 1; height *:= height0;
    powera *:= a;
    MLine cat:= [ powergen * powera : powergen in powersgen ];
    M := Transpose(Matrix(CC, [ MLine ]));

    // Now split and take an IntegralLeftKernel
    MSplit := SplitMatrix(M);
    Ker := IntegralLeftKernel(MSplit);
    for row in Rows(Ker) do
        // Height condition
        if &and[ Abs(c) le height : c in Eltseq(row) ] then
            f := &+[ &+[ row[i*degF + j + 1]*F.1^j : j in [0..(degF - 1)] ] * x^i : i in [0..degf] ];
            // Factor (to eliminate redundancy) and check
            Fac := Factorization(f);
            for tup in Fac do
                g := tup[1];
                gCC := EmbedAtInfinitePlace(g, RCC);
                if RR ! Abs(Evaluate(gCC, a)) lt RR`epscomp then
                    return g;
                end if;
            end for;
        end if;
    end for;
end while;

end intrinsic;


intrinsic RelativeMinimalPolynomialsPartial(gens::SeqEnum, F::Fld) -> SeqEnum
{Polynomializes matrices.}

pols := [ ];
for gen in gens do
    pols_new := [ RelativeMinimalPolynomial(c, F) : c in Eltseq(gen[1]) ];
    pols cat:= pols_new;
end for;
return pols;

end intrinsic;


intrinsic FractionalApproximation(a::FldComElt) -> FldRatElt
{Fractional approximation of a complex number a.}

CC := Parent(a); RR := RealField(CC);
M := Matrix(RR, [ [ 1 ], [ -Real(a) ] ]);
Ker := IntegralLeftKernel(M); q := Ker[1,1] / Ker[1,2];
if (RR ! Abs(q - a)) lt RR`epscomp then
    return q;
else
    error Error("LLL not does return a sufficiently good fraction");
end if;

end intrinsic;


intrinsic FractionalApproximation(a::FldReElt) -> FldRatElt
{Fractional approximation of a real number a.}

RR := Parent(a);
M := Matrix(RR, [ [ 1 ], [ -a ] ]);
K := IntegralLeftKernel(M); q := K[1,1] / K[1,2];
if (RR ! Abs(q - a)) lt RR`epscomp then
    return q;
else
    error Error("LLL not does return a sufficiently good fraction");
end if;

end intrinsic;


intrinsic AlgebraizeElementInRelativeField(a::FldComElt, K::Fld) -> .
{Finds an algebraic approximation of a as an element of K.}
// TODO: This assumes that the extension is at most double.

degK := Degree(K); R<x> := PolynomialRing(K);
F := BaseField(K); degF := Degree(F);
CC := Parent(a); RR := RealField(CC); prec := Precision(CC);

genK := CC ! Evaluate(K.1, K`iota : Precision := prec); genF := CC ! Evaluate(F.1, F`iota : Precision := prec);
powersgenK := [ genK^i : i in [0..(degK - 1)] ]; powersgenF := [ genF^i : i in [0..(degF - 1)] ];
MLine := &cat[ [ powergenF * powergenK : powergenF in powersgenF ] : powergenK in powersgenK ] cat [-a];
M := Transpose(Matrix(CC, [ MLine ]));

// Now split and take an IntegralLeftKernel
MSplit := SplitMatrix(M);
Ker := IntegralLeftKernel(MSplit);
for row in Rows(Ker) do
    den := row[#Eltseq(row)];
    if den ne 0 then
        sCC := &+[ &+[ row[i*degF + j + 1]*genF^j : j in [0..(degF - 1)] ] * genK^i : i in [0..(degK - 1)] ] / den;
        // Check correct to given precision
        if (RR ! Abs(sCC - a)) lt RR`epscomp then
            s := &+[ &+[ row[i*degF + j + 1]*F.1^j : j in [0..(degF - 1)] ] * K.1^i : i in [0..(degK - 1)] ] / den;
            return s;
        end if;
    end if;
end for;
error Error("LLL fails to algebraize element in ambient");

end intrinsic;


intrinsic AlgebraizeMatrixInRelativeField(A::., K::Fld) -> AlgMatElt
{Algebraizes a matrix.}

return Matrix([ [ AlgebraizeElementInRelativeField(c, K) : c in Eltseq(row) ] : row in Rows(A) ]);

end intrinsic;
