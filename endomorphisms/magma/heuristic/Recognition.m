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


function TestCloseToRoot(f, a)
/* Tests whether a is close to a root of f */

CC := Parent(a); RCC := PolynomialRing(CC);
fCC := EmbedAtInfinitePlace(f, RCC); tups := Roots(fCC);
for tup in tups do
    rt := tup[1];
    //vprint EndoFind : "Absolute value:", RealField(5) ! Abs(a - rt);
    //vprint EndoFind : "Evaluation:", ComplexField(5) ! Evaluate(fCC, rt);
    if Abs(a - rt) le CC`epscomp then
        return true;
    end if;
end for;
return false;

end function;


intrinsic RelativeMinimalPolynomial(a::FldComElt, F::Fld : Factor := false, UpperBound := Infinity()) -> RngUPolElt
{Returns a relative minimal polynomial of the complex number a with respect to
the stored infinite place of F. If Factor is set to true (it is false by
default), then factorizations are taken in some intermediate steps.}

degF := Degree(F); R<x> := PolynomialRing(F);
CC := Parent(a); RR := RealField(CC); prec := Precision(CC);

/* The next line is slightly inefficient, but it is not a bottleneck. It takes
 * the images in CC of a power basis of F over QQ. */
powersgen := [ CC ! Evaluate(F.1^i, F`iota : Precision := prec) : i in [0..(degF - 1)] ];

/* Create first entry corresponding to constant term */
MLine := [ ];
degf := 0;
powera := CC ! 1;
MLine cat:= [ powergen * powera : powergen in powersgen ];

/* Successively adding other entries to find relations */
while degf lt UpperBound do
    degf +:= 1;
    vprint EndoFind : "Trying degree", degf;
    powera *:= a;
    MLine cat:= [ powergen * powera : powergen in powersgen ];
    M := Transpose(Matrix(CC, [ MLine ]));

    /* Split and take an IntegralLeftKernel */
    MSplit := HorizontalSplitMatrix(M);
    Ker, test_ker := IntegralLeftKernel(MSplit);
    /* NOTE: We only consider the first element for now */
    if test_ker then
        for row in [ Rows(Ker)[1] ] do
            vprint EndoFind : "Row:", row;
            test_height := &and[ Abs(c) le CC`height_bound : c in Eltseq(row) ];
            /* TODO: Uncomment if desired */
            //if true then
            if test_height then
                f := &+[ &+[ row[i*degF + j + 1]*F.1^j : j in [0..(degF - 1)] ] * x^i : i in [0..degf] ];
                if not Factor then
                    if TestCloseToRoot(f, a) then
                        return f, true;
                    end if;
                else
                    Fac := Factorization(f);
                    for tup in Fac do
                        g := tup[1];
                        if TestCloseToRoot(g, a) then
                            return g, true;
                        end if;
                    end for;
                end if;
            end if;
        end for;
    end if;
end while;
return 0, false;

end intrinsic;


intrinsic RelativeMinimalPolynomialsMatrices(gens::SeqEnum, F::Fld : UpperBound := Infinity()) -> SeqEnum
{Returns a relative minimal polynomials of the entries of the matrices gens
with respect to the stored infinite place of F.}

test := true;
pols := [ ];
for gen in gens do
    pols_new := [ ];
    for c in Eltseq(gen[1]) do
        pol_new, test_pol_new := RelativeMinimalPolynomial(c, F : UpperBound := UpperBound);
        test and:= test_pol_new;
        Append(~pols_new, pol_new);
    end for;
    pols cat:= pols_new;
end for;
return pols, test;

end intrinsic;


intrinsic FractionalApproximation(a::FldComElt) -> FldRatElt
{Returns a fractional approximation of the complex number a.}

CC := Parent(a); RR := RealField(CC);
M := Matrix(RR, [ [ 1 ], [ -Real(a) ] ]);
Ker, test_ker := IntegralLeftKernel(M);
if not test_ker then
    return Rationals() ! 0, false;
end if;
q := Ker[1,1] / Ker[1,2];
if (RR ! Abs(q - a)) lt RR`epscomp then
    return q, true;
else
    return Rationals() ! 0, false;
end if;

end intrinsic;


intrinsic FractionalApproximation(a::FldReElt) -> FldRatElt
{Returns a fractional approximation of the real number a.}
/* TODO: Copies previous function, which is stupid */

RR := Parent(a);
M := Matrix(RR, [ [ 1 ], [ -a ] ]);
K, test_ker := IntegralLeftKernel(M);
if not test_ker then
    return Rationals() ! 0, false;
end if;
q := K[1,1] / K[1,2];
if (RR ! Abs(q - a)) lt RR`epscomp then
    return q, true;
else
    return Rationals() ! 0, false;
end if;

end intrinsic;


intrinsic FractionalApproximationMatrix(A::.) -> .
{Returns a fractional approximation of the matrix A.}

test := true;
rows_alg := [ ];
for row in Rows(A) do
    row_alg := [ ];
    for c in Eltseq(row) do
        q, test_q := FractionalApproximation(c);
        test and:= test_q;
        Append(~row_alg, q);
    end for;
    Append(~rows_alg, row_alg);
end for;
return Matrix(rows_alg), test;

end intrinsic;


intrinsic AlgebraizeElementInRelativeField(a::FldComElt, K::Fld) -> .
{Returns an approximation of the complex number a as an element of K. This
assumes that K is at most a double extension of QQ.}

/* An alternative, more stable way would be to find the corresponding roots in
 * K and check which one does the trick. This would be much slower, however. */

degK := Degree(K); R<x> := PolynomialRing(K);
F := BaseField(K); degF := Degree(F);
CC := Parent(a); RR := RealField(CC); prec := Precision(CC);

genK := CC ! Evaluate(K.1, K`iota : Precision := prec); genF := CC ! Evaluate(F.1, F`iota : Precision := prec);
powersgenK := [ CC ! Evaluate(K.1^i, K`iota : Precision := prec) : i in [0..(degK - 1)] ];
powersgenF := [ CC ! Evaluate(F.1^i, F`iota : Precision := prec) : i in [0..(degF - 1)] ];
MLine := &cat[ [ powergenF * powergenK : powergenF in powersgenF ] : powergenK in powersgenK ] cat [-a];
M := Transpose(Matrix(CC, [ MLine ]));

/* Split and take an IntegralLeftKernel */
MSplit := HorizontalSplitMatrix(M);
Ker, test_ker := IntegralLeftKernel(MSplit);
/* NOTE: We only consider the first element for now */
if test_ker then
    for row in [ Rows(Ker)[1] ] do
        test_height := &and[ Abs(c) le CC`height_bound : c in Eltseq(row) ];
        /* Do not use height test for now */
        if true then
            den := row[#Eltseq(row)];
            if den ne 0 then
                sCC := &+[ &+[ row[i*degF + j + 1]*genF^j : j in [0..(degF - 1)] ] * genK^i : i in [0..(degK - 1)] ] / den;
                if (RR ! Abs(sCC - a)) lt RR`epscomp then
                    s := &+[ &+[ row[i*degF + j + 1]*F.1^j : j in [0..(degF - 1)] ] * K.1^i : i in [0..(degK - 1)] ] / den;
                    return s, true;
                end if;
            end if;
        end if;
    end for;
end if;
return 0, false;

end intrinsic;


intrinsic AlgebraizeMatrixInRelativeField(A::., K::Fld) -> AlgMatElt
{Returns approximations of the entries of A as elements of K. This assumes that
K is at most a double extension of QQ.}

test := true;
rows_alg := [ ];
for row in Rows(A) do
    row_alg := [ ];
    for c in Eltseq(row) do
        alpha, test_alpha := AlgebraizeElementInRelativeField(c, K);
        test and:= test_alpha;
        Append(~row_alg, alpha);
    end for;
    Append(~rows_alg, row_alg);
end for;
return Matrix(rows_alg), test;

end intrinsic;
