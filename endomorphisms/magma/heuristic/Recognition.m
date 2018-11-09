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
/* Tests whether complex number a is close to a root of polynomial f over
 * NumberFieldExtra */

CC := Parent(a); RCC := PolynomialRing(CC);
fCC := RCC ! EmbedAtInfinitePlacePolynomial(f); rts := [ tup[1] : tup in Roots(fCC) ];
for rt in rts do
    if Abs(a - rt) le CC`epscomp then
        return true;
    end if;
end for;
return false;

end function;


intrinsic MinimalPolynomialLLL(a::FldComElt, F::Fld : LowerBound := 1, UpperBound := Infinity(), StepSize := 1) -> RngUPolElt
{Returns a relative minimal polynomial of the complex number a with respect to
the stored infinite place of F.}

degF := Degree(F); R<x> := PolynomialRing(F); CC := F`CC;

/* The next line is slightly inefficient, but it is not a bottleneck. It takes
 * the images in CC of a power basis of F over QQ. */
powersgen := [ CC ! EvaluateExtra(F.1^i, F`iota) : i in [0..(degF - 1)] ];

/* Create first entry corresponding to constant term */
MLine := [ ];
degf := 0;
powera := CC ! 1;
MLine cat:= [ powergen * powera : powergen in powersgen ];

vprint EndoFind : "";
/* Successively adding other entries to find relations */
while degf lt UpperBound do
    degf +:= 1;
    powera *:= a;
    MLine cat:= [ powergen * powera : powergen in powersgen ];
    M := Transpose(Matrix(CC, [ MLine ]));

    if (degf ge LowerBound) and (degf mod StepSize eq 0) then
        /* Split and take an IntegralLeftKernel */
        MSplit := HorizontalSplitMatrix(M);
        Ker, test_ker := IntegralLeftKernel(MSplit : OneRow := true);
        /* We only consider the first row */
        if test_ker then
            row := Rows(Ker)[1];
            f := &+[ &+[ row[i*degF + j + 1]*F.1^j : j in [0..(degF - 1)] ] * x^i : i in [0..degf] ];
            if TestCloseToRoot(f, a) then
                return f;
            end if;
        end if;
    end if;
end while;
return "Failed to find minimal polynomial using LLL";

end intrinsic;


intrinsic FractionalApproximation(a::FldComElt) -> FldRatElt
{Returns a fractional approximation of the complex number a.}

CC := Parent(a); RR := RealField(CC);
M := Matrix(RR, [ [ 1 ], [ -Real(a) ] ]);
Ker, test_ker := IntegralLeftKernel(M : OneRow := true);
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

RR := Parent(a);
M := Matrix(RR, [ [ 1 ], [ -a ] ]);
K, test_ker := IntegralLeftKernel(M : OneRow := true);
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


intrinsic AlgebraizeElement(a::FldComElt, K::Fld : minpol := 0) -> .
{Returns an approximation of the complex number a as an element of K.}

CC := K`CC;
assert Precision(Parent(a)) ge Precision(CC);

if Type(minpol) eq RngIntElt then
    minpol := MinimalPolynomialLLL(a, RationalsExtra(Precision(CC)));
end if;
rts := RootsPari(minpol, K);

for rt in rts do
    rtCC := EvaluateExtra(rt, K`iota);
    if Abs(rtCC - CC ! a) le CC`epscomp then
        vprint EndoFind : "Root:", rt;
        return true, rt;
    end if;
end for;
return false, 0;

end intrinsic;


intrinsic AlgebraizeMatrix(M::., K::Fld) -> .
{Returns an approximation of the complex matrix M over K.}

rows := [ ];
for rowCC in Rows(M) do
    row := [ ];
    for cCC in Eltseq(rowCC) do
        test, c := AlgebraizeElement(cCC, K);
        if not test then
            return false, 0;
        end if;
        Append(~row, c);
    end for;
    Append(~rows, row);
end for;
return true, Matrix(rows);

end intrinsic;


intrinsic AlgebraizeElementLLL(a::FldComElt, K::Fld : minpol := 0) -> .
{Returns an approximation of the complex number a as an element of K.}

degK := Degree(K); R<x> := PolynomialRing(K);
F := BaseField(K); degF := Degree(F);
CC := Parent(a); RR := RealField(CC); prec := Precision(CC);

genK := CC ! EvaluateExtra(K.1, K`iota); genF := CC ! EvaluateExtra(F.1, F`iota);
powersgenK := [ CC ! EvaluateExtra(K.1^i, K`iota) : i in [0..(degK - 1)] ];
powersgenF := [ CC ! EvaluateExtra(F.1^i, F`iota) : i in [0..(degF - 1)] ];
MLine := &cat[ [ powergenF * powergenK : powergenF in powersgenF ] : powergenK in powersgenK ] cat [-a];
M := Transpose(Matrix(CC, [ MLine ]));

/* Split and take an IntegralLeftKernel */
MSplit := HorizontalSplitMatrix(M);
Ker, test_ker := IntegralLeftKernel(MSplit : OneRow := true);
/* NOTE: We only consider the first element for now */
if test_ker then
    for row in [ Rows(Ker)[1] ] do
        vprint EndoFind : "Row:", row;
        test_height := &and[ Abs(c) le CC`height_bound : c in Eltseq(row) ];
        if test_height then
            den := row[#Eltseq(row)];
            if den ne 0 then
                sCC := &+[ &+[ row[i*degF + j + 1]*genF^j : j in [0..(degF - 1)] ] * genK^i : i in [0..(degK - 1)] ] / den;
                if (RR ! Abs(sCC - a)) lt RR`epscomp then
                    s := &+[ &+[ row[i*degF + j + 1]*F.1^j : j in [0..(degF - 1)] ] * K.1^i : i in [0..(degK - 1)] ] / den;
                    return true, s;
                end if;
            end if;
        end if;
    end for;
end if;
return false, 0;

end intrinsic;


intrinsic MinimalPolynomialExtra(aCC::FldComElt, K::Fld : UpperBound := Infinity(), minpolQQ := 0) -> RngUPolElt
{Given a complex number aCC and a NumberFieldExtra K, finds the minimal polynomial of aCC over K.}

/* Use minimal polynomial over QQ */
CC := Parent(aCC); RCC := PolynomialRing(CC);
/* TODO: Deal with the horrible side effects that these declarations have */
//QQ := RationalsExtra(Precision(CC));
QQ := Rationals();
QQ`base := Rationals(); QQ`base_gen := QQ`base ! 1;
QQ`CC := K`CC; QQ`iota := QQ`CC ! 1;
if Type(minpolQQ) eq RngIntElt then
    f := MinimalPolynomialLLL(aCC, QQ : UpperBound := UpperBound);
else
    f := minpolQQ;
end if;

/* If this applies, then it is fast, so we consider this case first */
rts := RootsPari(f, K);
for rt in rts do
    rtCC := CC ! EvaluateExtra(rt, K`iota);
    if Abs(rtCC - aCC) le CC`epscomp then
        R := PolynomialRing(K);
        return R.1 - rt, f;
    end if;
end for;

gs := FactorizationPari(f, K);
for g in gs do
    gCC := EmbedAtInfinitePlacePolynomial(g);
    for tuprt in Roots(gCC) do
        rtCC := tuprt[1];
        //if Abs(rtCC - aCC) le 10^50*CC`epscomp then
        if Abs(rtCC - aCC) le CC`epscomp then
            return g, f;
        end if;
    end for;
end for;
error "Failed to find relative minimal polynomial";

end intrinsic;
