/***
 *  Recognizing complex numbers as algebraic numbers in relative fields
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


forward MinimalPolynomialLLL;
forward AlgebraizeElementLLL;


function ScalingFactor(aCC)
/* Simple rational factor such that aCC / fac has small absolute value */

return 1;
if Abs(aCC) le Parent(aCC)`epscomp then return 1; end if;

absa := Abs(aCC); fac := 1;
while absa gt 2 do fac *:= 2; absa /:= 2; end while;
while absa lt 1/2 do fac *:= 1/2; absa /:= 1/2; end while;
return fac;

absa := Abs(aCC);
if absa gt 1 then d := Floor(absa); else d := 1 / Floor(1 / absa); end if;
bCC := aCC / d; absb := Abs(bCC);
d *:= (Floor(4*absb) / 4);
return d;

absa := Abs(aCC);
if absa gt 1 then d := Floor(absa); else d := 1 / Floor(1 / absa); end if;
bCC := aCC / d; absb := Abs(bCC);
d *:= (Floor(4*absb) / 4);
return d;

end function;


function TestCloseToRoot(f, aCC);
/* Tests whether complex number aCC is close to a root of polynomial f over
 * NumberFieldExtra. One Newton should do it */

if IsZero(f) then return false; end if;
if Degree(f) eq 0 then return false; end if;
fCC := EmbedPolynomialExtra(f); csCC := Coefficients(fCC); CC := Parent(aCC);
abscsCC := Max([ Abs(c) : c in csCC ]);
absaCC := Max([1, Abs(aCC)^Degree(f) ]);
abs := abscsCC*absaCC*CC`height_bound;
//if not &+[ csCC[i + 1]*aCC^i : i in [0..Degree(f)] ] lt abs then
//    return false;
//end if;
aCCnew := aCC - (Evaluate(fCC, aCC)/Evaluate(Derivative(fCC), aCC));
if Abs(aCCnew - aCC) lt CC`epscomp then
    return true;
end if;
return false;

 /*
fCC := EmbedPolynomialExtra(f); rtsCC := [ tup[1] : tup in Roots(fCC) ];
for rtCC in rtsCC do
    abs := Abs(rtCC - aCC);
    if abs le Parent(aCC)`epscomp then
        vprint EndoFind, 3 : "";
        vprint EndoFind, 3 : "Distance to root:";
        vprint EndoFind, 3 : RealField(5) ! abs;
        return true;
    end if;
end for;
return false;
*/

end function;


function MinimalPolynomialLLL(aCC, K : UpperBound := Infinity())
/* Returns a relative minimal polynomial of the complex number aCC with respect
 * to the stored infinite place of K. */

nsavedrows := 2;
CCSmall := ComplexField(5);
vprint EndoFind, 3 : "";
vprint EndoFind, 3 : "Determining minimal polynomial using LLL for:";
vprint EndoFind, 3 : CCSmall ! aCC;
vprint EndoFind, 3 : "over:";
vprint EndoFind, 3 : K;

degK := Degree(K); R := PolynomialRing(K); genK := K.1; CCK := K`CC;
assert Precision(Parent(aCC)) ge Precision(CCK);

/* Scale and take the images of a power basis of K over QQ */
/* This is all over a complex field CC with higher precision than CCK */
faca := ScalingFactor(aCC); aCCsc := aCC / faca;
genCC := EmbedExtra(genK); CCiota := Parent(genCC);
facK := ScalingFactor(genCC); genK /:= facK; genCC /:= facK;
powersgenCC := [ EmbedExtra(genK^i) : i in [0..(degK - 1)] ];

/* Create first entry corresponding to constant term */
MLine := [ ]; poweraCCsc := 1; degf := 0;
MLine cat:= [ powergenCC * poweraCCsc : powergenCC in powersgenCC ];

/* Successively adding other entries to find relations */
savedrows := [* 0 : i in [1..nsavedrows] *];
while degf lt UpperBound do
    degf +:= 1; poweraCCsc *:= aCCsc;
    MLine cat:= [ powergenCC * poweraCCsc : powergenCC in powersgenCC ];
    M := Transpose(Matrix(CCK, [ MLine ]));
    Ker, test_ker := IntegralLeftKernel(M : CalcAlg := true);
    vprint EndoFind, 3 : "";
    vprint EndoFind, 3 : "Kernel during LLL minimal polynomial determination:";
    vprint EndoFind, 3 : Ker;

    if test_ker then
        row := Rows(Ker)[1]; seq := Eltseq(row); ht := Max([ Height(c) : c in seq ]);
        f := &+[ &+[ seq[i*degK + j + 1]*genK^j : j in [0..(degK - 1)] ] * R.1^i : i in [0..degf] ];
        f := Evaluate(f, R.1 / faca); f /:= LeadingCoefficient(f);
        if TestCloseToRoot(f, aCC) then
            if ht lt CCK`height_bound then
                vprint EndoFind, 3 : "";
                vprint EndoFind, 3 : f;
                vprint EndoFind, 3 : "done determining minimal polynomial using LLL.";
                return f;
            elif &and[ Type(entry) eq SeqEnum : entry in savedrows ] then
                testrepeat := true;
                for i in [1..nsavedrows] do
                    seqold := savedrows[i] cat [ 0 : j in [1..(nsavedrows - i + 1)*degK] ];
                    if not seq eq seqold then
                        testrepeat := false;
                        break;
                    end if;
                end for;
                if testrepeat then
                    vprint EndoFind, 3 : "";
                    vprint EndoFind, 3 : f;
                    vprint EndoFind, 3 : "done determining minimal polynomial using LLL.";
                    return f;
                end if;
            end if;

            for i in [2..nsavedrows] do
                savedrows[i - 1] := savedrows[i];
            end for;
            savedrows[nsavedrows] := seq;
        end if;
    end if;
end while;
return "Failed to find minimal polynomial using LLL";

end function;


function AlgebraizeElementLLL(aCC, K)
/* Returns an approximation of the complex number aCC as an element of K. */

genK := K.1; genCC := EmbedExtra(genK); CCiota := Parent(genCC);
facK := ScalingFactor(genCC); genK /:= facK; genCC /:= facK;
CCK := K`CC; prec := Precision(CCK);
assert Precision(Parent(aCC)) ge Precision(CCK);

degK := Degree(K);
MLine := [ EmbedExtra(genK^i) : i in [0..(degK - 1)] ] cat [ -aCC ];
M := Transpose(Matrix(CCK, [ MLine ]));

vprint EndoFind, 3 : "";
vprint EndoFind, 3 : "Algebraizing element...";
vprint EndoFind, 3 : "";

Ker, test_ker := IntegralLeftKernel(M : CalcAlg := true);
row := Rows(Ker)[1]; den := row[#Eltseq(row)];
/* There is no height test in this case: Be aware of the corresponding risk */
if den ne 0 then
    a := &+[ row[i + 1]*genK^i : i in [0..(degK - 1)] ] / den;
    if Abs(EmbedExtra(a) - aCC) le CCK`epscomp then
        vprint EndoFind, 3 : a;
        vprint EndoFind, 3 : "done algebraizing element.";
        return true, a;
    end if;
end if;
vprint EndoFind, 3 : "no element found.";
return false, 0;

end function;


intrinsic FractionalApproximation(aCC::FldComElt) -> FldRatElt
{Returns a fractional approximation of the complex number aCC.}

CC := Parent(aCC); RR := RealField(CC);
M := Matrix(RR, [ [ 1 ], [ -Real(aCC) ] ]);

Ker, test_ker := IntegralLeftKernel(M : CalcAlg := true);
if not test_ker then return Rationals() ! 0, false; end if;

q := Ker[1,1] / Ker[1,2];
if (RR ! Abs(q - aCC)) le RR`epscomp then
    return q, true;
else
    return Rationals() ! 0, false;
end if;

end intrinsic;


intrinsic FractionalApproximation(aRR::FldReElt) -> FldRatElt
{Returns a fractional approximation of the real number aRR.}

RR := Parent(aRR);
M := Matrix(RR, [ [ 1 ], [ -aRR ] ]);

K, test_ker := IntegralLeftKernel(M : CalcAlg := true);
if not test_ker then return Rationals() ! 0, false; end if;

q := K[1,1] / K[1,2];
if (RR ! Abs(q - aRR)) le RR`epscomp then
    return q, true;
else
    return Rationals() ! 0, false;
end if;

end intrinsic;


intrinsic FractionalApproximationMatrix(ACC::.) -> .
{Returns a fractional approximation of the matrix ACC.}

test := true;
rows_alg := [ ];
for row in Rows(ACC) do
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


function RationalReconstruction(r);
    CC := Parent(r); r1 := Real(r); r2 := Imaginary(r); p := Precision(r2);
    if r2 ne Parent(r2) ! 0 then e := Log(AbsoluteValue(r2)); else e := -p; end if;
    if -e lt p/2 then return false, 0; end if;
    best := 0; i := p div 10; b := BestApproximation(r1, 10^i);
    while b ne best and i le p do
        i +:= 5; best := b; b := BestApproximation(r1, 10^i);
    end while;
    if b ne best then return false, 0; end if;
    test := Abs(r - b) lt Parent(r)`epscomp;
    return test, b;
end function;


intrinsic AlgebraizeElementExtra(aCC::FldComElt, K::Fld : UseQQ := false, UseRatRec := true, minpol := 0) -> .
{Returns an approximation of the complex number aCC as an element of K. It calculates a minimal polynomial over QQ first.}

if UseRatRec and Type(K) eq FldRat then return RationalReconstruction(aCC); end if;
if not UseQQ then return AlgebraizeElementLLL(aCC, K); end if;

CCK := K`CC; CCiota := Parent(K`iota);
assert Precision(Parent(aCC)) ge Precision(CCK);

if Type(minpol) eq RngIntElt then minpol := MinimalPolynomialLLL(aCC, RationalsExtra(Precision(CCK))); end if;

vprint EndoFind, 3 : "";
vprint EndoFind, 3 : "Finding roots in field with Pari...";
rts := RootsPari(minpol, K);
vprint EndoFind, 3 : rts;
vprint EndoFind, 3 : "done finding roots in field with Pari.";

for rt in rts do if Abs(EmbedExtra(rt) - aCC) le CCK`epscomp then return true, rt; end if; end for;
return false, 0;

end intrinsic;


intrinsic AlgebraizeElementsExtra(LCC::SeqEnum, K::Fld : UseQQ := false, UseRatRec := true) -> .
{Returns an approximation of the elements of the list LCC over K.}

L := [ ];
for aCC in LCC do
    test, a := AlgebraizeElementExtra(aCC, K : UseQQ := UseQQ, UseRatRec := UseRatRec);
    if not test then return false, 0; end if;
    Append(~L, a);
end for;
return true, L;

end intrinsic;


intrinsic AlgebraizeMatrixExtra(MCC::., K::Fld : UseQQ := false, UseRatRec := true) -> .
{Returns an approximation of the complex matrix MCC over K.}

rows := [ ];
for rowCC in Rows(MCC) do
    row := [ ];
    for cCC in Eltseq(rowCC) do
        test, c := AlgebraizeElementExtra(cCC, K : UseQQ := UseQQ, UseRatRec := UseRatRec);
        if not test then return false, 0; end if;
        Append(~row, c);
    end for;
    Append(~rows, row);
end for;
return true, Matrix(rows);

end intrinsic;


intrinsic MinimalPolynomialExtra(aCC::FldComElt, K::Fld : UpperBound := Infinity(), UseQQ := true, Alg := false) -> .
{Given a complex number aCC and a NumberFieldExtra K, finds the minimal polynomial of aCC over K. The general version may be more stable than MinimalPolynomialLLL via the use of RootsPari. This minimal polynomial over QQ is the second return value.}

if Alg then
    test, a := AlgebraizeElementExtra(aCC, K);
    if test then R := PolynomialRing(K); return R.1 - a, 0; end if;
end if;
if not UseQQ then return MinimalPolynomialLLL(aCC, K : UpperBound := UpperBound), 0; end if;

CCK := K`CC; CCiota := Parent(K`iota);
assert Precision(Parent(aCC)) ge Precision(CCK);
gQQ := MinimalPolynomialLLL(aCC, RationalsExtra(Precision(CCK)) : UpperBound := UpperBound);
assert IsIrreducible(gQQ);
if Degree(gQQ) eq 1 then return gQQ, gQQ; end if;

/* First try faster algorithm for linear factors */
rts := RootsPari(gQQ, K);
R := PolynomialRing(K);
gKs := [ R.1 - rt : rt in rts ];
for gK in gKs do if TestCloseToRoot(gK, aCC) then return gK, gQQ; end if; end for;

/* Then try more general algorithm */
gKs := FactorizationPari(gQQ, K);
for gK in gKs do if TestCloseToRoot(gK, aCC) then return gK, gQQ; end if; end for;
error "Failed to find relative minimal polynomial";

end intrinsic;
