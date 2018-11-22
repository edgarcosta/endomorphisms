/***
 *  Cantor equation functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "Branches.m": PuiseuxLeadingExponent, InitializeImageBranch;
import "Initialize.m": InitializeCurve, ChangeTangentAction;
import "FractionalCRT.m": FractionalCRTQQ, RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;


forward CandidateFunctions;
forward FunctionValuesFromApproximations;
forward FunctionsFromApproximations;

forward CheckApproximation;
forward CheckCantor;

forward CantorFromMatrixAmbientGlobal;
forward CantorFromMatrixAmbientSplit;
forward CantorFromMatrixByDegree;

forward ChangeFunctions;
forward AbsoluteToRelative;
forward RelativeToAbsolute;


function CandidateFunctions(X, d)
/* Candidate numerators and denominators for Cantor functions */
// TODO: Use Riemann-Roch space instead

g := X`g; f := X`DEs[1]; R := X`RA;
x := R.1; y := R.2;
/* Change in hyperelliptic case for greater effectiveness: */
if Degree(f, x) lt Degree(f, y) then
    x := R.2; y := R.1;
end if;

if X`is_hyperelliptic or (g eq 1) then
    /* This case distinction is a bit silly */
    if Degree(f) mod 2 eq 1 then
        nums := [ x^i : i in [0..(d div 2)] ] cat [ x^i*y : i in [0..((d - 2*g - 1) div 2)] ];
    else
        nums := [ x^i : i in [0..(d div 2)] ] cat [ x^i*y : i in [0..((d - 2*g - 2) div 2)] ];
    end if;
    dens := [ x^i : i in [0..(d div 2)] ];
elif X`is_planar then
    nums := [ x^i*y^j : i in [0..d], j in [0..(Degree(f, y) - 1)] | i + j le d ];
    dens := nums;
end if;
return dens, nums;

end function;


function FunctionValuesFromApproximations(Y, Qs)
/* Evaluates Cantor functions in the branches Qs */

PR := Parent(Qs[1][1]); R<t> := PolynomialRing(PR);
pol_approx := &*[ t - Q[1] : Q in Qs ];
/* Start with trace and end with norm */
as_approx := Reverse(Coefficients(pol_approx)[1..Y`g]);
v := Matrix([ [ Q[2] : Q in Qs ] ]);
M := Transpose(Matrix([ [ Q[1]^i : i in [0..(Y`g - 1)] ] : Q in Qs ]));
w := v*M^(-1);
bs_approx := Reverse(Eltseq(w));
return as_approx cat bs_approx;

end function;


function FunctionsFromApproximations(X, Y, P, Qs, d)
/* Finds candidate functions in degree d from a developed branch */

I := ideal<X`RA | X`DEs[1]>;
dens, nums := CandidateFunctions(X, d);
fs_approx := FunctionValuesFromApproximations(Y, Qs);

fs := [ ];
for f_approx in fs_approx do
    ev_dens := [ -f_approx * Evaluate(den, P) : den in dens ];
    ev_nums := [ Evaluate(num, P) : num in nums ];
    evs := ev_dens cat ev_nums;
    prec := Floor(Minimum([ AbsolutePrecision(ev) : ev in evs ]));
    M := Matrix([ [ X`F ! Coefficient(ev, i) : i in [0..(prec - 1)] ] : ev in evs ]);
    Ker := Kernel(M);

    if Dimension(Ker) eq 0 then
        return false, [ ];

    else
        B := Basis(Ker);
        for b in B do
            v := Eltseq(b);
            f := &+[ v[i + #dens]*nums[i] : i in [1..#nums] ] / &+[ v[i]*dens[i] : i in [1..#dens] ];
            /* The next check should be superfluous if the precision is above a
             * small bound */
            //if (X`RA ! Numerator(X`KA ! f)) in I then
            //    return false, [ ];
            //end if;
            Append(~fs, X`KA ! f);
            break;
        end for;
    end if;
end for;
return true, fs;

end function;


function CheckApproximation(X, Y, P, Qs, fs)
/*
 * Verifies if the given functions approximate well
 */

g := Y`g;
as := fs[1..g]; bs := fs[(g + 1)..(2*g)];
for Q in Qs do
    if not IsWeaklyZero(Q[1]^g + &+[ Evaluate(as[i], P) * Q[1]^(g - i) : i in [1..g] ]) then
        return false;
    end if;
    if not IsWeaklyZero(Q[2]   - &+[ Evaluate(bs[i], P) * Q[1]^(g - i) : i in [1..g] ]) then
        return false;
    end if;
end for;
return true;

end function;


function CheckCantor(X, Y, fs)
/*
 * Verifies if the given functions satisfy the Cantor equations and hence
 * defines a morphism
 */

I := ideal<X`RA | X`DEs[1]>;
for cantor_eq in Y`cantor_eqs do
    if not (X`RA ! Numerator(X`KA ! Evaluate(cantor_eq, fs))) in I then
        return false;
    end if;
end for;
return true;

end function;


intrinsic CantorFromMatrixAmbientGlobal(X::Crv, P0:: Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^5, LowerBound := 1, UpperBound := Infinity()) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding Cantor morphism (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

InitializeCurve(X, P0); InitializeCurve(Y, Q0 : NonWP := true);
NormM := ChangeTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);

d := LowerBound;
Iterator := InitializedIterator(X, Y, NormM, 2*Y`g + 2);

while true do
    found, fs, Iterator := CantorFromMatrixByDegree(X, Y, Iterator, d : Margin := Margin);
    if found then
        return true, ChangeFunctions(X, Y, fs);
    end if;
    d +:= 1;
    if d gt UpperBound then
        return false, [];
    end if;
end while;

end intrinsic;


intrinsic CantorFromMatrixAmbientSplit(X::Crv, P0:: Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^5, LowerBound := 1, UpperBound := Infinity(), B := 300) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding Cantor morphism (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

/* We start at a suspected estimate and then increase degree until we find an appropriate divisor */
InitializeCurve(X, P0); InitializeCurve(Y, Q0 : NonWP := true);
NormM := ChangeTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);

vprintf EndoCheck, 3 : "Differential bases on factors:";
vprint EndoCheck, 3: X`NormB;
vprint EndoCheck, 3: Y`NormB;

/* Some global elements needed below */
F := X`F; OF := X`OF; RX := X`RA; KX := X`KA;
/* Bit more global margin just to be sure */

vprintf EndoCheck, 3 : "Initializing iterator...";
Iterator, f := InitializedIterator(X, Y, NormM, 2*Y`g + 2);
P := Iterator[1]; Qs := Iterator[2];
vprintf EndoCheck, 3 : "done.";

prs := [ ]; fss_red := [* *];
I := ideal<X`OF | 1>;

d := LowerBound;
while true do
    /* Find new prime */
    repeat
        pr, h := RandomSplitPrime(f, B);
    until not pr in prs;
    Append(~prs, pr); I *:= pr;
    vprintf EndoCheck : "Split prime over %o\n", #Codomain(h);

    /* Add corresponding data */
    X_red := ReduceCurveSplit(X, h); Y_red := ReduceCurveSplit(Y, h);
    NormM_red := ReduceMatrixSplit(NormM, h);

    Iterator_red := InitializedIterator(X_red, Y_red, NormM_red, 2*Y`g + 2);
    while true do
        found, fs_red, Iterator_red := CantorFromMatrixByDegree(X_red, Y_red, Iterator_red, d : Margin := Margin);
        /* If that does not work, give up and try one degree higher. Note that
         * d is initialized in the outer loop, so that we keep the degree that
         * works. */
        if found then
            break;
        end if;
        d +:= 1;
        if d gt UpperBound then
            return false, [];
        end if;
    end while;
    Append(~fss_red, fs_red);

    vprintf EndoCheck : "Fractional CRT... ";
    fs := [ ];
    for i:=1 to #fss_red[1] do
        num := RX ! 0;
        for mon in Monomials(Numerator(fss_red[1][i])) do
            exp := Exponents(mon);
            rs := [* MonomialCoefficient(Numerator(fss_red[j][i]), exp) : j in [1..#fss_red] *];
            num +:= FractionalCRTSplit(rs, prs : I := I) * Monomial(RX, exp);
        end for;
        den := RX ! 0;
        for mon in Monomials(Denominator(fss_red[1][i])) do
            exp := Exponents(mon);
            rs := [* MonomialCoefficient(Denominator(fss_red[j][i]), exp) : j in [1..#fss_red] *];
            den +:= FractionalCRTSplit(rs, prs) * Monomial(RX, exp);
        end for;
        Append(~fs, KX ! (num / den));
    end for;
    vprintf EndoCheck : "done.\n";

    vprintf EndoCheck : "Checking:\n";
    vprintf EndoCheck : "Step 1... ";
    test1 := CheckApproximation(X, Y, P, Qs, fs);
    vprintf EndoCheck : "done.\n";

    if test1 then
        vprintf EndoCheck : "Step 2... ";
        test2 := CheckCantor(X, Y, fs);
        vprintf EndoCheck : "done.\n";
        if test2 then
            vprintf EndoCheck : "Functions found!\n";
            return true, ChangeFunctions(X, Y, fs);
        end if;
    end if;
end while;

end intrinsic;


function CantorFromMatrixByDegree(X, Y, Iterator, d : Margin := 2^5)

vprintf EndoCheck, 2 : "Trying degree %o...\n", d;
dens, nums := CandidateFunctions(X, d);
n := #dens + #nums + Margin;
e := ExponentDenominator(Iterator[2][1][1]);
vprintf EndoCheck, 2 : "Number of digits in expansion: %o.\n", n*e;

/* Take non-zero image branch */
vprintf EndoCheck, 2 : "Expanding branches... ";
while true do
    P, Qs, _, _ := Explode(Iterator);
    prec := Minimum([ RelativePrecision(c) : c in P cat &cat(Qs) ]);
    if prec ge n then
        break;
    end if;
    Iterator := IterateIterator(Iterator);
end while;
P, Qs, _, _ := Explode(Iterator);
vprintf EndoCheck, 2 : "done.\n";

/* Fit a Cantor morphism to it */
vprintf EndoCheck, 2 : "Solving linear system... ";
test, fs := FunctionsFromApproximations(X, Y, P, Qs, d);
vprintf EndoCheck, 2 : "done.\n";

if test then
    vprintf EndoCheck, 2 : "Checking:\n";
    vprintf EndoCheck, 2 : "Step 1... ";
    test1 := CheckApproximation(X, Y, P, Qs, fs);
    vprintf EndoCheck, 2 : "done.\n";
    if test1 then
        vprintf EndoCheck, 2 : "Step 2... ";
        test2 := CheckCantor(X, Y, fs);
        vprintf EndoCheck, 2 : "done.\n";
        if test2 then
            vprintf EndoCheck, 2 : "Functions found!\n";
            return true, fs, Iterator;
        end if;
    end if;
end if;
return false, [], Iterator;

end function;


function ChangeFunctions(X, Y, fs)
/*
 * Change the functions on patches to rational functions on the original curves.
 * Especially relevant when considering elliptic curve factors.
 */

/* For now we only do this in genus 1: */
g := X`g; R := X`RA; K := X`KA;
if Y`g eq 1 then
subsX := [ K ! X`RA.1, K ! X`RA.2 ];
if X`is_hyperelliptic or (X`g eq 1) then
    if X`patch_index eq 3 then
        subsX := [ subsX[2] / subsX[1]^(g + 1), 1 / subsX[1] ];
    end if;
elif X`is_planar then
    if X`patch_index eq 2 then
        subsX := [ subsX[1] / subsX[2], 1 / subsX[2] ];
    elif X`patch_index eq 3 then
        subsX := [ subsX[2] / subsX[1], 1 / subsX[1] ];
    end if;
end if;
if X`unif_index eq 2 then
    subsX := [ subsX[2], subsX[1] ];
end if;
fs := [ X`KA ! Evaluate(f, subsX) : f in fs ];
end if;

// TODO: Transform to Y in general by thinking about how the Cantor equation
// changes. But this is nasty.
if Y`g eq 1 then
    fs := [ -fs[1], fs[2] ];
    if Y`unif_index eq 2 then
        fs := [ fs[2], fs[1] ];
    end if;
    if Y`patch_index eq 3 then
        fs := [ 1 / fs[2], fs[1] / fs[2]^2 ];
    end if;
end if;

/* Making denominators contains x only */
if X`is_hyperelliptic then
    if IsAffine(X) then
        A := X;
    else
        A := AffinePatch(X, 1);
    end if;
    S := PolynomialRing(X`F); T := PolynomialRing(S);
    DER := R ! DefiningEquations(A)[1]; DET := AbsoluteToRelative(DER, R, T);
    _, h := HyperellipticPolynomials(X); h := T ! S ! h; h := RelativeToAbsolute(h, R, T);

    fs_red := [ ];
    for f in fs do
        num := R ! Numerator(f); den := R ! Denominator(f);
        den_conj := Evaluate(den, [R.1, -R.2 - h ]);
        num *:= den_conj; den *:= den_conj;
        num := AbsoluteToRelative(num, R, T); den := AbsoluteToRelative(den, R, T);
        num := RelativeToAbsolute(num mod DET, R, T); den := RelativeToAbsolute(den mod DET, R, T);
        Append(~fs_red, num / den);
    end for;
    return fs_red;
end if;
return fs;

end function;


function AbsoluteToRelative(fR, R, T)
/*
 * R is a polynomial ring in two variables,
 * T the same ring but as a relative polynomial ring
 */

S := BaseRing(T);
fT := T ! 0;
for mon in Monomials(fR) do
    c := MonomialCoefficient(fR, mon); exp := Exponents(mon);
    fT +:= c * S.1^exp[1] * T.1^exp[2];
end for;
return fT;

end function;


function RelativeToAbsolute(fT, R, T)
/*
 * R is a polynomial ring in two variables,
 * T the same ring but as a relative polynomial ring
 */

S := BaseRing(T); h := hom<S -> R | [ R.1 ] >;
fR := R ! 0;
for mon in Monomials(fT) do
    c := MonomialCoefficient(fT, mon); exp := Exponents(mon);
    fR +:= h(c) * R.2^exp[1];
end for;
return fR;

end function;
