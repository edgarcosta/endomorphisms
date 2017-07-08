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


forward CantorEquations;

forward CandidateFunctions;
forward FunctionValuesFromApproximations;
forward FunctionsFromApproximations;

forward CheckApproximation;
forward CheckCantor;

forward CantorFromMatrix;
forward CantorFromMatrixSplit;
forward CantorFromMatrixSplitStepModP;

forward ChangeFunctions;
forward AbsoluteToRelative;
forward RelativeToAbsolute;


import "Divisor.m": InitializeCurve, ChangeTangentAction;
import "LocalInfo.m": PuiseuxRamificationIndex, InitializeImageBranch;
import "FractionalCRT.m": FractionalCRTQQ, RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;


function CantorEquations(X);
/* Gives the equations in the description a (x) = 0, y = b (x) */

g := X`g; f := X`DEs[1]; R := X`R;
F := BaseRing(R); S := PolynomialRing(F, 2*g); T<t> := PolynomialRing(S);
/* Names:
 * a1 is trace term before t^(g - 1), a_g is norm term before t^0,
 * b1 is term before t^(g - 1), bg is term before t^0
 */
varnames := [ Sprintf("a%o", i) : i in [1..g] ] cat [ Sprintf("b%o", i) : i in [1..g] ];
AssignNames(~S, varnames);

/* Start with trace and end with norm */
canpol := t^g + &+[ S.i * t^(g - i) : i in [1..g] ];
substpol := &+[ S.(g + i) * t^(g - i) : i in [1..g] ];
P := [t, substpol];
eqpol := Evaluate(f, P) mod canpol;
return Coefficients(eqpol);

end function;


function CandidateFunctions(X, d)
/* Candidate numerators and denominators for Cantor functions */

f := X`DEs[1]; R := X`R; x := X`x; y := X`y;
if X`is_hyperelliptic or X`g eq 1 then
    nums := [ x^i : i in [0..(d div 2)] ] cat [ x^i*y : i in [0..((d - X`g - 1) div 2)] ];
    dens := [ x^i : i in [0..(d div 2)] ];
    dens := nums;
elif X`is_planar then
    // TODO: Use Riemann-Roch space instead
    nums := [ x^i*y^j : i in [0..d], j in [0..(Degree(f, y) - 1)] | i + j le d ];
    dens := [ x^i : i in [0..d] ];
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

I := ideal<X`R | X`DEs[1]>;
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
            //if (X`R ! Numerator(X`K ! f)) in I then
            //    return false, [ ];
            //end if;
            Append(~fs, X`K ! f);
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

I := ideal<X`R | X`DEs[1]>;
for cantor_eq in Y`cantor_equations do
    if not (X`R ! Numerator(X`K ! Evaluate(cantor_eq, fs))) in I then
        return false;
    end if;
end for;
return true;

end function;


intrinsic CantorFromMatrix(X::Crv, P0:: Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity()) -> Sch
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding Cantor morphism (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

output := InitializeCurve(X, P0); output := InitializeCurve(Y, Q0);
NormM := ChangeTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);
e := PuiseuxRamificationIndex(NormM);

d := LowerBound;
while true do
    vprintf EndoCheck : "Trying degree %o...\n", d;
    dens, nums := CandidateFunctions(X, d);
    n := #dens + #nums + Margin;
    vprintf EndoCheck : "Number of digits in expansion: %o.\n", n;

    /* Take non-zero image branch */
    vprintf EndoCheck : "Expanding... ";
    P, Qs := ApproximationsFromTangentAction(X, Y, NormM, n*e);
    vprint EndoCheck, 3 : P, Qs;
    vprintf EndoCheck : "done.\n";

    /* Fit a Cantor morphism to it */
    vprintf EndoCheck : "Solving linear system... ";
    test, fs := FunctionsFromApproximations(X, Y, P, Qs, d);
    vprintf EndoCheck : "done.\n";

    if test then
        vprintf EndoCheck : "Checking:\n";
        vprintf EndoCheck : "Step 1... ";
        test1 := CheckApproximation(X, Y, P, Qs, fs);
        vprintf EndoCheck : "done.\n";
        if test1 then
            vprintf EndoCheck : "Step 2...\n";
            test2 := CheckCantor(X, Y, fs);
            vprintf EndoCheck : "done.\n";
            if test2 then
                vprintf EndoCheck : "Functions found!\n";
                return true, ChangeFunctions(X, Y, fs);
            end if;
        end if;
    end if;

    /* If that does not work, give up and try one degree higher: */
    d +:= 1;
    if d gt UpperBound then
        return false, "";
    end if;
end while;

end intrinsic;


intrinsic CantorFromMatrixSplit(X::Crv, P0:: Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), B := 300) -> Sch
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding Cantor morphism (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

output := InitializeCurve(X, P0); output := InitializeCurve(Y, Q0);
NormM := ChangeTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);
tjs0, f := InitializeImageBranch(NormM);
e := PuiseuxRamificationIndex(NormM);

/* Some global elements needed below */
gY := Y`g; F := X`F; rF := X`rF; OF := X`OF; BOF := X`BOF; RX := X`R; KX := X`K;
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, gY + 1);

ps_rts := [ ]; prs := [ ]; fss_red := [* *];
I := ideal<X`OF | 1>;
have_to_check := true;

d := LowerBound;
while true do
    /* Find new prime */
    repeat
        p_rt := RandomSplitPrime(f, B);
        p, rt := Explode(p_rt);
    until not p in [ tup[1] : tup in ps_rts ];
    Append(~ps_rts, p_rt);
    vprintf EndoCheck : "Split prime over %o\n", p;

    /* Add corresponding data */
    pr := ideal<OF | [ p, rF - rt ]>;
    Append(~prs, pr); I *:= pr;
    X_red := ReduceCurveSplit(X, p, rt); Y_red := ReduceCurveSplit(Y, p, rt);
    NormM_red := ReduceMatrixSplit(NormM, p, rt);
    BI := Basis(I);

    while true do
        found, fs_red := CantorFromMatrixSplitStepModP(X_red, Y_red, NormM_red, d : Margin := Margin, have_to_check := have_to_check);
        /* If that does not work, give up and try one degree higher. Note that
         * d is initialized in the outer loop, so that we keep the degree that
         * works. */
        if found then
            break;
        end if;
        d +:= 1;
        if d gt UpperBound then
            return false, "";
        end if;
    end while;
    have_to_check := false;
    Append(~fss_red, fs_red);

    vprintf EndoCheck : "Fractional CRT... ";
    fs := [ ];
    for i:=1 to #fss_red[1] do
        num := RX ! 0;
        for mon in Monomials(Numerator(fss_red[1][i])) do
            exp := Exponents(mon);
            rs := [* MonomialCoefficient(Numerator(fss_red[j][i]), exp) : j in [1..#fss_red] *];
            num +:= FractionalCRTSplit(rs, prs, OF, I, BOF, BI, F) * Monomial(RX, exp);
        end for;
        den := RX ! 0;
        for mon in Monomials(Denominator(fss_red[1][i])) do
            exp := Exponents(mon);
            rs := [* MonomialCoefficient(Denominator(fss_red[j][i]), exp) : j in [1..#fss_red] *];
            den +:= FractionalCRTSplit(rs, prs, OF, I, BOF, BI, F) * Monomial(RX, exp);
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


function CantorFromMatrixSplitStepModP(X_red, Y_red, NormM_red, d : Margin := 2^4, have_to_check := true)
/* Step mod p of the above */

vprintf EndoCheck : "Trying degree %o...\n", d;
dens_red, nums_red := CandidateFunctions(X_red, d);
n := #dens_red + #nums_red + Margin;
e := PuiseuxRamificationIndex(NormM_red);
vprintf EndoCheck : "Number of digits in expansion: %o.\n", n*e;

/* Take non-zero image branch */
vprintf EndoCheck, 2 : "Expanding... ";
P_red, Qs_red := ApproximationsFromTangentAction(X_red, Y_red, NormM_red, n*e);
vprint EndoCheck, 3 : P_red, Qs_red;
vprintf EndoCheck, 2 : "done.\n";

/* Fit a Cantor morphism to it */
vprintf EndoCheck, 2 : "Solving linear system... ";
test_red, fs_red := FunctionsFromApproximations(X_red, Y_red, P_red, Qs_red, d);
vprintf EndoCheck, 2 : "done.\n";

if test_red then
    vprintf EndoCheck, 2 : "Checking:\n";
    vprintf EndoCheck, 2 : "Step 1: ";
    if not have_to_check then
        vprintf EndoCheck, 2 : "done.\n";
        vprintf EndoCheck, 2 : "Functions found!\n";
        return true, fs_red;
    end if;
    test := CheckCantor(X_red, Y_red, fs_red);
    vprintf EndoCheck, 2 : "done.\n";
    if test then
        vprintf EndoCheck, 2 : "Functions found!\n";
        have_to_check := false;
        return true, fs_red;
    end if;
end if;
return false, [ ];

end function;


function ChangeFunctions(X, Y, fs)
/*
 * Change the functions on patches to rational functions on the original curves.
 * Especially relevant when considering elliptic curve factors.
 */

g := X`g; R := X`R; K := X`K;
subsX := [ K ! X`x, K ! X`y ];
if X`is_hyperelliptic then
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
fs := [ X`K ! Evaluate(f, subsX) : f in fs ];

// TODO: Transform to Y in general by thinking about how the Cantor equation
// changes.
// If Y has genus 1, then we add a minus for passage from Cantor to naive
if Y`g eq 1 then
    if Y`patch_index eq 1 then
        fs := [ -fs[1], fs[2] ];
    elif Y`patch_index eq 3 then
        fs := [ -1/fs[1], fs[2]/fs[1]^2 ];
    end if;
end if;

if X`is_hyperelliptic then
    if IsAffine(X) then
        A := X;
    else
        A := AffinePatch(X, 1);
    end if;
    S := PolynomialRing(X`F); T := PolynomialRing(S);
    DER := R ! DefiningEquations(A)[1]; DET := AbsoluteToRelative(DER, R, T);
    _, h := HyperellipticPolynomials(X); h := T ! S ! h; h := RelativeToAbsolute(h, R, T);

    // Making denominators contains x only
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
