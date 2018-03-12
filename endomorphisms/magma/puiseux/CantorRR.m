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


import "Branches.m": PuiseuxRamificationIndex, InitializeImageBranch;
import "Initialize.m": InitializeCurve, ChangeTangentAction;
import "FractionalCRT.m": FractionalCRTQQ, RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "RiemannRoch.m": RRBasis, RRGenerators, RREvaluations, GlobalGenerators;


forward FunctionValuesFromApproximations;
forward VectorsFromApproximations;

forward CheckApproximation;
forward FunctionsFromVectors;
forward CheckCantor;

forward CantorFromMatrixByDegree;
forward CantorFromMatrixRRGlobal;
forward CantorFromMatrixRRSplit;

forward ChangeFunctions;
forward AbsoluteToRelative;
forward RelativeToAbsolute;


function FunctionValuesFromApproximations(Y, Qs)
/* Evaluates Cantor functions in the branches Qs */

PR := Parent(Qs[1][1]); S<t> := PolynomialRing(PR);
pol_approx := &*[ t - Q[1] : Q in Qs ];
/* Start with trace and end with norm */
as_approx := Reverse(Coefficients(pol_approx)[1..Y`g]);
v := Matrix([ [ Q[2] : Q in Qs ] ]);
M := Transpose(Matrix([ [ Q[1]^i : i in [0..(Y`g - 1)] ] : Q in Qs ]));
w := v*M^(-1);
bs_approx := Reverse(Eltseq(w));
return as_approx cat bs_approx;

end function;


function VectorsFromApproximations(X, Y, P, Qs, d)
/* Finds candidate functions in degree d from a developed branch */

I := ideal<X`RU | X`DEs[1]>;
fs_approx := FunctionValuesFromApproximations(Y, Qs);

vs := [ ];
for f_approx in fs_approx do
    ev_nums := RREvaluations(X, d, P);
    ev_dens := [ -f_approx * ev : ev in ev_nums ];
    evs := ev_dens cat ev_nums;
    e := Maximum([ Maximum([ Denominator(Valuation(c - Coefficient(c, 0))) : c in Q ]) : Q in Qs ]);
    min := Minimum([ Valuation(ev) : ev in evs ]);
    max := Minimum([ AbsolutePrecision(ev) : ev in evs ]);
    M := Matrix([ [ X`F ! Coefficient(ev, i/e) : i in [(e*min)..(e*max - 1)] ] : ev in evs ]);
    Ker := Kernel(M);

    if Dimension(Ker) eq 0 then
        return false, [ ];

    else
        nums := RRBasis(X, d);
        dens := nums;
        B := Basis(Ker);
        for b in B do
            v := Eltseq(b);
            f := &+[ v[i + #dens]*nums[i] : i in [1..#nums] ] / &+[ v[i]*dens[i] : i in [1..#dens] ];
            if f ne 0 then
                Append(~vs, v);
                break;
            else
                return false, [ ];
            end if;
        end for;
    end if;
end for;
return true, vs;

end function;


function CheckApproximation(X, Y, vs, d, P, Qs)
/*
 * Verifies if the given vectors approximate well
 */

g := Y`g; evs := RREvaluations(X, d, P);
vsa := vs[1..g]; vsb := vs[(g + 1)..(2*g)];
for Q in Qs do
    evsa := [ ];
    for i in [1..#vsa] do
        den := &+[ vsa[i][j]*evs[j] : j in [1..#evs] ];
        num := &+[ vsa[i][#evs + j]*evs[j] : j in [1..#evs] ];
        Append(~evsa, num/den);
    end for;
    if not IsWeaklyZero(Q[1]^g + &+[ evsa[i] * Q[1]^(g - i) : i in [1..g] ]) then
        return false;
    end if;
    evsb := [ ];
    for i in [1..#vsb] do
        den := &+[ vsb[i][j]*evs[j] : j in [1..#evs] ];
        num := &+[ vsb[i][#evs + j]*evs[j] : j in [1..#evs] ];
        Append(~evsb, num/den);
    end for;
    if not IsWeaklyZero(Q[2]   - &+[ evsb[i] * Q[1]^(g - i) : i in [1..g] ]) then
        return false;
    end if;
end for;
return true;

end function;


function FunctionsFromVectors(X, Y, vs, d)
/*
 * Creates the actual functions
 */

g := Y`g; B := RRBasis(X, d);
vsa := vs[1..g]; vsb := vs[(g + 1)..(2*g)];
as := [ ]; bs := [ ];
for i in [1..#vsa] do
    den := &+[ vsa[i][j]*B[j] : j in [1..#B] ];
    num := &+[ vsa[i][#B + j]*B[j] : j in [1..#B] ];
    Append(~as, num/den);
end for;
for i in [1..#vsb] do
    den := &+[ vsb[i][j]*B[j] : j in [1..#B] ];
    num := &+[ vsb[i][#B + j]*B[j] : j in [1..#B] ];
    Append(~bs, num/den);
end for;
vprint EndoCheck, 4 : "Function basis before simplification:";
vprint EndoCheck, 4 : as cat bs;
return as cat bs;

end function;


function CheckCantor(X, Y, fs)
/*
 * Verifies if the given functions satisfy the Cantor equations and hence
 * defines a morphism
 */

for cantor_eq in Y`cantor_eqs do
    if not Evaluate(cantor_eq, fs) eq 0 then
        return false;
    end if;
end for;
return true;

end function;


function CantorFromMatrixByDegree(X, Y, NormM, d : Margin := 2^4)
/* Step mod p of the above */

vprintf EndoCheck, 2 : "Trying degree %o...\n", d;
// TODO: This number n can be tweaked
n := 2*d + Margin;
e := PuiseuxRamificationIndex(NormM);
vprintf EndoCheck, 2 : "Number of digits in expansion: %o.\n", n*e;

/* Take non-zero image branch */
vprintf EndoCheck, 2 : "Expanding... ";
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, n*e);
vprintf EndoCheck, 4 : "Base point:\n";
_<t> := Parent(P[1]);
_<r> := BaseRing(Parent(P[1]));
vprint EndoCheck, 4 : P;
vprintf EndoCheck, 4 : "Resulting branches:\n";
vprint EndoCheck, 4 : Qs;
vprint EndoCheck, 4 : BaseRing(Parent(P[1]));
vprintf EndoCheck, 2 : "done.\n";

/* Fit a Cantor morphism to it */
vprintf EndoCheck, 2 : "Solving linear system... ";
test, vs := VectorsFromApproximations(X, Y, P, Qs, d);
vprintf EndoCheck, 2 : "done.\n";

if test then
    // TODO: In some circumstances this step is optional, we could skip it with
    // a flag
    vprintf EndoCheck, 2 : "Checking:\n";
    vprintf EndoCheck, 2 : "Step 1... ";
    test1 := CheckApproximation(X, Y, vs, d, P, Qs);
    vprintf EndoCheck, 2 : "done.\n";
    if test1 then
        vprintf EndoCheck, 2 : "Step 2... ";
        fs := FunctionsFromVectors(X, Y, vs, d);
        test2 := CheckCantor(X, Y, fs);
        vprintf EndoCheck, 2 : "done.\n";
        if test2 then
            vprintf EndoCheck, 2 : "Functions found!\n";
            return true, fs, vs;
        end if;
    end if;
end if;
return false, [], [];

end function;


intrinsic CantorFromMatrixRRGlobal(X::Crv, P0:: Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity()) -> Sch
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding Cantor morphism (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

InitializeCurve(X, P0); InitializeCurve(Y, Q0);
X`RRgens := RRGenerators(X);
X`globgens, X`DEs_sub := GlobalGenerators(X);

vprintf EndoCheck, 3 : "Tangent matrix before change of basis: ";
vprint EndoCheck, 3 : M;
NormM := ChangeTangentAction(X, Y, M);
vprintf EndoCheck, 3 : "Tangent matrix after change of basis: ";
vprint EndoCheck, 3 : NormM;
NormM := Y`T * NormM * (X`T)^(-1);
tjs0, f := InitializeImageBranch(NormM);

d := LowerBound;
while true do
    found, fs, vs := CantorFromMatrixByDegree(X, Y, NormM, d : Margin := 2^4);
    if found then
        return true, ChangeFunctions(X, Y, fs);
    end if;
    /* If that does not work, give up and try one degree higher: */
    d +:= 1;
    if d gt UpperBound then
        return false, [];
    end if;
end while;

end intrinsic;


intrinsic CantorFromMatrixRRSplit(X::Crv, P0:: Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), B := 300) -> Sch
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding Cantor morphism (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

InitializeCurve(X, P0); InitializeCurve(Y, Q0);
X`RRgens := RRGenerators(X);
X`globgens, X`DEs_sub := GlobalGenerators(X);

vprintf EndoCheck, 3 : "Tangent matrix before change of basis: ";
vprint EndoCheck, 3 : M;
NormM := ChangeTangentAction(X, Y, M);
vprintf EndoCheck, 3 : "Tangent matrix after change of basis: ";
vprint EndoCheck, 3 : NormM;
NormM := Y`T * NormM * (X`T)^(-1);
tjs0, f := InitializeImageBranch(NormM);

/* Some global elements needed below */
gY := Y`g; F := X`F; rF := X`rF; OF := X`OF; BOF := X`BOF;
/* TODO: Add decent margin here, + 1 already goes wrong occasionally */
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, gY + 2);

ps_rts := [ ]; prs := [ ]; vss_red := [* *];
I := ideal<X`OF | 1>;

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
        found, fs_red, vs_red := CantorFromMatrixByDegree(X_red, Y_red, NormM_red, d : Margin := Margin);
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
    Append(~vss_red, vs_red);

    vprintf EndoCheck : "Fractional CRT... ";
    vs := [ ];
    for i in [1..#vss_red[1]] do
        v_reds := [* vs_red[i] : vs_red in vss_red *];
        v := [ FractionalCRTSplit([* v_red[j] : v_red in v_reds *], prs, OF, I, BOF, BI, F) : j in [1..#v_reds[1]] ];
        Append(~vs, v);
    end for;
    vprintf EndoCheck : "done.\n";

    vprintf EndoCheck : "Checking:\n";
    vprintf EndoCheck : "Step 1... ";
    test1 := CheckApproximation(X, Y, vs, d, P, Qs);
    vprintf EndoCheck : "done.\n";

    if test1 then
        vprintf EndoCheck : "Step 2... ";
        fs := FunctionsFromVectors(X, Y, vs, d);
        test2 := CheckCantor(X, Y, fs);
        vprintf EndoCheck : "done.\n";
        if test2 then
            vprintf EndoCheck : "Functions found!\n";
            return true, ChangeFunctions(X, Y, fs);
        end if;
    end if;
end while;

end intrinsic;


function ChangeFunctions(X, Y, fs)
/*
 * Change the functions on patches to rational functions on the original curves.
 * Especially relevant when considering elliptic curve factors.
 */

/* For now we only do this in genus 1: */
g := X`g; RA := X`RA; KA := X`KA;
if Y`g eq 1 then
    subsX := [ KA ! X`RA.1, KA ! X`RA.2 ];
    if X`is_hyperelliptic or (X`g eq 1) then
        if X`patch_index eq 3 then
            vprint EndoCheck, 3 : "Modifying functions for patch index of X";
            subsX := [ subsX[2] / subsX[1]^(g + 1), 1 / subsX[1] ];
        end if;
    elif X`is_planar then
        if X`patch_index eq 2 then
            vprint EndoCheck, 3 : "Modifying functions for patch index of X";
            subsX := [ subsX[1] / subsX[2], 1 / subsX[2] ];
        elif X`patch_index eq 3 then
            vprint EndoCheck, 3 : "Modifying functions for patch index of X";
            subsX := [ subsX[2] / subsX[1], 1 / subsX[1] ];
        end if;
    end if;
    if X`unif_index eq 2 then
        vprint EndoCheck, 3 : "Modifying functions for uniformizer index of X";
        subsX := [ subsX[2], subsX[1] ];
    end if;
    fs := [ Evaluate(X`KA ! f, subsX) : f in fs ];
end if;

// TODO: Transform to Y in general by thinking about how the Cantor equation
// changes. But this is nasty.
if Y`g eq 1 then
    fs := [ -fs[1], fs[2] ];
    if Y`unif_index eq 2 then
        vprint EndoCheck, 3 : "Modifying functions for uniformizer index of Y";
        fs := [ fs[2], fs[1] ];
    end if;
    if Y`patch_index eq 3 then
        vprint EndoCheck, 3 : "Modifying functions for patch index of Y";
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
    DER := RA ! DefiningEquations(A)[1]; DET := AbsoluteToRelative(DER, RA, T);
    _, h := HyperellipticPolynomials(X); h := T ! S ! h; h := RelativeToAbsolute(h, RA, T);

    fs_red := [ ];
    for f in fs do
        num := RA ! Numerator(f); den := RA ! Denominator(f);
        den_conj := Evaluate(den, [RA.1, -RA.2 - h ]);
        num *:= den_conj; den *:= den_conj;
        num := AbsoluteToRelative(num, RA, T); den := AbsoluteToRelative(den, RA, T);
        num := RelativeToAbsolute(num mod DET, RA, T); den := RelativeToAbsolute(den mod DET, RA, T);
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
