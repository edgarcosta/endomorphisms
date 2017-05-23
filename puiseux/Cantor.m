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
forward FunctionsCheck;

forward CantorMorphismFromMatrix;
forward ChangePatchFunctions;
forward AbsoluteToRelative;
forward RelativeToAbsolute;


import "Divisor.m": InitializeCurve, ChangePatchBasisOfDifferentials;
import "LocalInfo.m": PuiseuxRamificationIndex, InitializeImageBranch;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "FractionalCRT.m": FractionalCRTQQ;


function CantorEquations(X);

g := X`g; f := X`DEs[1]; R := X`R;
/* TODO: Generalize this step to arbitrary curves */
if (X`is_hyperelliptic or X`g eq 1) and X`patch_index eq 3 then
    f := Evaluate(f, [R.2, R.1]);
end if;
F := BaseRing(R); S := PolynomialRing(F, 2*g); T<t> := PolynomialRing(S);
//varnames := [ Sprintf("a%o", i) : i in [1..g] ] cat [ Sprintf("b%o", i) : i in [1..g] ];
//AssignNames(~S, varnames);
/* Start with trace and end with norm: */
canpol := t^g + &+[ S.i * t^(g - i) : i in [1..g] ];
substpol := &+[ S.(g + i) * t^(g - i) : i in [1..g] ];
eqpol := Evaluate(f, [t, substpol]) mod canpol;
return Coefficients(eqpol);

end function;


function CandidateFunctions(X, d)

f := X`DEs[1];
R := X`R; x := X`x; y := X`y;
if X`is_hyperelliptic or X`g eq 1 then
    nums := [ x^i : i in [0..(d div 2)] ] cat [ x^i*y : i in [0..((d - X`g - 1) div 2)] ];
    dens := [ x^i : i in [0..(d div 2)] ];
    dens := nums;
elif X`is_planar then
    nums := [ x^i*y^j : i in [0..d], j in [0..(Degree(f, y) - 1)] | i + j le d ];
    dens := [ x^i : i in [0..d] ];
    dens := nums;
end if;
return dens, nums;

end function;


function FunctionValuesFromApproximations(Y, Qs)

if (Y`is_hyperelliptic or Y`g eq 1) and Y`patch_index eq 3 then
    Qs := [ [ Q[2], Q[1] ] : Q in Qs ];
end if;
PR := Parent(Qs[1][1]); R<t> := PolynomialRing(PR);
pol_approx := &*[ t - Q[1] : Q in Qs ];
/* Start with trace and end with norm: */
as_approx := Reverse(Coefficients(pol_approx)[1..Y`g]);
v := Matrix([ [ Q[2] : Q in Qs ] ]);
M := Transpose(Matrix([ [ Q[1]^i : i in [0..(Y`g - 1)] ] : Q in Qs ]));
w := v*M^(-1);
bs_approx := Reverse(Eltseq(w));
return as_approx cat bs_approx;

end function;


function FunctionsFromApproximations(X, Y, P, Qs, d)

/*
// For debugging:
K<y,x> := X`K; R<y,x> := X`R;
fX := R ! X`DEs[1]; fY := R ! Y`DEs[1];
IX := ideal<R | fX>;
fs := [ y / (1 + x + 3*x^2)^2, x^2 / (1 + x + 3*x^2) ];
print "Well-defined?", R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX;
dx := K ! 1;
dy := K ! -Derivative(fX, 2)/Derivative(fX, 1);
ev := ((K ! Derivative(fs[2], 2))*dx + (K ! Derivative(fs[2], 1))*dy) / (K ! -Evaluate(Derivative(fY, 1)/2, fs));
M := Matrix(X`F, [[1,2,0]]);
print "Correct pullback?", R ! Numerator(K ! (ev + (M[1,3] + M[1,2]*x + M[1,1]*x^2)/y)) in IX;
Q := Qs[1];
ev := [ Evaluate(f, P) : f in fs ];
print Evaluate(Y`DEs[1], Q); print Evaluate(Y`DEs[1], [ Evaluate(f, P) : f in fs ]);
print Valuation(Q[1] - ev[1]); print Valuation(Q[2] - ev[2]);
print FractionalCRTQQ([Coefficient(Q[1], 2)], [Characteristic(X`F)]); print FractionalCRTQQ([Coefficient(ev[1], 2)], [Characteristic(X`F)]);
print FractionalCRTQQ([Coefficient(Q[1], 3)], [Characteristic(X`F)]); print FractionalCRTQQ([Coefficient(ev[1], 3)], [Characteristic(X`F)]);
*/

//print Evaluate(X`DEs[1], P); print Evaluate(Y`DEs[1], Qs[1]);
I := ideal<X`R | X`DEs[1]>;
dens, nums := CandidateFunctions(X, d);
fs_approx := FunctionValuesFromApproximations(Y, Qs);
fs := [ ];
test := true;
for f_approx in fs_approx do
    ev_dens := [ -f_approx * Evaluate(den, P) : den in dens ];
    ev_nums := [ Evaluate(num, P) : num in nums ];
    evs := ev_dens cat ev_nums;
    prec := Floor(Minimum([ AbsolutePrecision(ev) : ev in evs ]));
    M := Matrix([ [ X`F ! Coefficient(ev, i) : i in [0..(prec - 1)] ] : ev in evs ]);
    Ker := Kernel(M);
    if Dimension(Ker) eq 0 then
        //print "Dimension 0";
        test := false;
        fs := [];
        break;
        //Append(~fs, X`K ! 0);
    else
        found := false;
        B := Basis(Ker);
        for b in B do
            v := Eltseq(b);
            f := &+[ v[i + #dens]*nums[i] : i in [1..#nums] ] / &+[ v[i]*dens[i] : i in [1..#dens] ];
            //print f;
            //if not (X`R ! Numerator(X`K ! f)) in I then
                Append(~fs, X`K ! f);
                found := true;
                break;
            //end if;
        end for;
        if not found then
            test := false;
            fs := [];
            break;
            //Append(~fs, X`K ! 0);
        end if;
    end if;
end for;
return test, fs;

end function;


function FunctionsCheck(X, Y, fs)

I := ideal<X`R | X`DEs[1]>;
for cantor_eq in Y`cantor_equations do
    if not (X`R ! Numerator(X`K ! Evaluate(cantor_eq, fs))) in I then
        return false;
    end if;
end for;
return true;

end function;


intrinsic CantorMorphismFromMatrix(X::Crv, P0:: Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity()) -> Sch
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent
representation of a projection morphism on the standard basis of differentials,
returns a corresponding Cantor morphism (if it exists). The parameter Margin
specifies how many potentially superfluous terms are used in the development of
the branch, the parameter LowerBound specifies at which degree one starts to
look for a divisor, and the parameter UpperBound specifies where to stop.}

output := InitializeCurve(X, P0); output := InitializeCurve(Y, Q0);
NormM := ChangePatchBasisOfDifferentials(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);
e := PuiseuxRamificationIndex(NormM);

d := LowerBound;
while true do
    vprintf EndoCheck : "Trying degree %o...\n", d;
    dens, nums := CandidateFunctions(X, d);
    n := #dens + #nums + Margin;
    vprintf EndoCheck : "Number of digits in expansion: %o.\n", n;

    /* TODO: This does some work many times over.
     * On the other hand, an iterator also has its disadvantages because of superfluous coefficients. */
    vprintf EndoCheck : "Expanding... ";
    P, Qs := ApproximationsFromTangentAction(X, Y, NormM, n*e);
    vprint EndoCheck, 3 : P, Qs;
    vprintf EndoCheck : "done.\n";

    /* Fit a Cantor morphism to it: */
    vprintf EndoCheck : "Solving linear system... ";
    test, fs := FunctionsFromApproximations(X, Y, P, Qs, d);
    vprintf EndoCheck : "done.\n";

    if test then
        as := fs[1..Y`g]; bs := fs[(Y`g + 1)..(2*Y`g)];
        vprintf EndoCheck : "Checking:\n";
        vprintf EndoCheck : "Step 1... ";
        //test1 := &and[ IsWeaklyZero(Q[1]^g + &+[ Evaluate(as[i], P) * Q[1]^(g - i) : i in [1..g] ]) : Q in Qs ];
        //test2 := &and[ IsWeaklyZero(Q[2]   - &+[ Evaluate(bs[i], P) * Q[1]^(g - i) : i in [1..g] ]) : Q in Qs ];
        //vprintf EndoCheck : "done.\n";
        //if test1 and test2 then
            //vprintf EndoCheck : "Step 2...\n";
            test := FunctionsCheck(X, Y, fs);
            vprintf EndoCheck : "done.\n";
            if test then
                vprintf EndoCheck : "Functions found!\n";
                return true, ChangePatchFunctions(X, Y, fs);
            end if;
        //end if;
    end if;

    /* If that does not work, give up and try one degree higher: */
    d +:= 1;
    if d gt UpperBound then
        return false, "";
    end if;
end while;

end intrinsic;


intrinsic CantorMorphismFromMatrixSplit(X::Crv, P0:: Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), B := 300) -> Sch
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent
representation of a projection morphism on the standard basis of differentials,
returns a corresponding Cantor morphism (if it exists). The parameter Margin
specifies how many potentially superfluous terms are used in the development of
the branch, the parameter LowerBound specifies at which degree one starts to
look for a divisor, and the parameter UpperBound specifies where to stop.}

output := InitializeCurve(X, P0); output := InitializeCurve(Y, Q0);

/*
// For debugging:
K<y,x> := X`K; R<y,x> := X`R;
fX := R ! X`DEs[1]; fY := R ! Y`DEs[1];
print fX;
print M;
IX := ideal<R | fX>;
fs := [ y / (1 + x + 3*x^2)^2, x^2 / (1 + x + 3*x^2) ];
print "Well-defined?", R ! Numerator(K ! Evaluate(R ! fY, fs)) in IX;
dx := K ! 1;
dy := K ! -Derivative(fX, 2)/Derivative(fX, 1);
ev := ((K ! Derivative(fs[2], 2))*dx + (K ! Derivative(fs[2], 1))*dy) / (K ! -Evaluate(Derivative(fY, 1)/2, fs));
print "Correct pullback?", R ! Numerator(K ! (ev + (M[1,3] + M[1,2]*x + M[1,1]*x^2)/y)) in IX;
*/

NormM := ChangePatchBasisOfDifferentials(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);
tjs0, f := InitializeImageBranch(NormM);
e := PuiseuxRamificationIndex(NormM);

/* Some global elements needed below: */
gY := Y`g; F := X`F; rF := X`rF; OF := X`OF; BOF := X`BOF; RX := X`R; KX := X`K;
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, gY + 1);
/* TODO: This is for the test approximation test later on. It makes things a
 * bit less elegant. */
if (Y`is_hyperelliptic or Y`g eq 1) and Y`patch_index eq 3 then
    Qs := [ [ Q[2], Q[1] ] : Q in Qs ];
end if;

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

    /* Add corresponding data: */
    pr := ideal<OF | [ p, rF - rt ]>;
    Append(~prs, pr); I *:= pr;
    X_red := ReduceCurveSplit(X, p, rt); Y_red := ReduceCurveSplit(Y, p, rt);
    NormM_red := ReduceMatrixSplit(NormM, p, rt);
    BI := Basis(I);

    /* Uncomment for check on compatibility with reduction */
    //print CantorMorphismFromMatrix(X_red`U, X_red`P0, (X_red`T)^(-1) * M_red * X_red`T);

    while true do
        vprintf EndoCheck : "Trying degree %o...\n", d;
        dens_red, nums_red := CandidateFunctions(X_red, d);
        n := #dens_red + #nums_red + Margin;
        vprintf EndoCheck : "Number of digits in expansion: %o.\n", n*e;

        /* Take non-zero image branch: */
        /* TODO: This does some work many times over, but only the first time */
        vprintf EndoCheck, 2 : "Expanding... ";
        P_red, Qs_red := ApproximationsFromTangentAction(X_red, Y_red, NormM_red, n*e);
        vprint EndoCheck, 3 : P_red, Qs_red;
        vprintf EndoCheck, 2 : "done.\n";

        /* Fit a Cantor morphism to it: */
        vprintf EndoCheck, 2 : "Solving linear system... ";
        test_red, fs_red := FunctionsFromApproximations(X_red, Y_red, P_red, Qs_red, d);
        vprintf EndoCheck, 2 : "done.\n";

        if test_red then
            vprintf EndoCheck, 2 : "Checking:\n";
            vprintf EndoCheck, 2 : "Step 1: ";
            if not have_to_check then
                vprintf EndoCheck, 2 : "done.\n";
                vprintf EndoCheck, 2 : "Functions found!\n";
                break;
            end if;
            test := FunctionsCheck(X_red, Y_red, fs_red);
            vprintf EndoCheck, 2 : "done.\n";
            if test then
                vprintf EndoCheck, 2 : "Functions found!\n";
                have_to_check := false;
                break;
            end if;
        end if;

        /* If that does not work, give up and try one degree higher. */
        d +:= 1;
        if d gt UpperBound then
            return false, "";
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
    as := fs[1..gY]; bs := fs[(gY + 1)..(2*gY)];
    test1 := true;
    for Q in Qs do
        if not IsWeaklyZero(Q[1]^gY + &+[ Evaluate(as[i], P) * Q[1]^(gY - i) : i in [1..gY] ]) then
            test1 := false;
            break;
        end if;
        if not IsWeaklyZero(Q[2]    - &+[ Evaluate(bs[i], P) * Q[1]^(gY - i) : i in [1..gY] ]) then
            test1 := false;
            break;
        end if;
    end for;
    vprintf EndoCheck : "done.\n";

    if test1 then
        vprintf EndoCheck : "Step 2... ";
        //vprintf EndoCheck : "Step 2...\n";
        //vprintf EndoCheck : "Candidate functions:\n";
        //vprint EndoCheck : fs;
        test2 := FunctionsCheck(X, Y, fs);
        vprintf EndoCheck : "done.\n";
        if test2 then
            vprintf EndoCheck : "Functions found!\n";
            vprintf EndoCheck, 2 : "Degree: %o\n", d;
            vprintf EndoCheck, 2 : "Before changing patch: %o\n", fs;
            return true, ChangePatchFunctions(X, Y, fs);
        end if;
    end if;
end while;

end intrinsic;


function ChangePatchFunctions(X, Y, fs)
/* TODO: Generalize all steps to higher genus when the time comes */

R := X`R; K := X`K;
if X`is_hyperelliptic then
    if X`patch_index eq 3 then
        fs := [ X`K ! Evaluate(f, [ R.2/R.1^(X`g + 1), 1/R.1 ]) : f in fs ];
    end if;
else
    if X`patch_index eq 2 then
        fs := [ X`K ! Evaluate(f, [ R.1/R.2, 1/R.2 ]) : f in fs ];
    elif X`patch_index eq 3 then
        fs := [ X`K ! Evaluate(f, [ R.2/R.1, 1/R.1 ]) : f in fs ];
    end if;
end if;

/* Minus for passage from Cantor to naive */
if Y`g eq 1 then
    if Y`patch_index eq 1 then
        fs := [ -fs[1], fs[2] ];
    elif Y`patch_index eq 3 then
        fs := [ -1/fs[1], fs[2]/fs[1]^2 ];
    end if;
end if;

fs_red := [ ];
if IsAffine(X) then
    A := X;
else
    A := AffinePatch(X, 1);
end if;
S := PolynomialRing(X`F); T := PolynomialRing(S);
DER := R ! DefiningEquations(A)[1]; DET := AbsoluteToRelative(DER, R, T);
for f in fs do
    num := R ! Numerator(f); den := R ! Denominator(f);
    //if Y`patch_index eq 3 then
        den_conj := Evaluate(den, [R.1, -R.2]);
        num *:= den_conj; den *:= den_conj;
    //end if;
    num := AbsoluteToRelative(num, R, T); den := AbsoluteToRelative(den, R, T);
    num := RelativeToAbsolute(num mod DET, R, T); den := RelativeToAbsolute(den mod DET, R, T);
    Append(~fs_red, num / den);
end for;
return fs_red;

end function;


function AbsoluteToRelative(fR, R, T)

S := BaseRing(T);
fT := T ! 0;
for mon in Monomials(fR) do
    c := MonomialCoefficient(fR, mon); exp := Exponents(mon);
    fT +:= c * S.1^exp[1] * T.1^exp[2];
end for;
return fT;

end function;


function RelativeToAbsolute(fT, R, T)

S := BaseRing(T); h := hom<S -> R | [ R.1 ] >;
fR := R ! 0;
for mon in Monomials(fT) do
    c := MonomialCoefficient(fT, mon); exp := Exponents(mon);
    fR +:= h(c) * R.2^exp[1];
end for;
return fR;

end function;
