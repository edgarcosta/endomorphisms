/***
 *  Divisor functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


declare attributes Crv : is_hyperelliptic, is_planar, is_smooth, is_plane_quartic;
declare attributes Crv : unif, unif_index;
declare attributes Crv : g, U, P0, A, DEs;
declare attributes Crv : patch_index, R, x, y, K;
declare attributes Crv : F, rF, OF, BOF;
declare attributes Crv : OurB, NormB, T;
declare attributes Crv : initialized;
declare attributes Crv : cantor_equations;

declare verbose EndoCheck, 3;


forward InitializeCurve;
forward AlgebraicUniformizer;
forward OurBasisOfDifferentials;
forward ChangePatchBasisOfDifferentials;
forward NormalizedBasisOfDifferentials;

forward CandidateDivisors;
forward IrreducibleComponentsFromBranches;
forward IrreducibleComponentCheck;

forward DivisorFromMatrixSplit;


import "LocalInfo.m": DevelopPoint, InitializeImageBranch;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "Cantor.m": CantorEquations;


function InitializeCurve(X, P0)

if not assigned X`initialized then
    X`initialized := false;
end if;
if X`initialized then
    return 0;
end if;
X`is_hyperelliptic := IsHyperelliptic(X); X`is_planar := IsPlaneCurve(X); X`is_smooth := IsNonSingular(X);
X`g := Genus(X); X`is_plane_quartic := (X`is_planar) and (X`is_smooth) and (X`g eq 3);
if IsAffine(X) then
    X`U := X; X`P0 := P0; X`patch_index := 1;
else
    X`U, X`P0, X`patch_index := AffinePatch(X, P0);
end if;
X`A := Ambient(X`U); X`R := CoordinateRing(X`A);
if (X`is_hyperelliptic or X`g eq 1) and (X`patch_index eq 3) then
    X`x := X`R.2; X`y := X`R.1;
else
    X`x := X`R.1; X`y := X`R.2;
end if;
X`F := BaseRing(X`R);
if Type(X`F) eq FldRat then
    X`rF := 1;
    X`OF := Integers();
else
    X`rF := Denominator(X`F.1) * X`F.1;
    X`OF := Order([ X`rF^i : i in [0..Degree(X`F) - 1] ]);
end if;
X`BOF := Basis(X`OF); X`K := FieldOfFractions(X`R);
X`DEs := DefiningEquations(X`U);
X`unif, X`unif_index := AlgebraicUniformizer(X);
X`OurB := OurBasisOfDifferentials(X);
X`NormB, X`T := NormalizedBasisOfDifferentials(X);
if X`is_planar then
    X`cantor_equations := CantorEquations(X);
end if;
X`initialized := true;
return 0;

end function;


function AlgebraicUniformizer(X)
/*
 * Input:   A curve X.
 * Output:  A uniformizing element at P0
 *          and the corresponding index.
 */

fX := X`DEs[1]; x := X`x; y := X`y;
if Evaluate(Derivative(fX, y), X`P0) ne 0 then
    if X`R.1 eq x then
        return x, 1;
    else
        return x, 2;
    end if;
end if;
if X`R.1 eq y then
    return y, 1;
else
    return y, 2;
end if;

/* TODO: Previous functionality */
Gens := GeneratorsSequence(X`R);
M := Matrix([ [ Evaluate(Derivative(DE, gen), X`P0) : gen in Gens ] : DE in X`DEs ]);
/* Default is the first coordinate: */
i0 := 1;
for i in [1..#Gens] do
    if &and[ M[j, i] eq 0 : j in [1..#Rows(M)] ] then
        i0 := i;
    end if;
end for;
return Gens[i0], i0;

end function;


function OurBasisOfDifferentials(X);
/*
 * Input:   A curve X.
 * Output:  A basis of global differentials on X, represented by elements of
 *          the ambient.
 */

g := X`g;
R := X`R; x := X`x; y := X`y;
f := X`DEs[1];
s := MonomialCoefficient(f, y^2);
if g eq 0 then
    return [ ];
elif X`is_hyperelliptic or (g eq 1) then
    if g eq 1 then
        if X`unif eq x then
            return [ 2*s/Derivative(f, y) ];
        else
            return [ -2*s/Derivative(f, x) ];
        end if;
    end if;
    /* Otherwise the uniformizer will always be x */
    return [ 2*s*x^(i-1) / Derivative(f, y) : i in [1..g] ];
elif X`is_plane_quartic then
    /* TODO: Analyze sign */
    if X`unif_index eq 1 then
        return [ X`K ! (n / Derivative(f, 2)) : n in [R.1, R.2, 1] ];
    else
        return [ X`K ! (n / Derivative(f, 1)) : n in [R.1, R.2, 1] ];
    end if;
else
    /* TODO: Compatibility problems arise here */
    B := BasisOfDifferentialsFirstKind(X`U);
    du := Differential(AlgebraicUniformizer(X));
    return [ X`K ! (b / du) : b in B ];
end if;

end function;


function ChangePatchBasisOfDifferentials(X, Y, M)

/* TODO: Assumes that elliptic curve given as hyperelliptic curve */
if X`g eq 1 or X`is_hyperelliptic then
    if X`patch_index eq 3 then
        M := Matrix(X`F, [ Reverse([ -c : c in Eltseq(row)]) : row in Rows(M) ]);
    end if;
elif X`is_plane_quartic then
    if X`patch_index eq 2 then
        M := Matrix(X`F, [ [ row[1], row[3], row[2] ] : row in Rows(M) ]);
    elif X`patch_index eq 3 then
        M := Matrix(X`F, [ [ row[3], row[2], row[1] ] : row in Rows(M) ]);
    end if;
end if;

M := Transpose(M);
if Y`g eq 1 or Y`is_hyperelliptic then
    if Y`patch_index eq 3 then
        M := Matrix(Y`F, [ Reverse([ -c : c in Eltseq(row)]) : row in Rows(M) ]);
    end if;
elif Y`is_plane_quartic then
    if Y`patch_index eq 2 then
        M := Matrix(Y`F, [ [ row[1], row[3], row[2] ] : row in Rows(M) ]);
    elif Y`patch_index eq 3 then
        M := Matrix(Y`F, [ [ row[3], row[2], row[1] ] : row in Rows(M) ]);
    end if;
end if;

M := Transpose(M);
return M;

end function;


function NormalizedBasisOfDifferentials(X)
/*
 * Input:   A curve X.
 * Output:  A differential basis Bnorm that is normalized with respect to the uniformizing parameter,
 *          and a matrix T such that multiplication by T sends B to Bnorm.
 */

P := DevelopPoint(X, X`P0, X`g);
//print X`U;
//print X`unif;
BP := [ Evaluate(b, P) : b in X`OurB ];
T := Matrix([ [ Coefficient(BP[i], j - 1) : j in [1..X`g] ] : i in [1..X`g] ])^(-1);
NormB := [ &+[ T[i,j] * X`OurB[j] : j in [1..X`g] ] : i in [1..X`g] ];
return NormB, T;

end function;


function CandidateDivisors(X, Y, d)
/*
 * Input:   Two curves X and Y
 *          and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

gX := X`g; fX := X`DEs[1]; RX := X`R; xX := X`x; yX := X`y;
gY := Y`g; fY := Y`DEs[1]; RY := Y`R; xY := Y`x; yY := Y`y;
R := PolynomialRing(X`F, 2, "lex"); Rprod := PolynomialRing(X`F, 4, "lex");

if X`is_hyperelliptic then
    divsX := [ xX^i : i in [0..(d div 2)] ] cat [ xX^i*yX : i in [0..((d - gX - 1) div 2)] ];
    divsY := [ xY^i : i in [0..(d div 2)] ] cat [ xY^i*yY : i in [0..((d - gX - 1) div 2)] ];
    //divsY := [ xY^i : i in [0..gY] ] cat [ yY ];
    Reverse(~divsX); Reverse(~divsY);
elif X`is_planar then
    divsX := [ xX^i*yX^j : i in [0..d], j in [0..(Degree(fX, yX) - 1)] | i + j le d ];
    divsY := [ xY^i*yY^j : i in [0..gY], j in [0..(Degree(fY, yY) - 1)] ];
    divsY := [ xY^i : i in [0..gY] ] cat [ yY ];
    Reverse(~divsX); Reverse(~divsY);
end if;
//print divsX; print divsY;

divsX := [ R ! d : d in divsX ]; divsY := [ R ! d : d in divsY ];
//hX := hom<RX -> Rprod | [ Rprod.1, Rprod.2 ]>; hY := hom<RY -> Rprod | [ Rprod.3, Rprod.4 ]>;
hX := hom<RX -> Rprod | [ Rprod.3, Rprod.1 ]>; hY := hom<RY -> Rprod | [ Rprod.4, Rprod.2 ]>;
hs := [ hX, hY ];
CP := CartesianProduct([ divsX, divsY ]);
divs := [ &*[ hs[i](tup[i]) : i in [1..2] ] : tup in CP ];
divs := Reverse(Sort(divs));
//print divs;
return divs;

end function;


function IrreducibleComponentsFromBranches(X, Y, fs, P, Qs)
/*
 * Input:   Two curves X and Y,
 *          a basis of divisor equations fs,
 *          the precision n used when determining these,
 *          and branch expansions P and Qs.
 * Output:  The irreducible components corresponding that fit the given data.
 */

/* Recovering a linear system: */
e := Maximum([ Maximum([ Denominator(Valuation(c - Coefficient(c, 0))) : c in Q ]) : Q in Qs ]);
prec := Precision(Parent(P[1]));
M := [ ];
for f in fs do
    r := [ ];
    for Q in Qs do
        //ev := Evaluate(f, P cat Q);
        ev := Evaluate(f, [ Q[2], P[2], Q[1], P[1] ]);
        r cat:= [ Coefficient(ev, i/e) : i in [0..prec - X`g] ];
    end for;
    Append(~M, r);
end for;
M := Matrix(M);
B := Basis(Kernel(M));

/* Coerce back to ground field (possible because of echelon form): */
B := [ [ X`F ! c : c in Eltseq(b) ] : b in B ];

/* Corresponding equations: */
Rprod := Parent(fs[1]);
//hX := hom<X`R -> Rprod | [ Rprod.1, Rprod.2 ]>; hY := hom<Y`R -> Rprod | [ Rprod.3, Rprod.4 ]>;
hX := hom<X`R -> Rprod | [ Rprod.4, Rprod.2 ]>; hY := hom<Y`R -> Rprod | [ Rprod.3, Rprod.1 ]>;
eqs := [ &+[ b[i] * fs[i] : i in [1..#fs] ] : b in B ];
eqs := eqs cat [ hX(DE) : DE in X`DEs ] cat [ hY(DE) : DE in Y`DEs ];

vprintf EndoCheck, 2 : "Calculating Groebner basis... ";
GB := GroebnerBasis(ideal< Rprod | eqs >);
vprint EndoCheck, 3 : GB;
vprintf EndoCheck, 2 : "done.\n";
Append(~eqs, GB[#GB]);

/* Corresponding scheme: */
A := AffineSpace(Rprod);
S := Scheme(A, eqs);

return [ S ];

/* TODO: These steps may be a time sink and should be redundant, so we avoid
 *       them. They get eliminated as the degree increases anyway. */
Is := IrreducibleComponents(S);
return [ ReducedSubscheme(I) : I in Is ];

end function;


function IrreducibleComponentCheck(X, Y, I)
/*
 * Input:   An irreducible scheme I in X x Y.
 * Output:  Whether or not I intersects P0 x X with the correct multiplicity at
 *          P0 and nowhere else.
 */

A4 := Ambient(I); R4 := CoordinateRing(A4);
R2 := PolynomialRing(X`F, 2); A2 := AffineSpace(R2);
//h := hom< R4 -> R2 | [ X`P0[i] : i in [1..2] ] cat [ R2.i : i in [1..2] ] >;
h := hom< R4 -> R2 | [ R2.2, X`P0[2], R2.1, X`P0[1] ] >;
eqs2 := [ h(eq4) : eq4 in DefiningEquations(I) ];
S := Scheme(A2, eqs2);
if Dimension(S) eq 0 then
    if Degree(ReducedSubscheme(S)) eq 1 then
        if Degree(S) eq Y`g then
            /* TODO: This is potentially slightly unsafe but delivers a big speedup */
            //if Dimension(I) eq 1 then
                return true;
            //end if;
        end if;
    end if;
end if;
return false;

end function;


intrinsic DivisorFromMatrix(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity()) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent
representation of a projection morphism on the standard basis of differentials,
returns a corresponding divisor (if it exists). The parameter Margin specifies
how many potentially superfluous terms are used in the development of the
branch, the parameter LowerBound specifies at which degree one starts to look
for a divisor, and the parameter UpperBound specifies where to stop.}

output := InitializeCurve(X, P0); output := InitializeCurve(Y, Q0);
NormM := ChangePatchBasisOfDifferentials(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);

d := LowerBound;
while true do
    vprintf EndoCheck : "Trying degree %o...\n", d;
    fs := CandidateDivisors(X, Y, d);
    n := #fs + Margin;
    vprintf EndoCheck : "Number of digits in expansion: %o.\n", n;

    /* Take non-zero image branch: */
    vprintf EndoCheck : "Expanding... ";
    P, Qs := ApproximationsFromTangentAction(X, Y, NormM, n);
    vprint EndoCheck, 3 : P, Qs;
    vprintf EndoCheck : "done.\n";

    /* Fit a divisor to it: */
    vprintf EndoCheck : "Solving linear system... ";
    ICs := IrreducibleComponentsFromBranches(X, Y, fs, P, Qs);
    vprintf EndoCheck : "done.\n";

    for S in ICs do
        DEs := DefiningEquations(S);
        vprintf EndoCheck : "Checking:\n";
        vprintf EndoCheck : "Step 1... ";
        //test1 := &and[ &and[ IsWeaklyZero(Evaluate(DE, P cat Q)) : Q in Qs ] : DE in DEs ];
        //vprintf EndoCheck : "done.\n";
        //if test1 then
            //vprintf EndoCheck : "Step 2... ";
            test2 := IrreducibleComponentCheck(X, Y, S);
            vprintf EndoCheck : "done.\n";
            if test2 then
                vprintf EndoCheck : "Divisor found!\n";
                return true, S;
            end if;
        //end if;
    end for;

    /* If that does not work, give up and try one degree higher: */
    d +:= 1;
    if d gt UpperBound then
        return false, "";
    end if;
end while;

end intrinsic;


intrinsic DivisorFromMatrixSplit(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), B := 300) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent
representation of a projection morphism on the standard basis of differentials,
returns a corresponding divisor (if it exists). The parameter Margin specifies
how many potentially superfluous terms are used in the development of the
branch, the parameter LowerBound specifies at which degree one starts to look
for a divisor, and the parameter UpperBound specifies where to stop.}

/* We start at a suspected estimate and then increase degree until we find an appropriate divisor: */
output := InitializeCurve(X, P0); output := InitializeCurve(Y, Q0);
NormM := ChangePatchBasisOfDifferentials(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);
tjs0, f := InitializeImageBranch(NormM);
F := X`F; rF := X`rF; OF := X`OF; BOF := X`BOF;
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, X`g);
Rprod := PolynomialRing(X`F, 4, "lex");

ps_rts := [ ]; prs := [ ]; DEss_red := [* *];
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
    pr := ideal<X`OF | [ p, rF - rt ]>;
    Append(~prs, pr); I *:= pr;
    X_red := ReduceCurveSplit(X, p, rt); Y_red := ReduceCurveSplit(Y, p, rt);
    NormM_red := ReduceMatrixSplit(NormM, p, rt);
    BI := Basis(I);

    /* Uncomment for check on compatibility with reduction */
    //print DivisorFromMatrix(X_red`U, X_red`P0, (X_red`T)^(-1) * NormM_red * X_red`T);

    done := false;
    while true do
        vprintf EndoCheck : "Trying degree %o...\n", d;
        fs_red := CandidateDivisors(X_red, Y_red, d);
        n := #fs_red + Margin;
        vprintf EndoCheck : "Number of digits in expansion: %o.\n", n;

        /* Take non-zero image branch: */
        vprintf EndoCheck, 2 : "Expanding... ";
        P_red, Qs_red := ApproximationsFromTangentAction(X_red, Y_red, NormM_red, n);
        vprint EndoCheck, 3 : P_red, Qs_red;
        vprintf EndoCheck, 2 : "done.\n";

        /* Fit a divisor to it: */
        vprintf EndoCheck, 2 : "Solving linear system... ";
        ICs_red := IrreducibleComponentsFromBranches(X_red, Y_red, fs_red, P_red, Qs_red);
        vprintf EndoCheck, 2 : "done.\n";

        for S_red_it in ICs_red do
            vprintf EndoCheck, 2 : "Checking:\n";
            vprintf EndoCheck, 2 : "Step 1... ";
            if not have_to_check then
                vprintf EndoCheck, 2 : "done.\n";
                vprintf EndoCheck, 2 : "Divisor found!\n";
                done := true;
                S_red := S_red_it;
                break;
            end if;
            test := IrreducibleComponentCheck(X_red, Y_red, S_red_it);
            vprintf EndoCheck, 2 : "done.\n";
            if test then
                vprintf EndoCheck, 2 : "Divisor found!\n";
                done := true;
                S_red := S_red_it;
                have_to_check := false;
                break;
            end if;
        end for;

        if done then
            break;
        end if;

        /* If that does not work, give up and try one degree higher.
         * Note that d is initialized in the outer loop,
         * so that we keep the degree that works. */
        d +:= 1;
        if d gt UpperBound then
            return false, "";
        end if;
    end while;
    Append(~DEss_red, DefiningEquations(S_red));

    vprintf EndoCheck : "Fractional CRT... ";
    DEs := [ ];
    for i:=1 to #DEss_red[1] do
        DE := Rprod ! 0;
        for mon in Monomials(DEss_red[1][i]) do
            exp := Exponents(mon);
            rs := [* *];
            for j:=1 to #DEss_red do
                Rprod_red := Parent(DEss_red[j][1]);
                Append(~rs, MonomialCoefficient(DEss_red[j][i], Monomial(Rprod_red, exp)));
            end for;
            DE +:= FractionalCRTSplit(rs, prs, OF, I, BOF, BI, F) * Monomial(Rprod, exp);
        end for;
        Append(~DEs, DE);
    end for;
    vprintf EndoCheck : "done.\n";

    vprintf EndoCheck : "Checking:\n";
    vprintf EndoCheck : "Step 1... ";
    test1 := true;
    for DE in DEs do
        test1_int := true;
        for Q in Qs do
            //if not IsWeaklyZero(Evaluate(DE, P cat Q)) then
            if not IsWeaklyZero(Evaluate(DE, [ Q[2], P[2], Q[1], P[1] ])) then
                test1_int := false;
                break;
            end if;
        end for;
        if not test1_int then
            test1 := false;
            break;
        end if;
    end for;
    vprintf EndoCheck : "done.\n";

    if test1 then
        S := Scheme(AffineSpace(Rprod), DEs);
        vprintf EndoCheck : "Step 2... ";
        //vprintf EndoCheck : "Step 2...\n";
        //vprintf EndoCheck : "Candidate divisor:\n";
        //vprint EndoCheck : S;
        test2 := IrreducibleComponentCheck(X, Y, S);
        //test2 := true;
        vprintf EndoCheck : "done.\n";
        if test2 then
            vprintf EndoCheck : "Divisor found!\n";
            return true, S;
        end if;
    end if;
end while;

end intrinsic;
