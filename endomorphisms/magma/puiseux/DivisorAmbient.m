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


import "Branches.m": InitializeImageBranch, DevelopPoint;
import "Conventions.m": ExtractHomomorphismsRing, VariableOrder, ExtractPoints;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "Initialize.m": InitializeCurve, ChangeTangentAction;
import "RiemannRoch.m": RRBasis, RREvaluations, ProductEvaluations, GlobalProductBasis, GlobalProductBasisAlt;


forward CandidateDivisors;
forward IrreducibleComponentsFromBranches;
forward CheckEquations;
forward CheckIrreducibleComponent;

forward DivisorFromMatrix;
forward DivisorFromMatrixSplit;
forward DivisorFromMatrixByDegree;


function CandidateDivisors(X, Y, d)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

gX := X`g; fX := X`DEs[1]; RX := X`RA;
xX := RX.1; yX := RX.2;
/* Change in hyperelliptic case for greater effectiveness: */
if Degree(fX, RX.1) eq 2 then
    xX := RX.2; yX := RX.1;
end if;
gY := Y`g; fY := Y`DEs[1]; RY := Y`RA;
xY := RY.1; yY := RY.2;
/* Change in hyperelliptic case for greater effectiveness: */
if Degree(fY, RY.1) eq 2 then
    xY := RY.2; yY := RY.1;
end if;

if X`is_hyperelliptic then
    divsX := [ xX^i : i in [0..(d div 2)] ] cat [ xX^i*yX : i in [0..((d - gX - 1) div 2)] ];
    divsY := [ xY^i : i in [0..(d div 2)] ] cat [ xY^i*yY : i in [0..((d - gX - 1) div 2)] ];
    divsY := [ xY^i : i in [0..gY] ] cat [ yY ];
    //Reverse(~divsX); Reverse(~divsY);
elif X`is_planar then
    divsX := [ xX^i*yX^j : i in [0..d] , j in [0..(Degree(fX, yX) - 1)] | i + j le d  ];
    divsY := [ xY^i*yY^j : i in [0..gY], j in [0..(Degree(fY, yY) - 1)] | i + j le gY ];
    divsY := [ xY^i : i in [0..gY] ] cat [ yY ];
    //Reverse(~divsX); Reverse(~divsY);
end if;

hX, hY := ExtractHomomorphismsRing(X, Y); hs := [ hX, hY ];
CP := [ [* divX, divY *] : divX in divsX, divY in divsY ];
divs := [ &*[ hs[i](tup[i]) : i in [1..2] ] : tup in CP ];
//divs := Reverse(Sort(divs));
return divs;

end function;


function IrreducibleComponentsFromBranches(X, Y, fs, P, Qs : DivPP1 := false)
/*
 * Input:   Two curves X and Y,
 *          a basis of divisor equations fs,
 *          the precision n used when determining these,
 *          and branch expansions P and Qs.
 * Output:  The irreducible components that fit the given data.
 */

/* Recovering a linear system */
e := Maximum([ Maximum([ Denominator(Valuation(c - Coefficient(c, 0))) : c in Q ]) : Q in Qs ]);
prec := Precision(Parent(P[1]));
M := [ ];
for f in fs do
    r := [ ];
    for Q in Qs do
        seq := ExtractPoints(X, Y, P, Q);
        ev := Evaluate(f, seq);
        r cat:= [ Coefficient(ev, i/e) : i in [0..prec - X`g] ];
    end for;
    Append(~M, r);
end for;
B := Basis(Kernel(Matrix(M)));
/* Coerce back to ground field (possible because of echelon form) */
B := [ [ X`F ! c : c in Eltseq(b) ] : b in B ];

/* Corresponding equations */
hX, hY := ExtractHomomorphismsRing(X, Y);
Rprod := Codomain(hX);
eqs := [ Rprod ! (&+[ b[i] * fs[i] : i in [1..#fs] ]) : b in B ];
eqs := eqs cat [ hX(DE) : DE in X`DEs ] cat [ hY(DE) : DE in Y`DEs ];

if DivPP1 then
    Rprod := Codomain(hX);
    GB := GroebnerBasis(ideal< Rprod | eqs >);
    Append(~eqs, GB[#GB]);
end if;

/* Corresponding scheme */
A := AffineSpace(Rprod);
S := Scheme(A, eqs);
return [ S ];

/* TODO: These steps may be a time sink and should be redundant, so we avoid
 *       them. They get eliminated as the degree increases anyway. */
return [ ReducedSubscheme(I) : I in IrreducibleComponents(S) ];

end function;


function CheckEquations(X, Y, P, Qs, DEs)

for DE in DEs do
    for Q in Qs do
        seq := ExtractPoints(X, Y, P, Q);
        if not IsWeaklyZero(Evaluate(DE, seq)) then
            return false;
        end if;
    end for;
end for;
return true;

end function;


function CheckIrreducibleComponent(X, Y, I)
/*
 * Input:   An irreducible scheme I in X x Y.
 * Output:  Whether or not I intersects P0 x X with the correct multiplicity at
 *          P0 and nowhere else.
 */

A4 := Ambient(I); R4 := CoordinateRing(A4);
R2 := PolynomialRing(X`F, 2); A2 := AffineSpace(R2);
seq := [ X`P0[1], X`P0[2], R2.1, R2.2 ];
varord := VariableOrder();
seq := [ seq[varord[i]] : i in [1..#varord] ];
h := hom< R4 -> R2 | seq >;
eqs2 := [ h(eq4) : eq4 in DefiningEquations(I) ];
S := Scheme(A2, eqs2);

if Dimension(S) eq 0 then
    if Degree(ReducedSubscheme(S)) eq 1 then
        if Degree(S) eq Y`g then
            //if Dimension(I) eq 1 then
                return true;
            //end if;
        end if;
    end if;
end if;
return false;

end function;


intrinsic DivisorFromMatrix(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), DivPP1 := false) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

InitializeCurve(X, P0); InitializeCurve(Y, Q0);
NormM := ChangeTangentAction(X, Y, M);
vprintf EndoCheck, 3 : "Tangent representation:\n";
vprint EndoCheck, 3 : NormM;
NormM := Y`T * NormM * (X`T)^(-1);
vprintf EndoCheck, 3 : "Normalized tangent representation:\n";
vprint EndoCheck, 3 : NormM;

d := LowerBound;
while true do
    found, S := DivisorFromMatrixByDegree(X, Y, NormM, d : Margin := 2^4, DivPP1 := DivPP1, have_to_check := true);
    if found then
        return true, S;
    end if;
    /* If that does not work, give up and try one degree higher */
    d +:= 1;
    if d gt UpperBound then
        return false, [];
    end if;
end while;

end intrinsic;


intrinsic DivisorFromMatrixSplit(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), DivPP1 := false, B := 300) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

/* We start at a suspected estimate and then increase degree until we find an appropriate divisor */
InitializeCurve(X, P0); InitializeCurve(Y, Q0);
NormM := ChangeTangentAction(X, Y, M);
vprintf EndoCheck, 3 : "Tangent representation:\n";
vprint EndoCheck, 3 : NormM;
NormM := Y`T * NormM * (X`T)^(-1);
vprintf EndoCheck, 3 : "Normalized tangent representation:\n";
vprint EndoCheck, 3 : NormM;
tjs0, f := InitializeImageBranch(NormM);

/* Some global elements needed below */
F := X`F; rF := X`rF; OF := X`OF; BOF := X`BOF;
Rprod := PolynomialRing(X`F, 4, "lex");
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, X`g);

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

    /* Add corresponding data */
    pr := ideal<X`OF | [ p, rF - rt ]>;
    Append(~prs, pr); I *:= pr;
    X_red := ReduceCurveSplit(X, p, rt); Y_red := ReduceCurveSplit(Y, p, rt);
    NormM_red := ReduceMatrixSplit(NormM, p, rt);
    BI := Basis(I);

    while true do
        found, S_red := DivisorFromMatrixByDegree(X_red, Y_red, NormM_red, d : Margin := Margin, DivPP1 := DivPP1, have_to_check := have_to_check);
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
    have_to_check := false;
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
    /* Note that P and Qs are calculated at the beginning of this function */
    test1 := CheckEquations(X, Y, P, Qs, DEs);
    vprintf EndoCheck : "done.\n";

    if test1 then
        vprintf EndoCheck : "Step 2... ";
        S := Scheme(AffineSpace(Rprod), DEs);
        test2 := CheckIrreducibleComponent(X, Y, S);
        vprintf EndoCheck : "done.\n";
        if test2 then
            vprintf EndoCheck : "Divisor found!\n";
            return true, S;
        end if;
    end if;
end while;

end intrinsic;


function DivisorFromMatrixByDegree(X, Y, NormM, d : Margin := 2^4, DivPP1 := false, have_to_check := true)

vprintf EndoCheck, 2 : "Trying degree %o...\n", d;
fs := CandidateDivisors(X, Y, d);
//fsNew := CandidateDivisorsNew(X, Y, d);
n := #fs + Margin;
//n := 450;
vprintf EndoCheck, 2 : "Number of terms in expansion: %o.\n", n;

/* Take non-zero image branch */
vprintf EndoCheck, 2 : "Expanding... ";
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, n);
vprintf EndoCheck, 4 : "Base point:\n";
_<t> := Parent(P[1]);
_<r> := BaseRing(Parent(P[1]));
vprint EndoCheck, 4 : P;
vprintf EndoCheck, 4 : "Resulting branches:\n";
vprint EndoCheck, 4 : Qs;
vprint EndoCheck, 4 : BaseRing(Parent(P[1]));
vprintf EndoCheck, 2 : "done.\n";

/* Fit a divisor to it */
vprintf EndoCheck, 2 : "Solving linear system... ";
ICs := IrreducibleComponentsFromBranches(X, Y, fs, P, Qs : DivPP1 := DivPP1);
//ICsNew := IrreducibleComponentsFromBranchesNew(X, Y, fsNew, P, Qs : DivPP1 := DivPP1);
vprintf EndoCheck, 2 : "done.\n";

for S in ICs do
    DEs := DefiningEquations(S);
    vprintf EndoCheck, 2 : "Checking:\n";
    vprintf EndoCheck, 2 : "Step 1... ";
    test1 := CheckEquations(X, Y, P, Qs, DEs);
    vprintf EndoCheck, 2 : "done.\n";
    if test1 then
        vprintf EndoCheck, 2 : "Step 2... ";
        test2 := CheckIrreducibleComponent(X, Y, S);
        vprintf EndoCheck, 2 : "done.\n";
        if test2 then
            vprintf EndoCheck, 2 : "Divisor found!\n";
            return true, S;
        end if;
    end if;
end for;
return false, [ ];

end function;
