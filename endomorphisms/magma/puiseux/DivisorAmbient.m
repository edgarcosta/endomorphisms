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


/* TODO: These algorithms will be redundant soon */
import "Branches.m": InitializeImageBranch, DevelopPoint;
import "Conventions.m": ExtractHomomorphismsRing, VariableOrder, ExtractPoints;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "Initialize.m": InitializeCurve, ChangeTangentAction;
import "RiemannRoch.m": RRBasis, RREvaluations, ProductEvaluations, GlobalProductBasis, GlobalProductBasisAlt;


forward CandidateDivisors;
forward IrreducibleComponentsFromBranches;
forward CheckEquations;
forward CheckIrreducibleComponent;

forward DivisorFromMatrixAmbientGlobal;
forward DivisorFromMatrixAmbientSplit;
forward DivisorFromMatrixByDegree;


function CandidateDivisors(X, Y, d)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

gX := X`g; fX := X`DEs[1]; RX := X`RA;
xX := RX.1; yX := RX.2;
/* Change in elliptic case for greater effectiveness: */
if Degree(fX, RX.1) eq 2 then
    xX := RX.2; yX := RX.1;
end if;

gY := Y`g; fY := Y`DEs[1]; RY := Y`RA;
xY := RY.1; yY := RY.2;
/* Change in elliptic case for greater effectiveness: */
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

function _PuiseuxCoefficients(f, e, range)
    return [ Coefficient(f, i/e) : i in range ];
    // Q[i] == Coefficient(f, (v + i - 1)/d);
    Q, v, d := Coefficients(f);
    coeffs := [
        jd mod e ne 0 select 0 else (
            (idx ge 1 and idx le L) select Q[idx] else 0 where idx:=jd div e + 1 - v)
        where jd := j*d
                                                              :  j in range]
    where L := #Q;
    return coeffs;
end function;
    /*
    if slow ne coeffs then
        slow;
        coeffs;
        assert false;
    end if;
    assert slow eq coeffs;
    */

function IrreducibleComponentsFromBranches(X, Y, fs, P, Qs : Margin := 2^4)
/*
 * Input:   Two curves X and Y,
 *          a basis of divisor equations fs,
 *          the precision n used when determining these,
 *          and branch expansions P and Qs.
 * Output:  The irreducible components that fit the given data.
 */

/* Recovering a linear system */
e := Maximum(&cat[ [ ExponentDenominator(c) : c in Q ] : Q in Qs ]);
prec := Precision(Parent(P[1]));
vprintf EndoCheck, 3 : "Build M...";
//SetProfile(true);
vtime EndoCheck, 3:
M := [ &cat[ _PuiseuxCoefficients(ev, e, range)
    where ev := Evaluate(f, PQ) : PQ  in PQs] : f in fs]
    where PQs := [ExtractPoints(P, Q) : Q in Qs]
    where range:=[0 .. prec - 1 - Margin];
//SetProfile(false);
//ProfilePrintByTotalTime(ProfileGraph());
/* 
M := [ ];
for f in fs do
    r := [ ];
    for Q in Qs do
        seq := ExtractPoints(X, Y, P, Q);
        ev := Evaluate(f, ExtractPoints(X, Y, P, Q));
        r cat:= [ Coefficient(ev, i/e) : i in [0..prec - 1 - Margin] ];
    end for;
    Append(~M, r);
end for;
*/
vprintf EndoCheck, 2 : "Really Solving linear system %o x %o...", #M, #M[1];
vtime EndoCheck, 2:
B := Basis(Kernel(Matrix(M)));
/* Coerce back to ground field (possible because of echelon form) */
vprintf EndoCheck, 3 : "Coerce back to ground field...";
vtime EndoCheck, 3:
B := [ [ X`F ! c : c in Eltseq(b) ] : b in B ];

/* Corresponding equations */
vprintf EndoCheck, 3 : "ExtractHomomorphismsRing...";
vtime EndoCheck, 3:
hX, hY, hxs, hxsinv := ExtractHomomorphismsRing(X, Y);
Rprod := Codomain(hX);
vprintf EndoCheck, 3 : "eqs1...";
vtime EndoCheck, 3:
eqs := [ Rprod ! (&+[ b[i] * fs[i] : i in [1..#fs] ]) : b in B ];
vprintf EndoCheck, 3 : "eqs2...";
vtime EndoCheck, 3:
eqs := eqs cat [ hX(DE) : DE in X`DEs ] cat [ hY(DE) : DE in Y`DEs ];

if X`is_hyperelliptic and Y`is_hyperelliptic then
    Iprod := ideal< Rprod | eqs >;
    varord := VariableOrder;
    Bxs := Basis(EliminationIdeal(Iprod, { varord[1], varord[3] }));
    if (#Bxs ne 0) then
        gcdx := hxsinv(GCD([ hxs(b) : b in Bxs ]));
        eqs cat:= [ gcdx ];
    else
        eqs := [ Rprod ! 0 ];
    end if;
end if;

/* Corresponding scheme */
A := AffineSpace(Rprod); S := Scheme(A, eqs);
return [ S ];
return [ ReducedSubscheme(I) : I in IrreducibleComponents(S) ];

end function;


function CheckEquations(X, Y, P, Qs, DEs)

vprint EndoCheck, 3 : "";
vprint EndoCheck, 3 : "Check approximations zero:";
for DE in DEs do
    for Q in Qs do
        ev := Evaluate(DE, ExtractPoints(P, Q));
        vprint EndoCheck, 3 : ev;
        if not IsWeaklyZero(ev) then
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
varord := VariableOrder;
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


intrinsic DivisorFromMatrixAmbientGlobal(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^5, LowerBound := 1, UpperBound := Infinity()) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials (all over the same base ring), returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

require BaseRing(X) eq BaseRing(Y) : "Expected both curves to be defined over the same ring";
require BaseRing(X) eq BaseRing(M) : "Expected the tangent represeantion to be defined over the same ring as the curves";
require X eq Curve(P0) : "The point P0 is supposed to be a point in X";
require Y eq Curve(Q0) : "The point Q0 is supposed to be a point in Y";

/* The order below matters because of our initialization conventions */
vprintf EndoCheck, 3: "InitializeCurve(Y, Q0 : AssertNonWP := true)...";
vtime EndoCheck, 3:
InitializeCurve(Y, Q0 : AssertNonWP := true);
vprintf EndoCheck, 3: "InitializeCurve(X, P0)...";
vtime EndoCheck, 3:
InitializeCurve(X, P0);
/* Correct for patches and uniformization indices */
NormM := ChangeTangentAction(X, Y, M);
/* Correct for echelonization */
NormM := Y`T * NormM * (X`T)^(-1);

d := LowerBound;
Iterator := InitializedIterator(X, Y, NormM, X`g + 3);
while true do
    found, S, Iterator := DivisorFromMatrixByDegree(X, Y, Iterator, d : Margin := Margin);
    if found then
        return true, S;
    end if;
    if d ge UpperBound then
        return false, [];
    end if;
    d := Max(2*d, UpperBound);
end while;

end intrinsic;


intrinsic DivisorFromMatrixAmbientSplit(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^5, LowerBound := 1, UpperBound := Infinity(), B := 300) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials (all over the same base ring), returns a corresponding divisor (if it exists), constructed via a CRT approach.
    The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

require BaseRing(X) eq BaseRing(Y) : "Expected both curves to be defined over the same ring";
require BaseRing(X) eq BaseRing(M) : "Expected the tangent represeantion to be defined over the same ring as the curves, consider using Correspondence with Al:=\"Divisor\"";
require X eq Curve(P0) : "The point P0 is supposed to be a point in X";
require Y eq Curve(Q0) : "The point Q0 is supposed to be a point in Y";
/* The order below matters because of our initialization conventions */
vprintf EndoCheck, 3: "InitializeCurve(Y, Q0 : AssertNonWP := true)...";
vtime EndoCheck, 3:
InitializeCurve(Y, Q0 : AssertNonWP := true);
vprintf EndoCheck, 3: "InitializeCurve(X, P0)...";
vtime EndoCheck, 3:
InitializeCurve(X, P0);
/* Correct for patches and uniformization indices */
NormM := ChangeTangentAction(X, Y, M);
/* Correct for echelonization */
NormM := Y`T * NormM * (X`T)^(-1);

/* Some global elements needed below */
F := X`F;
OF := X`OF;
Rprod := PolynomialRing(X`F, 4, "lex");
Iterator, f := InitializedIterator(X, Y, NormM, X`g + 3);
P, Qs := Explode(Iterator);

prs := [ ]; DEss_red := [* *];
I := ideal<X`OF | 1>;

d := LowerBound;
while true do
    /* Find new prime */
    repeat
        pr, h := RandomSplitPrime(f, B);
    until not pr in prs;
    Append(~prs, pr); I *:= pr;
    vprint EndoCheck : "";
    vprint EndoCheck : "Split prime over", #Codomain(h);

    /* Add corresponding data */
    X_red := ReduceCurveSplit(X, h); Y_red := ReduceCurveSplit(Y, h);
    NormM_red := ReduceMatrixSplit(NormM, h);

    Iterator_red := InitializedIterator(X_red, Y_red, NormM_red, X`g + 3);
    while true do
        found, S_red, Iterator_red := DivisorFromMatrixByDegree(X_red, Y_red, Iterator_red, d : Margin := Margin);
        /* If that does not work, give up and try one degree higher. Note that
         * d is initialized in the outer loop, so that we keep the degree that
         * works. */
        if found then
            break;
        end if;
        if d ge UpperBound then
            return false, [];
        end if;
        d := Max(2*d, UpperBound);
    end while;
    Append(~DEss_red, DefiningEquations(S_red));

    vprint EndoCheck : "";
    vprint EndoCheck : "Fractional CRT...";
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
            DE +:= FractionalCRTSplit(rs, prs : I := I) * Monomial(Rprod, exp);
        end for;
        Append(~DEs, DE);
    end for;
    vprint EndoCheck : "done.";

    vprint EndoCheck : "";
    vprint EndoCheck : "Checking:";
    vprint EndoCheck : "Step 1... ";
    /* Note that P and Qs are calculated at the beginning of this function */
    test1 := CheckEquations(X, Y, P, Qs, DEs);
    vprint EndoCheck : "done.";

    if test1 then
        vprint EndoCheck : "Step 2...";
        S := Scheme(AffineSpace(Rprod), DEs);
        test2 := CheckIrreducibleComponent(X, Y, S);
        vprint EndoCheck : "done.";
        if test2 then
            vprint EndoCheck : "";
            vprint EndoCheck : "Divisor found!";
            return true, S;
        end if;
    end if;
end while;

end intrinsic;


function DivisorFromMatrixByDegree(X, Y, Iterator, d : Margin := 2^5)

vprint EndoCheck, 2 : "";
vprint EndoCheck, 2 : "Trying degree", d;
fs := CandidateDivisors(X, Y, d);
n := #fs + Margin;
vprint EndoCheck, 2 : "Number of terms in expansion:", n;

vprint EndoCheck, 2 : "Expanding branches...";
while true do
    P, Qs, _, _ := Explode(Iterator); Pold := P; Qsold := Qs;
    prec := Precision(Parent(Qs[1][1]));
    if prec ge n then
        break;
    end if;
    Iterator := IterateIterator(Iterator);
    P, Qs, _, _ := Explode(Iterator);
    /* Check that iteration does what it is supposed to do */
    assert &and[ IsWeaklyZero(P[j] - Pold[j]) : j in [1..#P] ];
    assert &and[ &and[ IsWeaklyZero(Qs[i][j] - Qsold[i][j]) : j in [1..#Qs[i] ] ] : i in [1..#Qs] ];
end while;
vprint EndoCheck, 2 : "done.";

/* Fit a divisor to it */
vprint EndoCheck, 2 : "Solving linear system...";
vtime EndoCheck, 2:
ICs := IrreducibleComponentsFromBranches(X, Y, fs, P, Qs : Margin := Margin div 2);
vprint EndoCheck, 2 : "done.";

for S in ICs do
    DEs := DefiningEquations(S);
    vprint EndoCheck, 2 : "Checking:";
    vprint EndoCheck, 2 : "Step 1...";
    vtime EndoCheck, 2:
    test1 := CheckEquations(X, Y, P, Qs, DEs);
    vprint EndoCheck, 2 : "done.";
    if test1 then
        vprint EndoCheck, 2 : "Step 2...";
        vtime EndoCheck, 2:
        test2 := CheckIrreducibleComponent(X, Y, S);
        vprint EndoCheck, 2 : "done.";
        if test2 then
            vprint EndoCheck, 2 : "Divisor found!";
            return true, S, Iterator;
        end if;
    end if;
end for;
return false, [ ], Iterator;

end function;
