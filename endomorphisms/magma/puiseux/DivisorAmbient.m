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


function IrreducibleComponentsFromBranches(X, Y, fs, P, Qs)
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
M := [ ];
for f in fs do
    r := [ ];
    for Q in Qs do
        seq := ExtractPoints(X, Y, P, Q);
        ev := Evaluate(f, seq);
        /* Small buffer: we might just as well take 1 instead of X`g */
        r cat:= [ Coefficient(ev, i/e) : i in [0..prec - X`g] ];
    end for;
    Append(~M, r);
end for;
B := Basis(Kernel(Matrix(M)));
/* Coerce back to ground field (possible because of echelon form) */
B := [ [ X`F ! c : c in Eltseq(b) ] : b in B ];

/* Corresponding equations */
hX, hY, hxs, hxsinv, hys, hysinv := ExtractHomomorphismsRing(X, Y);
Rprod := Codomain(hX);
eqs := [ Rprod ! (&+[ b[i] * fs[i] : i in [1..#fs] ]) : b in B ];
eqs := eqs cat [ hX(DE) : DE in X`DEs ] cat [ hY(DE) : DE in Y`DEs ];

if X`is_hyperelliptic and Y`is_hyperelliptic then
    Iprod := ideal< Rprod | eqs >;
    varord := VariableOrder();
    Bxs := Basis(EliminationIdeal(Iprod, { varord[1], varord[3] }));
    //Bys := Basis(EliminationIdeal(Iprod, { varord[1], varord[4] }));
    if (#Bxs ne 0) then
        gcdx := hxsinv(GCD([ hxs(b) : b in Bxs ]));
        //gcdy := hysinv(GCD([ hys(b) : b in Bys ]));
        eqs cat:= [ gcdx ];
        //eqs := [ gcdx ] cat Reverse(eqs)[1..2] cat [ eqs[1] ];
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


intrinsic DivisorFromMatrixAmbientGlobal(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^5, LowerBound := 1, UpperBound := Infinity()) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

InitializeCurve(X, P0); InitializeCurve(Y, Q0 : NonWP := true);
NormM := ChangeTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);

d := LowerBound;
Iterator := InitializedIterator(X, Y, NormM, 2*Y`g + 1);
while true do
    found, S, Iterator := DivisorFromMatrixByDegree(X, Y, Iterator, d : Margin := Margin);
    if found then
        return true, S;
    end if;
    d +:= 1;
    if d gt UpperBound then
        return false, [];
    end if;
end while;

end intrinsic;


intrinsic DivisorFromMatrixAmbientSplit(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^5, LowerBound := 1, UpperBound := Infinity(), B := 300) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

/* We start at a suspected estimate and then increase degree until we find an appropriate divisor */
InitializeCurve(X, P0); InitializeCurve(Y, Q0 : NonWP := true);
NormM := ChangeTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);

/* Some global elements needed below */
F := X`F; OF := X`OF;
Rprod := PolynomialRing(X`F, 4, "lex");
/* Bit more global margin just to be sure */
Iterator, f := InitializedIterator(X, Y, NormM, 2*Y`g + 1);
P := Iterator[1]; Qs := Iterator[2];

prs := [ ]; DEss_red := [* *];
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

    Iterator_red := InitializedIterator(X_red, Y_red, NormM_red, 2*Y`g + 1);
    while true do
        found, S_red, Iterator_red := DivisorFromMatrixByDegree(X_red, Y_red, Iterator_red, d : Margin := Margin);
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
            DE +:= FractionalCRTSplit(rs, prs : I := I) * Monomial(Rprod, exp);
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


function DivisorFromMatrixByDegree(X, Y, Iterator, d : Margin := 2^5)

vprintf EndoCheck, 2 : "Trying degree %o...\n", d;
fs := CandidateDivisors(X, Y, d);
n := #fs + Margin;
vprintf EndoCheck, 2 : "Number of terms in expansion: %o.\n", n;

vprintf EndoCheck, 2 : "Expanding branches... ";
while true do
    P, Qs, _, _ := Explode(Iterator);
    prec := Precision(Parent(Qs[1][1]));
    /* TODO: Multiply by ramification index like this? */
    //prec := Minimum([ RelativePrecision(c) : c in P cat &cat(Qs) ]);
    if prec ge n then
        break;
    end if;
    Iterator := IterateIterator(Iterator);
end while;
P, Qs, _, _ := Explode(Iterator);
vprintf EndoCheck, 2 : "done.\n";

/* Fit a divisor to it */
vprintf EndoCheck, 2 : "Solving linear system... ";
ICs := IrreducibleComponentsFromBranches(X, Y, fs, P, Qs);
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
            return true, S, Iterator;
        end if;
    end if;
end for;
return false, [ ], Iterator;

end function;
