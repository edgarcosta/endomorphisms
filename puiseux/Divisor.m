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
declare attributes Crv : g, U, P0, A, DEs;
declare attributes Crv : patch_index, unif_index;
declare attributes Crv : R, K;
declare attributes Crv : F, rF, OF, BOF;
declare attributes Crv : OurB, NormB, T;
declare attributes Crv : initialized;
declare attributes Crv : cantor_equations;

declare verbose EndoCheck, 4;


forward InitializeCurve;
forward OurAffinePatch;
forward AlgebraicUniformizerIndex;
forward OurBasisOfDifferentials;
forward ChangeTangentAction;
forward NormalizedBasisOfDifferentials;

forward VariableOrder;
forward ExtractPoints;
forward ExtractHomomorphismsRing;

forward CandidateDivisors;
forward IrreducibleComponentsFromBranches;
forward CheckEquations;
forward CheckIrreducibleComponent;

forward DivisorFromMatrix;
forward DivisorFromMatrixSplit;
forward DivisorFromMatrixByDegree;

forward ExtractHomomorphismsField;
forward CandidateDivisorsNew;
forward IrreducibleComponentsFromBranchesNew;
forward DivisorFromMatrixByDegreeNew;


import "LocalInfo.m": DevelopPoint, InitializeImageBranch;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "Cantor.m": CantorEquations;


procedure InitializeCurve(X, P0)
/* This is written in a silly way; it is essentially a procedure */

if not assigned X`initialized then
    X`initialized := false;
end if;
if X`initialized then
    return;
end if;

vprintf EndoCheck, 3 : "Curve:\n";
vprint EndoCheck, 3 : X;
X`is_hyperelliptic := IsHyperelliptic(X); X`is_planar := IsPlaneCurve(X); X`is_smooth := IsNonSingular(X);
X`g := Genus(X); X`is_plane_quartic := (X`is_planar) and (X`is_smooth) and (X`g eq 3);
if not X`is_planar then
    error "Please give your curve in planar form";
end if;

/* Find affine patch: */
X`U, X`DEs, X`P0, X`patch_index := OurAffinePatch(X, P0);
X`A := Ambient(X`U); X`R := CoordinateRing(X`A); X`K := FieldOfFractions(X`R);

/* Put uniformizer first: */
X`unif_index := AlgebraicUniformizerIndex(X);
if X`unif_index eq 2 then
    X`DEs := [ X`R ! Evaluate(DE, [ (X`R).2, (X`R).1 ]) : DE in X`DEs ];
    X`U := Scheme(X`A, X`DEs);
    X`P0 := X`U ! [ X`P0[2], X`P0[1] ];
end if;

vprintf EndoCheck, 3 : "Index of affine patch: ";
vprint EndoCheck, 3 : X`patch_index;
vprintf EndoCheck, 3 : "Index of uniformizer: ";
vprint EndoCheck, 3 : X`unif_index;
vprintf EndoCheck, 3 : "Affine patch:\n";
vprint EndoCheck, 3 : X`U;
vprintf EndoCheck, 3 : "Point:\n";
vprint EndoCheck, 3 : X`P0;

/* Construct equation order */
X`F := BaseRing(X`R);
if Type(X`F) eq FldRat then
    X`rF := 1;
    X`OF := Integers();
    X`BOF := Basis(X`OF);
elif Type(X`F) eq FldNum then
    X`rF := Denominator(X`F.1) * X`F.1;
    X`OF := Order([ X`rF^i : i in [0..Degree(X`F) - 1] ]);
    X`BOF := Basis(X`OF);
end if;

X`OurB := OurBasisOfDifferentials(X);
X`NormB, X`T := NormalizedBasisOfDifferentials(X);
_<u,v> := Parent(X`OurB[1]);
vprintf EndoCheck, 3 : "Standard basis of differentials:\n";
vprint EndoCheck, 3 : X`OurB;
vprintf EndoCheck, 3 : "Normalized basis of differentials:\n";
vprint EndoCheck, 3 : X`NormB;

X`cantor_equations := CantorEquations(X);
vprintf EndoCheck, 3 : "Cantor equations:\n";
vprint EndoCheck, 3 : X`cantor_equations;
X`initialized := true;
return;

end procedure;


function OurAffinePatch(X, P0);

if IsAffine(X) then
    return X, DefiningEquations(X), P0, 1;

elif X`is_hyperelliptic or (X`g eq 1) then
    /* Probably this does nothing... but still */
    U, P0, patch_index := AffinePatch(X, P0);
    DEs := DefiningEquations(AffinePatch(X, 1));
    R := Parent(DEs[1]);
    d := 2*(X`g) + 2;
    if patch_index eq 3 then
        DEs := [ R ! (-R.2^d * Evaluate(DE, [ 1/R.2, R.1/(R.2^(d div 2)) ])) : DE in DEs ];
    end if;
    U := Curve(AffineSpace(R), DEs);
    return U, DEs, P0, patch_index;

elif X`is_planar then
    U, P0, patch_index := AffinePatch(X, P0);
    DEs := DefiningEquations(AffinePatch(X, 1));
    R := Parent(DEs[1]);
    d := Degree(DEs[1]);
    if patch_index eq 2 then
        DEs := [ R ! (R.2^d * Evaluate(DE, [ R.1/R.2, 1/R.2 ])) : DE in DEs ];
    elif patch_index eq 3 then
        DEs := [ R ! (R.2^d * Evaluate(DE, [ 1/R.2, R.1/R.2 ])) : DE in DEs ];
    end if;
    U := Curve(AffineSpace(R), DEs);
    return U, DEs, P0, patch_index;
end if;

end function;


function AlgebraicUniformizerIndex(X)
/*
 * Input:   A plane curve X.
 * Output:  Index of the uniformizer.
 */

if X`g eq 1 then
    fX := X`DEs[1]; R := X`R; P0 := X`P0;
    // Prefer coordinate on PP^1:
    /* NOTE: Do NOT neglect to take an Eltseq here; omitting it is deadly,
     * since evaluating x at (0, 0) can be 0 */
    if Degree(fX, R.2) eq 2 then
        if Evaluate(Derivative(fX, R.2), Eltseq(P0)) ne 0 then
            return 1;
        else
            return 2;
        end if;
    else
        if Evaluate(Derivative(fX, R.1), Eltseq(P0)) ne 0 then
            return 2;
        else
            return 1;
        end if;
    end if;

elif X`is_hyperelliptic then
    // In this case we always get the coordinate on PP^1, since we avoid
    // Weierstrass points.
    if X`patch_index eq 1 then
        return 1;
    else
        return 2;
    end if;

else
    // Here we do the usual test, without the preference of the elliptic case.
    fX := X`DEs[1]; R := X`R; P0 := X`P0;
    /* NOTE: Do NOT neglect to take an Eltseq here; omitting it is deadly,
     * since evaluating x at (0, 0) can be 0 */
    if Evaluate(Derivative(fX, R.2), Eltseq(P0)) ne 0 then
        return 1;
    else
        return 2;
    end if;
end if;

end function;


function OurBasisOfDifferentials(X)
/*
 * Input:   A curve X.
 * Output:  A basis of global differentials on X, represented by elements of
 *          the rational function field by using our choice of uniformizer
 */

g := X`g; R := X`R; u := X`R.1; v := X`R.2; f := X`DEs[1];
if g eq 0 then
    return [ ];

elif g eq 1 then
    /* Elliptic case: we use dx / 2y */
    if Degree(f, v) eq 2 then
        s := MonomialCoefficient(f, v^2);
        return [ s / Derivative(f, v) ];
    else
        s := MonomialCoefficient(f, u^2);
        return [ -s / Derivative(f, v) ];
    end if;

elif X`is_hyperelliptic then
    /* (Hyper)elliptic case: we use x^i dx / 2y */
    s := MonomialCoefficient(f, v^2);
    return [ s*u^(i-1) / Derivative(f, v) : i in [1..g] ];

elif X`is_plane_quartic then
    /* Plane quartic case: we use ({x,y,1} / (dF / dy)) dx */
    return [ X`K ! ( n / Derivative(f, v)) : n in [u, v, 1] ];

else
    error "OurBasisOfDifferentials not implemented yet for this curve";
end if;

end function;


function ChangeTangentAction(X, Y, M)
/*
 * Input:  Two curves X, Y and a representation M on the standard basis of
 *         differentials.
 * Output: Matrix for standard differentials on the patches of X and Y used.
 */

F := X`F;
/* M acts on the right, so to precompose with the operation on X we multiply on
 * the right; we modify the columns. */
M := Transpose(M);
if X`g eq 1 or X`is_hyperelliptic then
    if X`patch_index eq 3 then
        vprint EndoCheck, 3 : "Modifying tangent action for patch index of X";
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, Reverse(rows));
    end if;

elif X`is_plane_quartic then
    if X`patch_index eq 2 then
        vprint EndoCheck, 3 : "Modifying tangent action for patch index of X";
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, [ rows[1], rows[3], rows[2] ]);
    elif X`patch_index eq 3 then
        vprint EndoCheck, 3 : "Modifying tangent action for patch index of X";
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := Matrix(F, [ rows[2], rows[3], rows[1] ]);
    end if;
    if X`unif_index eq 2 then
        vprint EndoCheck, 3 : "Modifying tangent action for uniformizing index of X";
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, [ rows[2], rows[1], rows[3] ]);
    end if;
end if;
M := Transpose(M);

/* For y we postcompose, hence we have to modify columns; we therefore take a
 * transpose and go back */
if Y`g eq 1 or Y`is_hyperelliptic then
    if Y`patch_index eq 3 then
        vprint EndoCheck, 3 : "Modifying tangent action for uniformizing index of Y";
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, Reverse(rows));
    end if;

elif Y`is_plane_quartic then
    if Y`patch_index eq 2 then
        vprint EndoCheck, 3 : "Modifying tangent action for patch index of Y";
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, [ rows[1], rows[3], rows[2] ]);
    elif Y`patch_index eq 3 then
        vprint EndoCheck, 3 : "Modifying tangent action for patch index of Y";
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := Matrix(F, [ rows[2], rows[3], rows[1] ]);
    end if;
    if Y`unif_index eq 2 then
        vprint EndoCheck, 3 : "Modifying tangent action for uniformizing index of Y";
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, [ rows[2], rows[1], rows[3] ]);
    end if;
end if;

return M;

end function;


function NormalizedBasisOfDifferentials(X)
/*
 * Input:   A curve X.
 * Output:  A differential basis Bnorm that is normalized with respect to the uniformizing parameter,
 *          and a matrix T such that multiplication by T on the left sends B to Bnorm.
 */

P := DevelopPoint(X, X`P0, X`g);
BP := [ Evaluate(b, P) : b in X`OurB ];
T := Matrix([ [ Coefficient(BP[i], j - 1) : j in [1..X`g] ] : i in [1..X`g] ])^(-1);
NormB := [ &+[ T[i,j] * X`OurB[j] : j in [1..X`g] ] : i in [1..X`g] ];
return NormB, T;

end function;


function VariableOrder()
/*
 * The order in which x(P), y(P), x(Q), y(Q) and hence x1, y1, x2, y2 are used
 * in the product space. Note that because of the lexicographical ordering
 * variables that occur later are eliminated for first.
 */

/* x(P) to 4th comp, y(P) to 2nd comp, etc */
// TODO: Test better ones
//return [1, 2, 3, 4];
return [4, 2, 3, 1];

end function;


function ExtractPoints(X, Y, P, Q)
/* Reflects order in VariableOrder */

seq := [ P[1], P[2], Q[1], Q[2] ];
varord := VariableOrder();
return [ seq[varord[i]] : i in [1..#varord] ];

end function;


function ExtractHomomorphismsRing(X, Y)

RX := X`R; RY := Y`R;
varord := VariableOrder();
// TODO: Test other orderings
Rprod := PolynomialRing(X`F, 4, "lex");
seqX := [ Rprod.Index(varord, i) : i in [1..2] ];
seqY := [ Rprod.Index(varord, i) : i in [3..4] ];
hX := hom<RX -> Rprod | seqX >;
hY := hom<RY -> Rprod | seqY >;
return [ hX, hY ];

end function;


function CandidateDivisors(X, Y, d)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

gX := X`g; fX := X`DEs[1]; RX := X`R;
xX := RX.1; yX := RX.2;
/* Change in hyperelliptic case for greater effectiveness: */
if Degree(fX, RX.1) eq 2 then
    xX := RX.2; yX := RX.1;
end if;
gY := Y`g; fY := Y`DEs[1]; RY := Y`R;
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

hs := ExtractHomomorphismsRing(X, Y);
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
hX, hY := Explode(ExtractHomomorphismsRing(X, Y));
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


function ExtractHomomorphismsField(X, Y)

KX := X`K; KY := Y`K;
varord := VariableOrder();
// TODO: Test other orderings
Kprod := RationalFunctionField(X`F, 4);
seqX := [ Kprod.Index(varord, i) : i in [1..2] ];
seqY := [ Kprod.Index(varord, i) : i in [3..4] ];
hX := hom<KX -> Kprod | seqX >;
hY := hom<KY -> Kprod | seqY >;
return [ hX, hY ];

end function;


function CandidateDivisorsNew(X, Y, d)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

gX := X`g; fX := X`DEs[1]; RX := X`R; KX := X`K; P := X`P0;
gY := Y`g; fY := Y`DEs[1]; RY := Y`R; KY := Y`K; Q := Y`P0;
V, phiV := RiemannRochSpace((d + 2*gX)*Divisor(P));
V, phiV := RiemannRochSpace((d + gX)*Divisor(P));
divsX := [ KX ! phiV(v) : v in Basis(V) ];
W, phiW := RiemannRochSpace(3*gY*Divisor(Q));
W, phiW := RiemannRochSpace(2*gY*Divisor(Q));
divsY := [ KY ! phiW(w) : w in Basis(W) ];
divsX := [ divsX[#divsX] ] cat divsX[1..(#divsX - 1)];
divsY := [ divsY[#divsY] ] cat divsY[1..(#divsY - 1)];

hs := ExtractHomomorphismsField(X, Y);
CP := [ [* divX, divY *] : divX in divsX, divY in divsY ];
divs := [ &*[ hs[i](tup[i]) : i in [1..2] ] : tup in CP ];
return divs;

end function;


function IrreducibleComponentsFromBranchesNew(X, Y, fs, P, Qs : DivPP1 := false)
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
evss := [ [ Evaluate(f, ExtractPoints(X, Y, P, Q)) : Q in Qs ] : f in fs ];
min := Minimum([ Valuation(ev) : ev in &cat(evss) ]);
max := Minimum([ AbsolutePrecision(ev) : ev in &cat(evss) ]);
M := Matrix([ &cat[ [ Coefficient(ev, i/e) : i in [(e*min)..(e*max - 10)] ] : ev in evs ] : evs in evss ]);
min := -9;
max := 300/2;
M := Matrix([ &cat[ [ Coefficient(ev, i/e) : i in [(e*min)..(e*max)] ] : ev in evs ] : evs in evss ]);
B := Basis(Kernel(M));
/* Coerce back to ground field (possible because of echelon form) */
B := [ [ X`F ! c : c in Eltseq(b) ] : b in B ];
print "Dimension:";
print #B;
return 0;

/* Corresponding equations */
hX, hY := Explode(ExtractHomomorphismsRing(X, Y));
Rprod := Codomain(hX);
Aprod := AffineSpace(Rprod);
eqs := [ &+[ b[i] * fs[i] : i in [1..#fs] ] : b in B ];
eqs := [ Rprod ! Numerator(&+[ b[i] * fs[i] : i in [1..#fs] ]) : b in B ];
S := Scheme(Aprod, [ hX(DE) : DE in X`DEs ] cat [ hY(DE) : DE in Y`DEs ] cat eqs);

I := DefiningIdeal(S);
varord := VariableOrder();
J1 := ideal<Rprod | [ hX(DE) : DE in X`DEs ] cat [ Rprod.Index(varord, i + 2) - (X`P0)[i] : i in [1..2] ]>;
J2 := ideal<Rprod | [ hY(DE) : DE in Y`DEs ] cat [ Rprod.Index(varord, i) - (Y`P0)[i] : i in [1..2] ]>;
vprintf EndoCheck, 3 : "Calculating colon ideals... ";
repeat
    Iold := I;
    I := ColonIdeal(I, J1);
until I eq Iold;
repeat
    Iold := I;
    I := ColonIdeal(I, J2);
until I eq Iold;
vprintf EndoCheck, 3 : "done.\n";
S := Scheme(Aprod, I);
S := ReducedSubscheme(S);
vprint EndoCheck, 3 : "Dimension:";
vprint EndoCheck, 3 : Dimension(S);
vprint EndoCheck, 3 : "Degree and dimensions to factors:";
vprint EndoCheck, 3 : BiDimDeg(X, X, S);
/*
for I in IrreducibleComponents(S) do
    print BiDimDeg(X, X, I);
end for;
*/

/*
if DivPP1 then
*/

/* Corresponding scheme */
eqs := [ &+[ b[i] * fs[i] : i in [1..#fs] ] : b in B ];
Kprod := FieldOfFractions(Rprod);
Sprod := Scheme(Aprod, [ hX(DE) : DE in X`DEs ] cat [ hY(DE) : DE in Y`DEs ]);
Ss := [ ];
for f in eqs do
    Scl := ProjectiveClosure(Sprod);
    Rcl := CoordinateRing(Ambient(Scl));
    Ccl := FieldOfFractions(Rcl);
    h := hom< Kprod -> Ccl | [ Ccl.1, Ccl.2, Ccl.3, Ccl.4 ] >;
    h := hom< Rprod -> Rcl | [ Rcl.1, Rcl.2, Rcl.3, Rcl.4 ] >;
    S_num := Divisor(Scl, h(Numerator(f)));
    S_den := Divisor(Scl, h(Denominator(f)));
    S := S_num - S_den;
    Append(~Ss, S);
end for;
S := GCD(Ss);
vprint EndoCheck, 3 : Support(SignDecomposition(S));
vprint EndoCheck, 3 : "Dimension:";
vprint EndoCheck, 3 : Dimension(S);
vprint EndoCheck, 3 : BiDimDeg(X, X, S);

return 0;
return [ S ];

end function;


function DivisorFromMatrixByDegreeNew(X, Y, NormM, d : Margin := 2^4, DivPP1 := false, have_to_check := true)

vprintf EndoCheck, 2 : "Trying degree %o...\n", d;
fs := CandidateDivisorsNew(X, Y, d);
n := #fs + Margin;
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
ICs := IrreducibleComponentsFromBranchesNew(X, Y, fs, P, Qs : DivPP1 := DivPP1);
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


intrinsic BiDimDeg(X::., Y::., I::.) -> .
{Gives dimension and degree of both projections.}

A4 := Ambient(I); R4 := CoordinateRing(A4);
R2 := PolynomialRing(X`F, 2); A2 := AffineSpace(R2);
seq1 := [ X`P0[1], X`P0[2], R2.1, R2.2 ];
seq2 := [ R2.1, R2.2, X`P0[1], X`P0[2] ];
varord := VariableOrder();
seq1 := [ seq1[varord[i]] : i in [1..#varord] ];
seq2 := [ seq2[varord[i]] : i in [1..#varord] ];
h1 := hom< R4 -> R2 | seq1 >;
h2 := hom< R4 -> R2 | seq2 >;
eqs1 := [ h1(eq4) : eq4 in DefiningEquations(I) ];
eqs2 := [ h2(eq4) : eq4 in DefiningEquations(I) ];
S1 := Scheme(A2, eqs1);
S2 := Scheme(A2, eqs2);
return [ [ Dimension(S1), Dimension(S2) ], [ Degree(S1), Degree(S2) ] ];

end intrinsic;


intrinsic Bidegree(X::., Y::., I::.) -> .
{Gives degree of both projections.}

A4 := Ambient(I); R4 := CoordinateRing(A4);
R2 := PolynomialRing(X`F, 2); A2 := AffineSpace(R2);
seq1 := [ X`P0[1], X`P0[2], R2.1, R2.2 ];
seq2 := [ R2.1, R2.2, X`P0[1], X`P0[2] ];
varord := VariableOrder();
seq1 := [ seq1[varord[i]] : i in [1..#varord] ];
seq2 := [ seq2[varord[i]] : i in [1..#varord] ];
h1 := hom< R4 -> R2 | seq1 >;
h2 := hom< R4 -> R2 | seq2 >;
eqs1 := [ h1(eq4) : eq4 in DefiningEquations(I) ];
eqs2 := [ h2(eq4) : eq4 in DefiningEquations(I) ];
S1 := Scheme(A2, eqs1);
S2 := Scheme(A2, eqs2);
return [ Degree(S1), Degree(S2) ];

f1 := map< Curve(I) -> Curve(X`U) | [ R4.varord[1], R4.varord[2] ]>;
f2 := map< Curve(I) -> Curve(X`U) | [ R4.varord[3], R4.varord[4] ]>;
f1 := ProjectiveClosure(f1); f2 := ProjectiveClosure(f2);
return [ Degree(f1), Degree(f2) ];

end intrinsic;
