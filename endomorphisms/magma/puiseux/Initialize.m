/***
 *  Curve initialization and conventions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "Branches.m": DevelopPoint;
import "Cantor.m": CantorEquations;


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


procedure InitializeCurve(X, P0)

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
/* M acts on the right, so to precompose with the operation on X we modify the
 * columns. */
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

/* M acts on the right, so to postcompose with the operation on Y we modify the
 * rows. */
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
