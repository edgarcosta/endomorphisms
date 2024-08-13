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


declare attributes Crv : is_hyperelliptic, is_planar, is_smooth, is_plane_quartic;
declare attributes Crv : g, U, P0, A, DEs, DEs_sub;
declare attributes Crv : patch_index, unif_index;
declare attributes Crv : RA, KA, RU, KU;
declare attributes Crv : F, rF, OF, BOF;
declare attributes Crv : OurB, NormB, echelon_exps, T;
declare attributes Crv : cantor_eqs;
declare attributes Crv : RRgens, globgens;
declare attributes Crv : initialized;

declare verbose EndoCheck, 4;


forward OurAffinePatch;
forward AlgebraicUniformizerIndex;
forward OurBasisOfDifferentials;
forward NormalizedBasisOfDifferentials;
forward CantorEquations;
forward ChangeTangentAction;

forward InitializeCurve;


function OurAffinePatch(X, P0);
/* The patches are chosen in such a way that we can later use the same differentials */
/* NOTE: Too much here happens by hand... but the minus signs are arguably important */

if IsAffine(X) then
    return X, DefiningEquations(X), P0, 1;

elif X`is_hyperelliptic or (X`g eq 1) then
    U, P0, patch_index := AffinePatch(X, P0); P0 := Eltseq(P0);
    DEs := DefiningEquations(AffinePatch(X, 1));
    RA := Parent(DEs[1]);
    d := 2*(X`g) + 2;
    if patch_index eq 3 then
        /* Make sure that second variable defines hyperelliptic extension with
         * same sign as usual */
        DEs := [ RA ! (RA.1^d * Evaluate(DE, [ 1/RA.1, RA.2/(RA.1^(d div 2)) ])) : DE in DEs ];
        P0 := Reverse(P0);
    end if;
    /* Check correctness of above */
    //print U; print DEs;
    U := Curve(AffineSpace(RA), DEs);
    return U, DEs, U ! P0, patch_index;

elif X`is_planar then
    U, P0, patch_index := AffinePatch(X, P0); P0 := Eltseq(P0);
    DEs := DefiningEquations(AffinePatch(X, 1));
    RA := Parent(DEs[1]);
    if patch_index eq 2 then
        DEs := [ RA ! (RA.2^Degree(DE) * Evaluate(DE, [ RA.1/RA.2, 1/RA.2 ])) : DE in DEs ];
    elif patch_index eq 3 then
        DEs := [ RA ! (RA.2^Degree(DE) * Evaluate(DE, [ 1/RA.2, RA.1/RA.2 ])) : DE in DEs ];
    end if;
    /* Check correctness of above */
    //print U; print DEs;
    U := Curve(AffineSpace(RA), DEs);
    return U, DEs, U ! P0, patch_index;
end if;

end function;


function AlgebraicUniformizerIndex(X)
/*
 * Input:   A plane curve X.
 * Output:  Index of the uniformizer.
 */

if X`is_hyperelliptic or (X`g eq 1) then
    fX := X`DEs[1]; RA := X`RA; P0 := X`P0;
    // Prefer coordinate on PP^1:
    /* NOTE: Do NOT neglect to take an Eltseq here; omitting it is deadly,
     * since evaluating x at (0, 0) can be unequal to 0 */
    if Degree(fX, RA.2) eq 2 then
        if Evaluate(Derivative(fX, RA.2), Eltseq(P0)) ne 0 then return 1; else return 2; end if;
    else
        if Evaluate(Derivative(fX, RA.1), Eltseq(P0)) ne 0 then return 2; else return 1; end if;
    end if;

else
    // Here we do the usual test, without the preference of the hyperelliptic case.
    fX := X`DEs[1]; RA := X`RA; P0 := X`P0;
    /* NOTE: Do NOT neglect to take an Eltseq here; omitting it is deadly,
     * since evaluating x at (0, 0) can be unequal to 0 */
    if Evaluate(Derivative(fX, RA.2), Eltseq(P0)) ne 0 then return 1; else return 2; end if;
end if;

end function;


function OurBasisOfDifferentials(X)
/*
 * Input:   A curve X.
 * Output:  A basis of global differentials on X, represented by elements of
 *          the rational function field by using our choice of uniformizer,
 *          which at the point of invocation is the first coordinate.
 */
/* NOTE: This function is only invoked AFTER the curve is in its correct patch */

g := X`g; RA := X`RA; u := RA.1; v := RA.2; f := X`DEs[1];
if g eq 0 then
    return [ ];

elif (g eq 1) or X`is_hyperelliptic then
    /* Elliptic case: we use dx / 2y */
    /* Usual version where the uniformizing coordinate u is that on PP^1 */
    if Degree(f, v) eq 2 then
        s := MonomialCoefficient(f, v^2);
        return [  s*u^(i-1) / Derivative(f, v) : i in [1..g] ];
    else
        /* If coordinate on PP^1 did not work, then this is the expression in
         * the new uniformizer. This only happens in the elliptic case, meaning
         * that we do not have to adapt differentials then */
        s := MonomialCoefficient(f, u^2);
        return [ -s*v^(i-1) / Derivative(f, v) : i in [1..g] ];
    end if;

elif X`is_plane_quartic then
    /* Plane quartic case: we use ({x,y,1} / (dF / dy)) dx */
    return [ X`KA ! (n / Derivative(f, v)) : n in [u, v, 1] ];

else
    x := X`KU ! (X`RU).1;
    return [ X`KA ! (b / Differential(x)) : b in BasisOfDifferentialsFirstKind(X`U) ];
end if;

end function;


function NormalizedBasisOfDifferentials(X : AssertNonWP := false)
/*
 * Input:   A curve X.
 * Output:  A differential basis Bnorm that is normalized with respect to the uniformizing parameter,
 *          and a matrix T such that multiplication by T on the left sends B to Bnorm.
 */

P := DevelopPoint(X, X`P0, 2*X`g + 1);
BP := [ Evaluate(b, P) : b in X`OurB ];
M := Matrix([ [ Coefficient(BP[i], j - 1) : j in [1..(2*X`g + 1)] ] : i in [1..X`g] ]);
E, T := EchelonForm(M); c := #Rows(Transpose(E));
echelon_exps := [ Minimum([ i : i in [1..c] | Eltseq(row)[i] ne 0 ]) : row in Rows(E) ];
if AssertNonWP and (not echelon_exps eq [1..X`g]) then
    error "Base point is Weierstrass";
end if;
NormB := [ &+[ T[i,j] * X`OurB[j] : j in [1..X`g] ] : i in [1..X`g] ];
return NormB, T, echelon_exps;

end function;


function CantorEquations(X);
/* Gives the equations in the description a (x) = 0, y = b (x)
 * May not use usual parameter if uniformizer differs */
/* TODO: That last convention is a bit suboptimal */

g := X`g; f := X`DEs[1]; F := X`F;
S := PolynomialRing(F, 2*g); T<t> := PolynomialRing(S);
/* Names:
 * a1 is trace term before t^(g - 1), a_g is norm term before t^0,
 * b1 is term before t^(g - 1), bg is term before t^0 */
varnames := [ Sprintf("a%o", i) : i in [1..g] ] cat [ Sprintf("b%o", i) : i in [1..g] ];
AssignNames(~S, varnames);

/* Start with trace and end with norm */
canpol := t^g + &+[ S.i * t^(g - i) : i in [1..g] ];
substpol := &+[ S.(g + i) * t^(g - i) : i in [1..g] ];
P := [t, substpol];
eqpol := Evaluate(f, P) mod canpol;
return Coefficients(eqpol);

end function;


procedure InitializeCurve(X, P0 : AssertNonWP := false)

if not assigned X`initialized then X`initialized := false; end if;
if X`initialized then return; end if;
if AssertNonWP then assert not IsWeierstrassPlace(Place(X ! P0)); end if;

X`is_hyperelliptic := Type(X) eq CrvHyp; X`is_planar := IsPlaneCurve(X); X`is_smooth := IsNonSingular(X);
X`g := Genus(X); X`is_plane_quartic := (X`is_planar) and (X`is_smooth) and (X`g eq 3);
if not X`is_planar then error "Please give your curve in hyperelliptic or planar form"; end if;

/* Find affine patch, then put uniformizer first by switching coordinates */
X`U, X`DEs, X`P0, X`patch_index := OurAffinePatch(X, P0);
X`A := Ambient(X`U); RA<u,v> := CoordinateRing(X`A); X`RA := RA; X`KA := FieldOfFractions(RA);
X`unif_index := AlgebraicUniformizerIndex(X);
if X`unif_index eq 2 then
    X`DEs := [ X`RA ! Evaluate(DE, [ (X`RA).2, (X`RA).1 ]) : DE in X`DEs ];
    X`U := Curve(Scheme(X`A, X`DEs));
    X`P0 := X`U ! [ X`P0[2], X`P0[1] ];
end if;
X`A := Ambient(X`U); RA<u,v> := CoordinateRing(X`A); X`RA := RA; X`KA := FieldOfFractions(RA);
X`RU := CoordinateRing(X`U); X`KU := FunctionField(X`U); X`F := BaseRing(X`RU);

vprint EndoCheck, 3 : "";
vprint EndoCheck, 3 : "Index of affine patch:";
vprint EndoCheck, 3 : X`patch_index;
vprint EndoCheck, 3 : "Index of uniformizer:";
vprint EndoCheck, 3 : X`unif_index;
vprint EndoCheck, 3 : "Affine patch:";
vprint EndoCheck, 3 : X`U;
vprint EndoCheck, 3 : "Base point:";
vprint EndoCheck, 3 : X`P0;

/* Construct equation order */
if Type(X`F) eq FldRat then
    X`OF := Integers();
elif Type(X`F) eq FldNum then
    X`OF := EquationOrder(X`F);
end if;

X`OurB := OurBasisOfDifferentials(X);
X`NormB, X`T, X`echelon_exps := NormalizedBasisOfDifferentials(X : AssertNonWP := AssertNonWP);
_<u,v> := Parent(X`OurB[1]);
vprint EndoCheck, 3 : "";
vprint EndoCheck, 3 : "Standard basis of differentials:";
vprint EndoCheck, 3 : X`OurB;
vprint EndoCheck, 3 : "Normalized basis of differentials:";
vprint EndoCheck, 3 : X`NormB;

X`cantor_eqs := CantorEquations(X);
vprint EndoCheck, 2 : "";
vprint EndoCheck, 2 : "Cantor equations:";
vprint EndoCheck, 2 : X`cantor_eqs;

X`initialized := true;
return;

end procedure;


function ChangeTangentAction(X, Y, M)
/*
 * Input:  Two curves X, Y and a representation M on the standard basis of
 *         differentials.
 * Output: Matrix for standard differentials on the patches of X and Y used.
 */
/* NOTE: Only invoked latter to correct for basis chosen in OurBasisOfDifferentials.
   TODO: Apply it here instead, if possible. */

F := X`F;
/* M acts on the right, so to precompose with the operation on X we modify the
 * columns. */
M := Transpose(M);
if X`g eq 1 or X`is_hyperelliptic then
    if X`patch_index eq 3 then
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, Reverse(rows));
    end if;

elif X`is_plane_quartic then
    if X`patch_index eq 2 then
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, [ rows[1], rows[3], rows[2] ]);
    elif X`patch_index eq 3 then
        rows := [ Eltseq(row) : row in Rows(M) ];
        M :=  Matrix(F, [ rows[2], rows[3], rows[1] ]);
    end if;
    if X`unif_index eq 2 then
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, [ rows[2], rows[1], rows[3] ]);
    end if;
end if;
M := Transpose(M);

/* M acts on the right, so to postcompose with the operation on Y we modify the
 * rows. */
if Y`g eq 1 or Y`is_hyperelliptic then
    if Y`patch_index eq 3 then
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, Reverse(rows));
    end if;

elif Y`is_plane_quartic then
    if Y`patch_index eq 2 then
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, [ rows[1], rows[3], rows[2] ]);
    elif Y`patch_index eq 3 then
        rows := [ Eltseq(row) : row in Rows(M) ];
        /* This line is essentially a coincidence, since the corresponding
         * inverse happens to coincide with a transpose */
        M :=  Matrix(F, [ rows[2], rows[3], rows[1] ]);
    end if;
    if Y`unif_index eq 2 then
        rows := [ Eltseq(row) : row in Rows(M) ];
        M := -Matrix(F, [ rows[2], rows[1], rows[3] ]);
    end if;
end if;

return M;

end function;
