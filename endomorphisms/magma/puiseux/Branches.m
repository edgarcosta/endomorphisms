/***
 *  Computation of local expansions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

/* Note that in all of this M is the matrix that represents the action on H^0
 * (X, omega_X)^* by left multiplication. */


forward LiftPuiseuxSeries;

forward InitializeMatrix;
forward PuiseuxLeadingExponent;
forward InitializeImageBranch;

forward DevelopPoint;

forward InitializeLift;
forward CreateLiftIteratorFunction;
forward InitializedIterator;
forward IterateIterator;


/* NOTE: Every generator of a Puiseux series ring should be coerced because of
 *       infinite precision issues */

function LiftPuiseuxSeries(f, PR, e)
/*
 * Input:  A Puiseux series f, a Puiseux ring PR, and an exponent denominator e.
 * Output: The coercion of f to PR, if all exponents of f have denominator e.
 */

/* Absolute precision for f is value of N in O (t^N) in expression for f */
L := [ Coefficient(f, i/e)*PR.1^(i/e) : i in [0..e*AbsolutePrecision(f) - 1] ];
if #L eq 0 then
    return PR ! 0;
end if;
return PR ! (&+L);

end function;


function PuiseuxLeadingExponent(M, echelon_exps)
/*
 * Input:   after normalizing to an upper triangular form,
 *          and the echelon exponents for the differentials of the
 *          corresponding curve.
 * Output:  The ramification index of the corresponding Puiseux expansions and
 *          the corresponding entries of M.
 */
/* NOTE: Key equality: if t = s^e, then t^i dt = s^(e i + e - 1) ds up to multiple */

rowsM := Rows(M);
gY := #rowsM; gX := #Eltseq(rowsM[1]);

/* In what follows, we that the point on Y is not Weierstrass.
 * i-th entry gives t^(i - 1) dt, pulling back to s^(e i - 1) ds which has to
 * equal x[j0] - 1. This implies e = x[j0]/i. */
exps := [ ];
for i in [1..gY] do
    row := Eltseq(rowsM[i]);
    j0 := Minimum([ j : j in [1..gX] | row[j] ne 0 ]);
    exp := echelon_exps[j0]/i;
    Append(~exps, exp);
end for;

/* Now reverse the above by seeing if there is a j, necessarily the smallest,
 * for which x[j] = e*i. This is stupid but who cares, the dimension is small. */
exp0 := Minimum(exps);
vals := [ ];
for i in [1..gY] do
    row := Eltseq(rowsM[i]);
    j0 := Minimum([ j : j in [1..gX] | row[j] ne 0 ]);
    if echelon_exps[j0] eq exp0*i then
        Append(~vals, (1/exp0)*M[i, j0]);
    else
        Append(~vals, BaseRing(Parent(M)) ! 0);
    end if;
end for;
return exp0, vals;

end function;


function InitializeImageBranch(M, echelon_exps)
/*
 * Input:   A matrix M that represents an endomorphism,
 *          after normalizing to an upper triangular form,
 *          and the echelon exponents for the differentials of the
 *          corresponding curve.
 * Output:  The leading coefficients of the corresponding Puiseux expansions
 *          and the polynomial that gives their field of definition.
 */
/* NOTE: This needs genericity assumptions */

/* Recovering old invariants: */
F := Parent(M[1,1]);
gY := #Rows(M);
exp, vals := PuiseuxLeadingExponent(M, echelon_exps);

/* If Y has genus 1 then we know more about the start of the development since
 * there is no need to desymmetrize */
if gY eq 1 then
    r := Eltseq(Rows(M)[1]);
    PF := PuiseuxSeriesRing(F, #r + 1); tF := PF.1;
    RF := PolynomialRing(F); xF := RF.1;
    return [ &+[ (r[n] / n) * tF^n : n in [1..#r] ] + O(tF^(#r + 1)) ], xF - 1;
end if;

/* Normalized equations (depend only on the matrix): */
A := AffineSpace(F, gY); RA := CoordinateRing(A);
eqs := [ ];
for n in [1..gY] do
    /* The next line implictly uses that Y has echelon exponents [1..g];
     * it inspects the coefficients in front of the pullbacks 
     * t^i dt = s^(e i - * 1) ds when using t_i = c_i s^e. */
    powersum := &+[ RA.i^n : i in [1..gY] ];
    Append(~eqs, powersum - vals[n]);
end for;
S := Scheme(A, eqs);

G := GroebnerBasis(ideal<RA | eqs>);
RF := PolynomialRing(F);
hc := [ RF ! 0 : i in [1..gY] ]; hc[#hc] := RF.1; h := hom<RA -> RF | hc>;
/* By symmetry, this extension always suffices */
K := SplittingField(h(G[#G]));
SK := BaseExtend(S, K);
P := Eltseq(Points(SK)[1]);

// TODO: See if this way to circumvent precision problems ever gives trouble
assert exp le 1;
PK := PuiseuxSeriesRing(K, 2);
wK := PK.1^exp;
return [ P[i] * wK + O(wK^2) : i in [1..gY] ], h(G[#G]);

end function;


function DevelopPoint(X, P0, n)
/*
 * Input:   An algebraic relation f between two variables,
 *          a point P0 that satisfies this,
 *          and the requested number of digits n.
 * Output:  A corresponding development of both components to O (n).
 *
 * The relation f has to be non-singular when developing in y.
 * The x-coordinate x0 of P can be specified as a constant element or as a
 * Puiseux series. In the latter case, the value for y is determined directly;
 * in the former, we consider x0 + t and find the corresponding value of y.
 *
 * (So this function does two things. Yes, that is bad practice...)
 */

if Type(P0[1]) in [ RngSerPuisElt, RngSerPowElt ] then
    F := BaseRing(Parent(P0[1]));
else
    F := Parent(P0[1]);
end if;
f := X`DEs[1]; df := Derivative(f, 2);
/* Note that x0 and y0 can still be Puiseux series */
x0 := P0[1]; y0 := P0[2];

/* Take constant terms and return them if n = 0 */
PR := PuiseuxSeriesRing(F, 1);
x := PR ! Coefficient(PR ! x0, 0); y := PR ! Coefficient(PR ! y0, 0);
if n eq 0 then
    return [x, y];
end if;

log := 0;
while log le Ceiling(Log(2, n)) - 1 do
    prec := Minimum(2^(log + 1), n);
    /* prec + 1 might seem more logical, but that messes up the final term,
     * which then requires a new iteration to be made */
    PR := PuiseuxSeriesRing(F, prec);
    if x0 in F then
        e := 1;
        x := x0 + PR.1;
    else
        e := ExponentDenominator(x0);
        x := LiftPuiseuxSeries(x0, PR, e);
    end if;
    y := LiftPuiseuxSeries(y, PR, e);
    y +:= -Evaluate(f, [x, y])/Evaluate(df, [x, y]);
    log +:= 1;
end while;
return [x, y];

end function;


function InitializeLift(X, Y, M)
/*
 * Input:   Curves X and Y and a normalized matrix M.
 * Output:  The very first terms in the development of P and the corresponding
 *          branches Q_j. Note that this result can contain some superfluous terms.
 */

P0 := X`P0; Q0 := Y`P0;
tjs0, f := InitializeImageBranch(M, X`echelon_exps);
PR := Parent(tjs0[1]);

/* Creating P */
P := [ PR ! P0[1], PR ! P0[2] ];
P[1] +:= PR.1;
// TODO: +1 should be good enough but OK
P := [ PR ! c : c in DevelopPoint(X, P, X`g + 3) ];

/* Creating Qs */
Qs := [ [ PR ! Q0[1], PR ! Q0[2] ] : i in [1..Y`g] ];
for i in [1..Y`g] do
    Qs[i][1] +:= tjs0[i];
end for;
/* Make sure precision buffer is always large enough to see first term */
Qs := [ [ PR ! c : c in DevelopPoint(Y, Qj, X`g + 3) ] : Qj in Qs ];

/* Fill out small terms */
IterateLift := CreateLiftIteratorFunction(X, Y, M);
while true do
    Pnew, Qsnew := IterateLift(P, Qs, X`g + 3);
    if Pnew eq P and Qsnew eq Qs then
        P := Pnew; Qs := Qsnew; break;
    end if;
    P := Pnew; Qs := Qsnew;
end while;
return P, Qs, f;

end function;


function CreateLiftIteratorFunction(X, Y, M)
/*
 * Input:   Curves X and Y and a normalized matrix M.
 * Output:  An iterator that refines the Puiseux expansion upon application.
 */
// NOTE: It is not absolutely necessary to use normalized differentials here.

fX := X`DEs[1]; dfX := Derivative(fX, X`RA.2); BX := X`NormB; gX := X`g;
fY := Y`DEs[1]; dfY := Derivative(fY, Y`RA.2); BY := Y`NormB; gY := Y`g;
e := Denominator(PuiseuxLeadingExponent(M, X`echelon_exps));

    function IterateLift(P, Qs, n);
    /*
     * Input:   Points P, branches Q, and a precision n.
     * Output:  A refinement to twice the current precision or n, whichever is
     *          smaller.
     */

    /* Create ring of higher precision: */
    K := BaseRing(Parent(P[1]));
    prec := Minimum(2*Precision(Parent(P[1])), n);
    PR := PuiseuxSeriesRing(K, prec);

    /* Lift Qs: */
    Qs := [ [ LiftPuiseuxSeries(c, PR, e) : c in Qj ] : Qj in Qs ];
    for i in [1..gY] do
        h := -Evaluate(fY, Qs[i])/Evaluate(dfY, Qs[i]);
        Qs[i][2] +:= h;
    end for;

    /* Calculate P to higher precision: */
    P := [ LiftPuiseuxSeries(c, PR, 1) : c in P ];
    P[1] := Coefficient(P[1], 0) + PR.1;
    h := -Evaluate(fX, P)/Evaluate(dfX, P);
    P[2] +:= h;

    /*
    vprint EndoCheck, 3: "Before lift:";
    vprint EndoCheck, 3: Qs;
    vprint EndoCheck, 3: P;
    vprint EndoCheck, 3: "After lift:";
    vprint EndoCheck, 3: Qs;
    vprint EndoCheck, 3: P;
    */

    /* Calculate LHS: */
    dtjs := [ Derivative(Qj[1]) : Qj in Qs ];
    BQs := Matrix([ [ Evaluate(BY[i], Qs[j]) : j in [1..gY] ] : i in [1..gY] ]);
    F_ev := Matrix([ [ Integral(&+[ BQs[i,j] * dtjs[j] : j in [1..gY] ]) : i in [1..gY] ] ]);
    DF_ev := Transpose(BQs);

    /*
    vprint EndoCheck, 3: "Determinant:";
    vprint EndoCheck, 3: Determinant(DF_ev);
    */
    if not IsInvertible(DF_ev) then
        error "Jacobian of Hensel lift is not invertible, so Puiseux lift fails";
    end if;

    /* Calculate RHS: */
    BP := [ Evaluate(BX[j], P) : j in [1..gX] ];
    G_ev := Matrix(PR, [ [ Integral(&+[ M[i,j] * BP[j] : j in [1..gX] ]) : i in [1..gY] ] ]);

    /* Calculate Hensel correction: */
    H := -(F_ev - G_ev) * DF_ev^(-1);

    /* Calculate Qs to higher precision: */
    for i in [1..gY] do
        Qs[i][1] +:= H[1,i];
        h := -Evaluate(fY, [Qs[i][1], Qs[i][2]])/Evaluate(dfY, [Qs[i][1], Qs[i][2]]);
        Qs[i][2] +:= h;
    end for;
    return P, Qs;

    end function;

return IterateLift;

end function;


intrinsic InitializedIterator(X::Crv, Y::Crv, M::., n::RngIntElt : MaxPrec := Infinity()) -> .
{Given curves X and Y, a matrix M that gives the tangent representation of a
homomorphism of Jacobians on the normalized basis of differentials, and an
integer n, returns a development of the branches to precision at least O(n)
plus or minus a negligible amount of digits. This development can be iterated
further. A minimal polynomial of the required field extension is also
returned. An optional input MaxPrec is to bound the precision in cases where we
know such a bound.}

e := Denominator(PuiseuxLeadingExponent(M, X`echelon_exps));
P, Qs, f := InitializeLift(X, Y, M);
IterateLift := CreateLiftIteratorFunction(X, Y, M);
/* Foolproof: continue until stabilization (costs one more calculation for a
 * lot more certainty) */
while true do
    Pnew, Qsnew := IterateLift(P, Qs, n);
    if Pnew eq P and Qsnew eq Qs then
        P := Pnew; Qs := Qsnew; break;
    end if;
    P := Pnew; Qs := Qsnew;
end while;
return [* P, Qs, IterateLift, MaxPrec *], f;

end intrinsic;


intrinsic IterateIterator(Iterator::List) -> .
{Applies Iterator to add maximal possible precision that does not get lost
later on.}
// TODO: Build in corresponding asserts

P, Qs, IterateLift, MaxPrec := Explode(Iterator);
e := Maximum(&cat[ [ ExponentDenominator(c) : c in Q ] : Q in Qs ]);
L := [ c : c in &cat(Qs) | not IsWeaklyZero(c) ];
prec := Minimum([ AbsolutePrecision(c) - 1/e : c in L ] cat [ MaxPrec ]);
/* May want to include bound here too, but for now that is useless */
P, Qs := IterateLift(P, Qs, Infinity());
/* This coercion seems slightly inefficient, but we take it */
PR := PuiseuxSeriesRing(BaseRing(Parent(P[1])), Integers() ! (2*((e*prec) - 1) + 1));
P := [ PR ! c : c in P ]; Qs := [ [ PR ! c : c in Q ] : Q in Qs ];
return [* P, Qs, IterateLift, MaxPrec *];

end intrinsic;
