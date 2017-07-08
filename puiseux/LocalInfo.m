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


forward LiftPuiseuxSeries;

forward InitializeMatrix;
forward PuiseuxRamificationIndex;
forward InitializeImageBranch;

forward DevelopPoint;

forward InitializeLift;
forward CreateLiftIterator;
forward ApproximationsFromTangentAction;


function LiftPuiseuxSeries(f, PR, e)

L := [ Coefficient(f, i/e)*PR.1^(i/e) : i in [0..e*AbsolutePrecision(f) - 1] ];
if #L eq 0 then
    return PR ! 0;
end if;
return PR ! (&+L);

end function;


function PuiseuxRamificationIndex(M)
/*
 * Input:   A matrix M that represents an endomorphism,
 *          after normalizing to an upper triangular form.
 * Output:  The ramification index of the corresponding Puiseux expansions.
 */

g := #Rows(M);
e := g;
while true do
    n := g div e;
    test := &or[ M[e*i, i] ne 0 : i in [1..n] ];
    if test then
        return e;
    end if;
    e := e - 1;
    /* Stopping criterion in case of ramification: */
    if e eq 0 then
        return 1;
    end if;
end while;

end function;


function InitializeImageBranch(M)
/*
 * Input:   A matrix M that represents an endomorphism
 *          after normalizing to an upper triangular form.
 * Output:  The leading coefficients of the corresponding Puiseux expansions
 *          and the polynomial that gives their field of definition.
 */

/* Recovering old invariants: */
F := Parent(M[1,1]);
g := #Rows(M);
e := PuiseuxRamificationIndex(M);

if g eq 1 then
    r := Eltseq(Rows(M)[1]);
    RF := PolynomialRing(F); xF := RF.1; PF := PowerSeriesRing(F, #r + 1); tF := PF.1;
    return [ &+[ (r[n] / n) * tF^n : n in [1..#r] ] + O(tF^(#r + 1)) ], xF - 1;
end if;

/* Normalized equations (depend only on the matrix): */
A := AffineSpace(F, g); RA := CoordinateRing(A);
eqs := [ ];
for n in [1..g] do
    powersum := &+[ RA.i^n : i in [1..g] ];
    if n mod e eq 0 then
        Append(~eqs, powersum - e * M[n, n div e]);
    else
        Append(~eqs, powersum);
    end if;
end for;
S := Scheme(A, eqs);

/* The upcoming steps are taken to avoid the use of an algebraic closure */
RF := PolynomialRing(F);
hc := [ RF!0 : i in [1..g] ]; hc[#hc] := RF.1; h := hom<RA -> RF | hc>;

/* By symmetry, this extension always suffices */
G := GroebnerBasis(ideal<RA | eqs>);
if not IsFinite(F) then
    //K := GaloisSplittingField(h(G[#G]));
    K := SplittingField(h(G[#G]));
else
    K := SplittingField(h(G[#G]));
end if;

/* Extending and evaluating: */
SK := BaseExtend(S, K);
P := Eltseq(Points(SK)[1]);
/* Power series ring is used if possible for efficiency: */
if e eq 1 then
    PK := PowerSeriesRing(K, 2);
    wK := PK.1;
else
    PK := PuiseuxSeriesRing(K, 2);
    wK := PK.1^(1/e);
end if;
return [ P[i] * wK + O(wK^2) : i in [1..g] ], h(G[#G]);

end function;


function DevelopPoint(X, P0, n)
/*
 * Input:   An algebraic relation f between two variables,
 *          a point P0 that satisfies this,
 *          and the requested number of digits n.
 * Output:  A corresponding development of both components.
 *
 * The relation f has to be non-singular when developing in y,
 * and the x-coordinate can be specified as a Puiseux series. */

f := X`DEs[1];
if Type(P0[1]) in [ RngSerPuisElt, RngSerPowElt ] then
    F := BaseRing(Parent(P0[1]));
else
    F := Parent(P0[1]);
end if;
df := Derivative(f, 2);
x0 := P0[1]; y0 := P0[2];

PR := PuiseuxSeriesRing(F, 1);
x := PR ! Coefficient(PR ! x0, 0); y := PR ! Coefficient(PR ! y0, 0);
if n eq 0 then
    return [x, y];
end if;

log := 0;
while log le Ceiling(Log(2, n)) - 1 do
    prec := Minimum(2^(log + 1), n);
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

P0 := X`P0; Q0 := Y`P0;
e := PuiseuxRamificationIndex(M);
tjs0 := InitializeImageBranch(M);
PR := Parent(tjs0[1]);

/* Creating P */
P := [ PR ! P0[1], PR ! P0[2] ];
P[1] +:= PR.1;
P := [ PR ! c : c in DevelopPoint(X, P, X`g + 1) ];

/* Creating Qs */
Qs := [ [ PR ! Q0[1], PR ! Q0[2] ] : i in [1..Y`g] ];
for i in [1..Y`g] do
    Qs[i][1] +:= tjs0[i];
end for;
Qs := [ [ PR ! c : c in DevelopPoint(Y, Qj, X`g + 1) ] : Qj in Qs ];
return P, Qs;

end function;


function CreateLiftIterator(X, Y, M)

fX := X`DEs[1]; dfX := Derivative(fX, X`y); BX := X`NormB; gX := X`g;
fY := Y`DEs[1]; dfY := Derivative(fY, Y`y); BY := Y`NormB; gY := Y`g;

e := PuiseuxRamificationIndex(M);

    function Iterate(P, Qs, n);

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

    /* Calculate LHS: */
    dtjs := [ Derivative(Qj[1]) : Qj in Qs ];
    BQs := Matrix([ [ Evaluate(BY[i], Qs[j]) : j in [1..gY] ] : i in [1..gY] ]);
    F_ev := Matrix([ [ Integral(&+[ BQs[i,j] * dtjs[j] : j in [1..gY] ]) : i in [1..gY] ] ]);
    DF_ev := Transpose(BQs);

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

return Iterate;

end function;


intrinsic ApproximationsFromTangentAction(X::Crv, Y::Crv, M::., n::RngIntElt) -> Tup, Tup
{Given curves X and Y, a matrix M that gives the tangent representation of a
homomorphism of Jacobians on the normalized basis of differentials, and an
integer n, returns a development of the branches to precision at least O(n).}

e := PuiseuxRamificationIndex(M);
P, Qs := InitializeLift(X, Y, M);
IterateLift := CreateLiftIterator(X, Y, M);
/* TODO: Determine exact bound needed here */
for i:=1 to Ceiling(Log(2, n + e + 1)) do
    P, Qs := IterateLift(P, Qs, n + e + 1);
end for;
return P, Qs;

end intrinsic;
