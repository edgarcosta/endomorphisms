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
import "Conventions.m": ExtractHomomorphismsRing, VariableOrder;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "Initialize.m": InitializeCurve, ChangeTangentAction;
import "RiemannRoch.m": RRBasis, RRGenerators, RREvaluations, GlobalGenerators;
import "RiemannRoch.m": ProductEvaluations, GlobalProductBasis, GlobalProductBasisAlt;


forward InfinitesimalEquationVectors;

forward CheckVanishing;
forward CheckMultiplicityAtPoint;
forward CheckMultiplicityAwayFromPoint;
forward CheckMultiplicity;
forward GlobalScheme;
forward CheckDimension;

forward DivisorFromMatrixByDegree;
forward DivisorFromMatrixRRGlobal;
forward DivisorFromMatrixRRSplit;


function InfinitesimalEquationVectors(X, Y, P, Qs, d)
/*
 * Input:   Two curves X and Y,
 *          a basis of divisor equations fs,
 *          the precision n used when determining these,
 *          and branch expansions P and Qs.
 * Output:  The irreducible components that fit the given data.
 */

e := Maximum([ Maximum([ Denominator(Valuation(c - Coefficient(c, 0))) : c in Q ]) : Q in Qs ]);
evss := ProductEvaluations(X, Y, d, P, Qs);
min := Minimum([ Valuation(ev) : ev in &cat(evss) ]);
max := Minimum([ AbsolutePrecision(ev) : ev in &cat(evss) ]);
// TODO: Remove this safety margin, and the branch overcalculation
//max := max - 1;
M := Matrix([ &cat[ [ X`F ! Coefficient(evs[i], j/e) : j in [(e*min)..(e*max - 1)] ] : evs in evss ] : i in [1..#evss[1]] ]);
return [ Eltseq(b) : b in Basis(Kernel(M)) ];

end function;


function CheckVanishing(X, Y, d, vs, P, Qs)
// Check if points lie on scheme
// TODO: Remove this check where it is superfluous

evss := ProductEvaluations(X, Y, d, P, Qs);
for evs in evss do
    for v in vs do
        if not IsWeaklyZero(&+[ v[i]*evs[i] : i in [1..#evs] ]) then
            return false;
        end if;
    end for;
end for;
return true;

end function;


function CheckMultiplicityAtPoint(X, Y, d, vs : Margin := 2^6)
// Checks for multiplicity of vertical intersection

vprint EndoCheck, 4 : "CheckMultiplicityAtPoint...";
// TODO: Keep refining precision
precP := d + X`g + 1 + Margin; precQ := 2*Y`g + 1 + Margin;
P := DevelopPoint(X, X`P0, precP); Q := DevelopPoint(Y, Y`P0, precQ);
K<piP> := Parent(P[1]); L<piQ> := Parent(Q[1]);
/* We make a relative extension and coerce to it */
S<t> := PuiseuxSeriesRing(K);
xs := RREvaluations(X, d + X`g + 1, P); ys := RREvaluations(Y, 2*Y`g + 1, Q);
ys_hom := [ ];
for y in ys do
    coeffs, offset := Coefficients(y);
    Append(~ys_hom, &+[ coeffs[i]*t^(i + offset - 1) : i in [1..#coeffs] ]);
end for;
evs := [ x*y : x in xs, y in ys_hom ];
eqs := [ &+[ v[i]*evs[i] : i in [1..#evs] ] : v in vs ];
eqs := [ piP^(-Minimum([ Valuation(c) : c in Coefficients(q) ])) * t^(-Valuation(q)) * q : q in eqs ];
for q in eqs do
    coeffs := Coefficients(q);
    n := d*(Y`g) + Margin;
    vals := [ Valuation(coeffs[i]) : i in [1..(Y`g + 1)] ];
    vprint EndoCheck, 4 : "Valuations of coefficients:";
    vprint EndoCheck, 4 : vals;
    min, ind := Minimum(vals);
    if ind eq (Y`g + 1) then
        if &and[ Valuation(coeffs[i]) gt min : i in [1..Y`g] ] then
            return true;
        end if;
    end if;
end for;
return false;

end function;


function CheckMultiplicityAwayFromPoint(X, Y, d, vs : Margin := 2^6)
// Checks for multiplicity of vertical intersection

vprint EndoCheck, 4 : "CheckMultiplicityAwayFromPoint...";
// TODO: Keep refining precision; seem to need Degree (res) / e
prec := d + X`g + Margin; P := DevelopPoint(X, X`P0, prec); K<pi> := Parent(P[1]);
R<u,v> := PolynomialRing(K, 2); S<t> := PolynomialRing(K); h := hom< R -> S | [t, 1] >;

xs := RREvaluations(X, d + X`g + 1, P); ys := [ Y`KA ! b : b in RRBasis(Y, 2*(Y`g) + 1) ];
den := LCM([ Denominator(y) : y in ys ]); ys := [ R ! Y`RA ! (den * y) : y in ys ];
evs := [ x*y : x in xs, y in ys ]; eqs := [ &+[ v[i]*evs[i] : i in [1..#evs] ] : v in vs ];
eqs := [ q / LeadingCoefficient(q) : q in eqs ];
g := R ! Y`DEs[1]; g /:= LeadingCoefficient(g);

gcd := Resultant(g, eqs[1], v); gcd /:= LeadingCoefficient(gcd);
i := 1;
repeat
    vprint EndoCheck, 4 : "Number of elements tried:";
    vprint EndoCheck, 4 : i;
    vprint EndoCheck, 4 : "GCD:";
    vprint EndoCheck, 4 : gcd;
    coeffs := Coefficients(gcd); test := true;
    for coeff in coeffs[2..#coeffs] do
        if Valuation(coeff) le 0 then
            test := false;
        end if;
    end for;
    if test then
        return true;
    end if;
    i +:= 1;
    if i le #eqs then
        res := Resultant(g, eqs[i], v); res /:= LeadingCoefficient(res);
        gcd := GCD(gcd, res); gcd /:= LeadingCoefficient(gcd);
    end if;
until i gt #eqs;
return false;

end function;


function CheckMultiplicity(X, Y, d, vs)
// Checks for multiplicity of vertical intersection

if not CheckMultiplicityAtPoint(X, Y, d, vs) then
    return false;
end if;
return CheckMultiplicityAwayFromPoint(X, Y, d, vs);

end function;


function GlobalScheme(X, Y, d, vs)
// Find equations in affine space from vector with elements of fraction field

hX, hY := ExtractHomomorphismsRing(X, Y); RAXY := Codomain(hX);
B := [ RAXY ! b : b in GlobalProductBasis(X, Y, d) ];
eqs := [ &+[ v[i]*B[i] : i in [1..#B] ] : v in vs ];
fs := [ hX(X`DEs_sub[1]), hY(Y`DEs_sub[1]) ];
return Scheme(AffineSpace(RAXY), eqs cat fs);

/*
B := ProductBasis(X, Y, d);
den := LCM([ Denominator(b) : b in B ]);
RAXY := Integers(Parent(B[1]));
B := [ RAXY ! (den * b) : b in B ];
eqs := [ &+[ v[i]*B[i] : i in [1..#B] ] : v in vs ];
hX, hY := ExtractHomomorphismsRing(X, Y);
fs := [ RAXY ! hX(X`DEs[1]), RAXY ! hY(Y`DEs[1]) ];
return Scheme(AffineSpace(RAXY), eqs cat fs);
*/

end function;


function CheckDimension(X, Y, d, vs);
// Check if dimension equals 1

D := GlobalScheme(X, Y, d, vs);
RAXY := Parent(DefiningEquations(D)[1]);
varord := VariableOrder(); var := RAXY.varord[1];
// TODO: This step could be costly; can it be sped up via a trick?
if EliminationIdeal(DefiningIdeal(D), { var }) eq ideal< RAXY | 0 > then
    return true, D, vs;
end if;
return false, [ ], [ ];

end function;


function DivisorFromMatrixByDegree(X, Y, NormM, d : Margin := 2^6)

vprintf EndoCheck, 2 : "Trying degree %o...\n", d;
n := d*(Y`g) + Margin;
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
vs := InfinitesimalEquationVectors(X, Y, P, Qs, d);
vprintf EndoCheck, 2 : "done.\n";
if #vs eq 0 then
    return false, [ ], [ ];
end if;

vprintf EndoCheck, 2 : "Checking:\n";
vprintf EndoCheck, 2 : "Multiplicity... ";
test2 := CheckMultiplicity(X, Y, d, vs);
vprintf EndoCheck, 2 : "done.\n";
if test2 then
    vprintf EndoCheck, 2 : "Dimension... ";
    test3, D := CheckDimension(X, Y, d, vs);

    /*
    print "D:";
    eqs := DefiningEquations(D);
    for q in eqs do
        for Q in Qs do
            print Valuation(Evaluate(q, [ Q[2]/Q[1]^3, P[2]/P[1]^3, 1/Q[1], 1/P[1] ]));
        end for;
    end for;
    print "Irreducible components:";
    for C in IrreducibleComponents(D) do
        print "Trying one...";
        eqs := DefiningEquations(C);
        for q in eqs do
            for Q in Qs do
                print Valuation(Evaluate(q, [ Q[2]/Q[1]^3, P[2]/P[1]^3, 1/Q[1], 1/P[1] ]));
            end for;
        end for;
    end for;
    */

    vprintf EndoCheck, 2 : "done.\n";
    if test3 then
        vprintf EndoCheck, 2 : "Divisor found!\n";
        return true, D, vs;
    end if;
end if;
return false, [ ], [ ];

end function;


intrinsic DivisorFromMatrixRRGlobal(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^6, LowerBound := 1, UpperBound := Infinity()) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

InitializeCurve(X, P0); InitializeCurve(Y, Q0);
X`RRgens := RRGenerators(X);
X`globgens, X`DEs_sub := GlobalGenerators(X);

NormM := ChangeTangentAction(X, Y, M);
vprintf EndoCheck, 3 : "Tangent representation:\n";
vprint EndoCheck, 3 : NormM;
NormM := Y`T * NormM * (X`T)^(-1);
vprintf EndoCheck, 3 : "Normalized tangent representation:\n";
vprint EndoCheck, 3 : NormM;

d := LowerBound;
while true do
    found, D, vs := DivisorFromMatrixByDegree(X, Y, NormM, d : Margin := Margin);
    if found then
        return true, D, vs;
    end if;
    /* If that does not work, give up and try one degree higher */
    d +:= 1;
    if d gt UpperBound then
        return false, [], [];
    end if;
end while;

end intrinsic;


intrinsic DivisorFromMatrixRRSplit(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^6, LowerBound := 1, UpperBound := Infinity(), B := 300) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

/* We start at a suspected estimate and then increase degree until we find an appropriate divisor */
InitializeCurve(X, P0); InitializeCurve(Y, Q0);
X`RRgens := RRGenerators(X);
X`globgens, X`DEs_sub := GlobalGenerators(X);

NormM := ChangeTangentAction(X, Y, M);
vprintf EndoCheck, 3 : "Tangent representation:\n";
vprint EndoCheck, 3 : NormM;
NormM := Y`T * NormM * (X`T)^(-1);
vprintf EndoCheck, 3 : "Normalized tangent representation:\n";
vprint EndoCheck, 3 : NormM;
tjs0, f := InitializeImageBranch(NormM);

/* Some global elements needed below */
F := X`F; rF := X`rF; OF := X`OF; BOF := X`BOF;
RAprod := PolynomialRing(X`F, 4, "lex");
// TODO: We cannot usually take X`g. Find out what this should be.
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, X`g + 2);

ps_rts := [ ]; prs := [ ]; vss_red := [* *];
I := ideal<X`OF | 1>;

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
        found, D_red, vs_red := DivisorFromMatrixByDegree(X_red, Y_red, NormM_red, d : Margin := Margin);
        /* If that does not work, give up and try one degree higher. Note that
         * d is initialized in the outer loop, so that we keep the degree that
         * works. */
        if found then
            break;
        end if;
        d +:= 1;
        if d gt UpperBound then
            return false, [], [];
        end if;
    end while;
    Append(~vss_red, vs_red);

    vprintf EndoCheck : "Fractional CRT... ";
    vs := [ ];
    for i in [1..#vss_red[1]] do
        v_reds := [* vs_red[i] : vs_red in vss_red *];
        v := [ FractionalCRTSplit([* v_red[j] : v_red in v_reds *], prs, OF, I, BOF, BI, F) : j in [1..#v_reds[1]] ];
        Append(~vs, v);
    end for;
    vprintf EndoCheck : "done.\n";
    print vs;

    vprintf EndoCheck : "Checking:\n";
    vprintf EndoCheck : "Vanishing... ";
    /* Note that P and Qs are calculated at the beginning of this function */
    test1 := CheckVanishing(X, Y, d, vs, P, Qs);
    vprintf EndoCheck : "done.\n";

    if test1 then
        vprintf EndoCheck : "Multiplicity... ";
        test2 := CheckMultiplicity(X, Y, d, vs);
        vprintf EndoCheck : "done.\n";
        if test2 then
            vprintf EndoCheck, 2 : "Dimension... ";
            test3, D := CheckDimension(X, Y, d, vs);
            vprintf EndoCheck, 2 : "done.\n";
            if test3 then
                vprintf EndoCheck, 2 : "Divisor found!\n";
                return true, D, vs;
            end if;
        end if;
    end if;
end while;

end intrinsic;
