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


/* May want to include bound here too, but for now that is useless */
import "Branches.m": InitializeImageBranch, DevelopPoint;
import "Conventions.m": ExtractHomomorphismsRing, VariableOrder;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "Initialize.m": InitializeCurve, ChangeTangentAction;
import "RiemannRoch.m": DegreeBound;
import "RiemannRoch.m": RRBasis, RRGenerators, RREvaluations;
import "RiemannRoch.m": ProductBasis, ProductEvaluations;
import "RiemannRoch.m": GlobalGenerators, GlobalBasis, GlobalProductBasis;


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


function InfinitesimalEquationVectors(X, Y, d, P, Qs)
/*
 * Input:   Two curves X and Y,
 *          a degree d,
 *          and branch expansions P and Qs.
 * Output:  The irreducible components that fit the given data.
 */

e := Maximum([ Maximum([ Denominator(Valuation(c - Coefficient(c, 0))) : c in Q ]) : Q in Qs ]);
evss := ProductEvaluations(X, Y, d, P, Qs);
min := Minimum([ Valuation(ev) : ev in &cat(evss) ]);
max := Minimum([ AbsolutePrecision(ev) : ev in &cat(evss) ]);
// TODO: Remove this safety margin, and any branch overcalculation
max -:= 10;
M := Matrix([ &cat[ [ X`F ! Coefficient(evs[i], j/e) : j in [(e*min)..(e*max - 1)] ] : evs in evss ] : i in [1..#evss[1]] ]);

/*
print "Denominator:", e;
print "Min:", min;
print "Max:", max;
print "Number of rows:", #Rows(M);
print "Number of columns:", #Rows(Transpose(M));
print "Dim Ker:", Dimension(Kernel(M));
*/

return [ Eltseq(b) : b in Basis(Kernel(M)) ];

end function;


function CheckVanishing(X, Y, d, vs, P, Qs)
// Check if points lie on scheme

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


// TODO: Currently not used, check degree bound equivalents
function CheckMultiplicityAtPoint(X, Y, d, vs : Margin := 2^8)
// Checks for multiplicity of vertical intersection

vprint EndoCheck, 3 : "";
vprint EndoCheck, 3 : "CheckMultiplicityAtPoint...";
precP := DegreeBound(X, d) + Margin; precQ := DegreeBound(Y, Y`g) + Margin;
P := DevelopPoint(X, X`P0, precP); Q := DevelopPoint(Y, Y`P0, precQ);
K<piP> := Parent(P[1]); L<piQ> := Parent(Q[1]);
/* We make a relative extension and coerce to it */
S<t> := PuiseuxSeriesRing(K);
xs := RREvaluations(X, DegreeBound(X, d), P); ys := RREvaluations(Y, DegreeBound(Y, Y`g), Q);
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
    vprint EndoCheck, 3 : "Valuations of coefficients:";
    vprint EndoCheck, 3 : vals;
    min, ind := Minimum(vals);
    if ind eq (Y`g + 1) then
        if &and[ Valuation(coeffs[i]) gt min : i in [1..Y`g] ] then
            return true;
        end if;
    end if;
end for;
return false;

end function;


function MakeDenseAndMonic(p)
// Removes leading non-zero coefficients in polynomials over Puiseux series
// rings and makes these monic

R := Parent(p); d := Degree(p);
dmp := R ! 0;
for i in [0..d] do
    c := Coefficient(p, i);
    if not IsWeaklyZero(c) then
        dmp +:= c*R.1^i;
    end if;
end for;
if Valuation(dmp) eq Infinity() then
    return dmp;
end if;
return dmp / LeadingCoefficient(dmp);

end function;


function GCDP(a, b)
// GCD of polynomials over Puiseux series rings

a := MakeDenseAndMonic(a); b := MakeDenseAndMonic(b);
va := Valuation(a); vb := Valuation(b);
if va eq Infinity() then
    return b;
end if;
if vb eq Infinity() then
    return a;
end if;
da := Degree(a); db := Degree(b);
if da le db then
    anew := a; bnew := b;
    danew := da; dbnew := db;
else
    anew := b; bnew := a;
    danew := db; dbnew := da;
end if;
R := Parent(anew);
while dbnew ge danew do
    bnew -:= Coefficient(bnew, dbnew)*R.1^(dbnew - danew)*anew;
    dbnew -:= 1;
end while;
return GCDP(anew, bnew);

end function;


function CheckMultiplicityAwayFromPoint(X, Y, d, vs : Margin := 2^8)
// Checks for multiplicity of vertical intersection
// TODO: Think more about precision

vprint EndoCheck, 3 : "Checking multiplicity away from point...";
prec := DegreeBound(X, d) + Margin; P := DevelopPoint(X, X`P0, prec); K<pi> := Parent(P[1]);
R<u,v> := PolynomialRing(K, 2); S<t> := PolynomialRing(K); h := hom< R -> S | [t, 1] >;

/* The next line must use the same bounds as elsewhere */
xs := RREvaluations(X, DegreeBound(X, d), P); ys := [ Y`KA ! b : b in RRBasis(Y, DegreeBound(Y, Y`g)) ];
den := LCM([ Denominator(y) : y in ys ]); ys := [ R ! Y`RA ! (den * y) : y in ys ];
evs := [ x*y : x in xs, y in ys ]; eqs := [ &+[ v[i]*evs[i] : i in [1..#evs] ] : v in vs ];
fY := R ! Y`DEs[1];

gcd := MakeDenseAndMonic(h(Resultant(fY, eqs[1], v)));
i := 1;
repeat
    vprint EndoCheck, 3 : "Number of elements tried:";
    vprint EndoCheck, 3 : i;
    vprint EndoCheck, 3 : "GCD:";
    vprint EndoCheck, 3 : gcd;
    d := Degree(gcd);
    test_away := true;
    /* Check g highest coefficients apart from the leading one */
    for n in [(d - Y`g)..(d - 1)] do
        if Valuation(Coefficient(gcd, n)) le 0 then
            test_away := false;
            break;
        end if;
    end for;
    if test_away then
        vprint EndoCheck, 3 : "Checking multiplicity at point...";
        test_at := true;
        for n in [0..(d - Y`g - 1)] do
            if not IsWeaklyZero(Coefficient(gcd, n)) then
                test_at := false;
                break;
            end if;
        end for;
        if test_at then
            return true;
        end if;
    end if;
    i +:= 1;
    if i le #eqs then
        res := MakeDenseAndMonic(h(Resultant(fY, eqs[i], v)));
        gcd := MakeDenseAndMonic(GCDP(gcd, res));
    end if;
until i gt #eqs;
return false;

end function;


function CheckMultiplicity(X, Y, d, vs)
// Checks for multiplicity of vertical intersection

return CheckMultiplicityAwayFromPoint(X, Y, d, vs);
//if not CheckMultiplicityAtPoint(X, Y, d, vs) then
//    return false;
//end if;
//return CheckMultiplicityAwayFromPoint(X, Y, d, vs);

end function;


function GlobalScheme(X, Y, d, vs)
// Find equations in affine space from vector with elements of fraction field

/* This version will have a lot of parasitic fibral components */
B := ProductBasis(X, Y, d);
den := LCM([ Denominator(b) : b in B ]);
RAXY := Integers(Parent(B[1]));
B := [ RAXY ! (den * b) : b in B ];
eqs := [ &+[ v[i]*B[i] : i in [1..#B] ] : v in vs ];
hX, hY := ExtractHomomorphismsRing(X, Y);
fs := [ RAXY ! hX(X`DEs[1]), RAXY ! hY(Y`DEs[1]) ];
return Scheme(AffineSpace(RAXY), eqs cat fs);

hX, hY := ExtractHomomorphismsRing(X, Y); RAXY := Codomain(hX);
B := [ RAXY ! b : b in GlobalProductBasis(X, Y, d) ];
eqs := [ &+[ v[i]*B[i] : i in [1..#B] ] : v in vs ];
fs := [ hX(X`DEs_sub[1]), hY(Y`DEs_sub[1]) ];
return Scheme(AffineSpace(RAXY), eqs cat fs);

end function;


function CheckDimension(X, Y, d, vs);
// Check if dimension equals 1

D := GlobalScheme(X, Y, d, vs);
// TODO: We ignore this for now because this check appears to be difficult,
// which is partly due to deficiencies in Magma.
return true, D;
RAXY := Parent(DefiningEquations(D)[1]);
varord := VariableOrder(); var := RAXY.varord[1];
// TODO: This step could be costly; can it be sped up via a trick?
if EliminationIdeal(DefiningIdeal(D), { var }) eq ideal< RAXY | 0 > then
    return true, D;
end if;
return false, [ ];

end function;


function DivisorFromMatrixByDegree(X, Y, NormM, d : Margin := 2^8)

vprint EndoCheck, 2 : "Trying degree", d;
/* Cardinality of basis of functions plus margin: */
n := (DegreeBound(X, d) + 1 - X`g)*(DegreeBound(Y, Y`g) + 1 - Y`g) + Margin;
vprint EndoCheck, 2 : "Number of terms in expansion:", n;

/* Take non-zero image branch */
vprint EndoCheck, 2 : "Expanding...";
P, Qs := InitializedIterator(X, Y, NormM, n);
_<t> := Parent(P[1]);
_<r> := BaseRing(Parent(P[1]));
vprint EndoCheck, 2 : "done.";

/* Fit a divisor to it */
vprint EndoCheck, 2 : "Solving linear system...";
vs := InfinitesimalEquationVectors(X, Y, d, P, Qs);
vprint EndoCheck, 2 : "done.";
if #vs eq 0 then
    return false, [ ], [ ];
end if;

vprint EndoCheck, 2 : "Checking:";
vprint EndoCheck, 2 : "Multiplicity...";
test_mult := CheckMultiplicity(X, Y, d, vs);
vprint EndoCheck, 2 : "done.";
if test_mult then
    vprint EndoCheck, 2 : "Dimension...";
    test_dim, D := CheckDimension(X, Y, d, vs);

    vprint EndoCheck, 2 : "done.";
    if test_dim then
        vprint EndoCheck, 2 : "";
        vprint EndoCheck, 2 : "Divisor found!";
        return true, D, vs;
    end if;
end if;
return false, [ ], [ ];

end function;


intrinsic DivisorFromMatrixRRGlobal(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^8, LowerBound := 1, UpperBound := Infinity()) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

InitializeCurve(X, P0); InitializeCurve(Y, Q0);
X`RRgens := RRGenerators(X);
//X`globgens, X`DEs_sub := GlobalGenerators(X);

NormM := ChangeTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);

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


intrinsic DivisorFromMatrixRRSplit(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^8, LowerBound := 1, UpperBound := Infinity(), B := 300) -> BoolElt, .
{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}

/* We start at a suspected estimate and then increase degree until we find an appropriate divisor */
InitializeCurve(X, P0); InitializeCurve(Y, Q0);
X`RRgens := RRGenerators(X);
//X`globgens, X`DEs_sub := GlobalGenerators(X);

NormM := ChangeTangentAction(X, Y, M);
NormM := Y`T * NormM * (X`T)^(-1);
tjs0, f := InitializeImageBranch(NormM);

/* Some global elements needed below */
F := X`F; OF := X`OF;
//RAprod := PolynomialRing(X`F, 4, "lex");
RAprod := PolynomialRing(X`F, 4);
// TODO: We cannot usually take X`g. Find out what this should be instead.
P, Qs := InitializedIterator(X, Y, NormM, X`g + 2);

prs := [ ]; vss_red := [* *];
I := ideal<X`OF | 1>;

d := LowerBound;
while true do
    /* Find new prime */
    repeat
        pr, h := RandomSplitPrime(f, B);
    until not pr in prs;
    Append(~prs, pr); I *:= pr;
    vprint EndoCheck : "Split prime over", #Codomain(h);

    /* Add corresponding data */
    X_red := ReduceCurveSplit(X, h); Y_red := ReduceCurveSplit(Y, h);
    NormM_red := ReduceMatrixSplit(NormM, h);

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

    vprint EndoCheck : "Fractional CRT...";
    vs := [ ];
    for i in [1..#vss_red[1]] do
        v_reds := [* vs_red[i] : vs_red in vss_red *];
        v := [ FractionalCRTSplit([* v_red[j] : v_red in v_reds *], prs) : j in [1..#v_reds[1]] ];
        Append(~vs, v);
    end for;
    vprint EndoCheck : "done.";

    vprint EndoCheck : "Checking:";
    vprint EndoCheck : "Vanishing...";
    /* Note that P and Qs are calculated at the beginning of this function */
    test_van := CheckVanishing(X, Y, d, vs, P, Qs);
    vprint EndoCheck : "done.";

    if test_van then
        vprint EndoCheck : "Multiplicity...";
        test_mult := CheckMultiplicity(X, Y, d, vs);
        vprint EndoCheck : "done.";
        if test_mult then
            vprint EndoCheck, 2 : "Dimension...";
            test_dim, D := CheckDimension(X, Y, d, vs);
            vprint EndoCheck, 2 : "done.";
            if test_dim then
                vprint EndoCheck, 2 : "Divisor found!";
                return true, D, vs;
            end if;
        end if;
    end if;
end while;

end intrinsic;
