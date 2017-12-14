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


import "Branches.m": InitializeImageBranch;
import "FractionalCRT.m": RandomSplitPrime, FractionalCRTSplit, ReduceMatrixSplit, ReduceCurveSplit;
import "Initialize.m": InitializeCurve, ChangeTangentAction, ExtractHomomorphismsRing, ExtractPoints, VariableOrder;
import "RiemannRoch.m": RRGenerators, RRBasis, RREvaluate;


forward MonomialEvaluations;
forward MonomialDivisors;
forward InfinitesimalEquationVectors;

forward CheckVanishing;
forward CheckMultiplicity;
forward GlobalGenerators;
forward GlobalScheme;
forward CheckDimension;

forward DivisorFromMatrixByDegree;
forward DivisorFromMatrix;
forward DivisorFromMatrixSplit;


function MonomialEvaluations(X, Y, P, Qs, d)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Evaluations for divisors of degree d coming from the ambient of X.
 */

xs := RREvaluate(X, P, d + X`g);
ys := &cat[ RREvaluate(Y, Q, 2*(Y`g)) : Q in Qs ];
return [ x*y : x in xs, y in ys ];

end function;


function MonomialDivisors(X, Y, d)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

P := X`P0; Q := Y`Q0;
xs := [ X`KA ! b : b in RRBasis(X, P, d + X`g) ];
ys := [ Y`KA ! b : b in RRBasis(Y, Q, 2*(Y`g)) ];
return [ x*y : x in xs, y in ys ];

end function;


function InfinitesimalEquationVectors(X, Y, P, Qs, d)
/*
 * Input:   Two curves X and Y,
 *          a basis of divisor equations fs,
 *          the precision n used when determining these,
 *          and branch expansions P and Qs.
 * Output:  The irreducible components that fit the given data.
 */

/* Recovering a linear system */
evs := MonomialEvaluations(X, Y, P, Qs, d);
e := Maximum([ Maximum([ Denominator(Valuation(c - Coefficient(c, 0))) : c in Q ]) : Q in Qs ]);
min := Minimum([ Valuation(ev) : ev in evs ]);
max := Minimum([ AbsolutePrecision(ev) : ev in evs ]);
M := Matrix([ [ X`F ! Coefficient(ev, i/e) : i in [(e*min)..(e*max - 1)] ] : ev in evs ]);
B := Basis(Kernel(M));
return [ Eltseq(b) : b in B ];

end function;


function CheckVanishing(X, Y, d, vs, P, Qs)
// Check if points lie on scheme

evs := MonomialEvaluations(X, Y, P, Qs, d);
for v in vs do
    if &+[ v[i]*evs[i] : i in [1..#evs] ] ne 0 then
        return false;
    end if;
end for;
return true;

end function;


function CheckMultiplicity(X, Y, d, vs)
// Checks for multiplicity of vertical intersection

return true;

end function;


function GlobalGenerators(X, Y)
// Find equations in affine space from vector with elements of fraction field

// Make product
if X`is_hyperelliptic then
    P := X`P0; g := X`g; f := X`DEs[1];
    RA<x,y> := X`RA; SA<u,v> := PolynomialRing(X`F, 2);
    u_sub := (1 / (x - P[1])); v_sub := y / (x - P[1])^(g + 1);
    x_sub := (1 / u) + P[1];   y_sub := v / u^(g + 1);
    f_sub := SA ! (u^(2*g + 2) * Evaluate(f, [x_sub, y_sub]));
    X_sub := Curve(AffineSpace(SA), f_sub);
    RX_sub := CoordinateRing(X_sub);
    KX_sub := FieldOfFractions(RX_sub);
    RA_sub := CoordinateRing(Ambient(X_sub));
    KA_sub := FieldOfFractions(RA_sub);
    gens := X`RRgens;
    //print [ Evaluate(X`KA ! gen, [x_sub, y_sub]) : gen in gens ];
    gens_sub := [ KX_sub ! Evaluate(X`KA ! gen, [x_sub, y_sub]) : gen in gens ];
    print X;
    gens_sub := [ KA_sub ! gen : gen in gens_sub ];
    print gens_sub;
    for gen in gens_sub do
        print Support(Divisor(X_sub, gen));
        print RX_sub ! KX_sub ! gen;
    end for;
    globX := [ f_sub ] cat gens_sub;
end if;
return globX;

end function;


function GlobalScheme(X, Y, d, vs)
// Find equations in affine space from vector with elements of fraction field

// Make product
print [ X`KA ! f : f in RRGenerators(X) ];
if X`is_hyperelliptic then
    A := X`A;
    f := X`DEs[1];
    eqs := [ f ];
end if;
return eqs;

end function;


function CheckDimension(X, Y, d, vs);
// Check if dimension equals 1

return true, 0;

end function;


//function DivisorFromMatrixByDegree(X, Y, NormM, d : Margin := 2^4)
//
//vprintf EndoCheck, 2 : "Trying degree %o...\n", d;
//// TODO: This number n can be tweaked
//n := d*(Y`g) + Margin;
//vprintf EndoCheck, 2 : "Number of terms in expansion: %o.\n", n;
//
///* Take non-zero image branch */
//vprintf EndoCheck, 2 : "Expanding... ";
//P, Qs := ApproximationsFromTangentAction(X, Y, NormM, n);
//vprintf EndoCheck, 4 : "Base point:\n";
//_<t> := Parent(P[1]);
//_<r> := BaseRing(Parent(P[1]));
//vprint EndoCheck, 4 : P;
//vprintf EndoCheck, 4 : "Resulting branches:\n";
//vprint EndoCheck, 4 : Qs;
//vprint EndoCheck, 4 : BaseRing(Parent(P[1]));
//vprintf EndoCheck, 2 : "done.\n";
//
///* Fit a divisor to it */
//vprintf EndoCheck, 2 : "Solving linear system... ";
//vs := InfinitesimalEquationVectors(X, Y, P, Qs, d);
//vprintf EndoCheck, 2 : "done.\n";
//
//vprintf EndoCheck, 2 : "Checking:\n";
//vprintf EndoCheck, 2 : "Step 1... ";
//test1 := CheckVanishing(X, Y, d, vs, P, Qs);
//vprintf EndoCheck, 2 : "done.\n";
//if test1 then
//    vprintf EndoCheck, 2 : "Step 2... ";
//    test2 := CheckMultiplicity(X, Y, d, vs);
//    vprintf EndoCheck, 2 : "done.\n";
//    if test2 then
//        vprintf EndoCheck, 2 : "Step 3... ";
//        test3, D := CheckDimension(X, Y, d, vs);
//        vprintf EndoCheck, 2 : "done.\n";
//        if test3 then
//            vprintf EndoCheck, 2 : "Divisor found!\n";
//            return true, D, vs;
//        end if;
//    end if;
//end if;
//return false, [ ];
//
//end function;
//
//
//intrinsic DivisorFromMatrix(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity()) -> BoolElt, .
//{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}
//
//InitializeCurve(X, P0); InitializeCurve(Y, Q0);
//NormM := ChangeTangentAction(X, Y, M);
//vprintf EndoCheck, 3 : "Tangent representation:\n";
//vprint EndoCheck, 3 : NormM;
//NormM := Y`T * NormM * (X`T)^(-1);
//vprintf EndoCheck, 3 : "Normalized tangent representation:\n";
//vprint EndoCheck, 3 : NormM;
//
//d := LowerBound;
//while true do
//    found, D, vs := DivisorFromMatrixByDegree(X, Y, NormM, d : Margin := Margin);
//    if found then
//        return true, D, vs;
//    end if;
//    /* If that does not work, give up and try one degree higher */
//    d +:= 1;
//    if d gt UpperBound then
//        return false, [], [];
//    end if;
//end while;
//
//end intrinsic;
//
//
//intrinsic DivisorFromMatrixSplit(X::Crv, P0::Pt, Y::Crv, Q0::Pt, M::. : Margin := 2^4, LowerBound := 1, UpperBound := Infinity(), B := 300) -> BoolElt, .
//{Given two pointed curves (X, P0) and (Y, Q0) along with a tangent representation of a projection morphism on the standard basis of differentials, returns a corresponding divisor (if it exists). The parameter Margin specifies how many potentially superfluous terms are used in the development of the branch, the parameter LowerBound specifies at which degree one starts to look for a divisor, and the parameter UpperBound specifies where to stop.}
//
///* We start at a suspected estimate and then increase degree until we find an appropriate divisor */
//InitializeCurve(X, P0); InitializeCurve(Y, Q0);
//NormM := ChangeTangentAction(X, Y, M);
//vprintf EndoCheck, 3 : "Tangent representation:\n";
//vprint EndoCheck, 3 : NormM;
//NormM := Y`T * NormM * (X`T)^(-1);
//vprintf EndoCheck, 3 : "Normalized tangent representation:\n";
//vprint EndoCheck, 3 : NormM;
//tjs0, f := InitializeImageBranch(NormM);
//
///* Some global elements needed below */
//F := X`F; rF := X`rF; OF := X`OF; BOF := X`BOF;
//RAprod := PolynomialRing(X`F, 4, "lex");
//P, Qs := ApproximationsFromTangentAction(X, Y, NormM, X`g);
//
//ps_rts := [ ]; prs := [ ]; vss_red := [* *];
//I := ideal<X`OF | 1>;
//
//d := LowerBound;
//while true do
//    /* Find new prime */
//    repeat
//        p_rt := RandomSplitPrime(f, B);
//        p, rt := Explode(p_rt);
//    until not p in [ tup[1] : tup in ps_rts ];
//    Append(~ps_rts, p_rt);
//    vprintf EndoCheck : "Split prime over %o\n", p;
//
//    /* Add corresponding data */
//    pr := ideal<X`OF | [ p, rF - rt ]>;
//    Append(~prs, pr); I *:= pr;
//    X_red := ReduceCurveSplit(X, p, rt); Y_red := ReduceCurveSplit(Y, p, rt);
//    NormM_red := ReduceMatrixSplit(NormM, p, rt);
//    BI := Basis(I);
//
//    while true do
//        found, D_red, vs_red := DivisorFromMatrixByDegree(X_red, Y_red, NormM_red, d : Margin := Margin);
//        /* If that does not work, give up and try one degree higher. Note that
//         * d is initialized in the outer loop, so that we keep the degree that
//         * works. */
//        if found then
//            break;
//        end if;
//        d +:= 1;
//        if d gt UpperBound then
//            return false, [];
//        end if;
//    end while;
//    Append(~vss_red, vs_red);
//
//    vprintf EndoCheck : "Fractional CRT... ";
//    vs := [ ];
//    for i in [1..#vss_red[1]] do
//        v_reds := [* vs_red[i] : vs_red in vss_red *];
//        v := [ FractionalCRTSplit([* v_red[j] : v_red in v_reds *], prs, OF, I, BOF, BI, F) : j in [1..#v_reds[1]] ];
//        Append(~vs, v);
//    end for;
//    vprintf EndoCheck : "done.\n";
//
//    vprintf EndoCheck : "Checking:\n";
//    vprintf EndoCheck : "Step 1... ";
//    /* Note that P and Qs are calculated at the beginning of this function */
//    test1 := CheckVanishing(X, Y, d, vs, P, Qs);
//    vprintf EndoCheck : "done.\n";
//
//    if test1 then
//        vprintf EndoCheck : "Step 2... ";
//        test2 := CheckMultiplicity(X, Y, d, vs);
//        vprintf EndoCheck : "done.\n";
//        if test2 then
//            vprintf EndoCheck, 2 : "Step 3... ";
//            test3, D := CheckDimension(X, Y, d, vs);
//            vprintf EndoCheck, 2 : "done.\n";
//            if test3 then
//                vprintf EndoCheck, 2 : "Divisor found!\n";
//                return true, D, vs;
//            end if;
//        end if;
//    end if;
//end while;
//
//end intrinsic;
