/***
 *  Riemann-Roch functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


import "Conventions.m": ExtractHomomorphismsField, ExtractHomomorphismsRing;


forward DegreeBound;

forward RRGenerators;
forward RRBasis;
forward RREvaluations;

forward ProductBasis;
forward ProductEvaluations;

forward GlobalGenerators;
forward GlobalBasis;
forward GlobalProductBasis;


function DegreeBound(X, d)

/* Heuristic version; this should suffices since we only care about a divisor
 * containing the base point and extra fibral contributions are all right */
return d + X`g + 1;
/* Base-point free version: */
return d + 2*X`g;
/* Very ample version: */
return d + 2*X`g + 1;

end function;


function RRGenerators(X)
// Gives all Riemann-Roch functions needed to determine the others as monomials
// in them

if assigned X`RRgens then
    return X`RRgens;
end if;
g := X`g; KU := X`KU; P0 := X`P0;
RRgens := [ KU ! 1 ];
// TODO: This does not work because not echelonized:
//V, phiV := RiemannRochSpace((2*g + 1)*Divisor(P0));
for d in [(g + 1)..(2*g + 1)] do
    V, phiV := RiemannRochSpace(d*Divisor(P0));
    for v in Basis(V) do
        f := phiV(v);
        zer, pol := SignDecomposition(Divisor(f));
        if Degree(pol) eq d then
            Append(~RRgens, KU ! f);
            break;
        end if;
    end for;
end for;
X`RRgens := RRgens;
return X`RRgens;

end function;


function RRBasis(X, d)
// Gives a basis for pole order up to d at a non-Weierstrass point

g := X`g; KU := X`KU; gens := X`RRgens;
if d le g then
    return [ gens[1] ];
end if;
if d le (2*g + 1) then
    return gens[1..(d + 1 - g)];
end if;
B := gens;
while #B lt (d + 1 - g) do
    Append(~B, B[2]*B[#B - g]);
end while;
return B;

end function;


function RREvaluations(X, d, P)
// Gives a basis for pole order up to d at a non-Weierstrass point, evaluated

g := X`g; gens := X`RRgens; KA := X`KA; F := Parent(P[1]);
if d le g then
    return [ Evaluate(KA ! gens[1], P) ];
end if;
if d le (2*g + 1) then
    return [ Evaluate(KA ! gen, P) : gen in gens[1..(d + 1 - g)] ];
end if;
evs := [ Evaluate(KA ! gen, P) : gen in gens ];
while #evs lt (d + 1 - g) do
    Append(~evs, evs[2]*evs[#evs - g]);
end while;
return evs;

end function;


function ProductBasis(X, Y, d)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

hX, hY := ExtractHomomorphismsField(X, Y);
xs := [ X`KA ! b : b in RRBasis(X, DegreeBound(X, d)) ];
ys := [ Y`KA ! b : b in RRBasis(Y, DegreeBound(Y, Y`g)) ];
return [ hX(x)*hY(y) : x in xs, y in ys ];

end function;


function ProductEvaluations(X, Y, d, P, Qs)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Evaluations for divisors of degree d coming from the ambient of X.
 */

xs := RREvaluations(X, DegreeBound(X, d), P);
yss := [ RREvaluations(Y, DegreeBound(Y, Y`g), Q) : Q in Qs ];
return [ [ x*y : x in xs, y in ys ] : ys in yss ];

end function;


function GlobalGenerators(X)
// Find equations in affine space from vector with elements of fraction field

if assigned X`globgens and assigned X`DEs_sub then
    return X`globgens, X`DEs_sub;
end if;
P := X`P0; g := X`g;
RA<x,y> := X`RA; SA<u,v> := PolynomialRing(X`F, 2);
if X`is_hyperelliptic then
    u_sub := (1 / (x - P[1])); v_sub := y / (x - P[1])^(g + 1);
    x_sub := (1 / u) + P[1];   y_sub := v / u^(g + 1);
    DEs_sub := [ SA ! (u^(2*g + 2) * Evaluate(DE, [x_sub, y_sub])) : DE in X`DEs ];
else
    u_sub := (1 / (x - P[1])); v_sub := y / (x - P[1]);
    x_sub := (1 / u) + P[1];   y_sub := v / u;
    DEs_sub := [ SA ! (u^Degree(DE) * Evaluate(DE, [x_sub, y_sub])) : DE in X`DEs ];
end if;
U_sub := Curve(AffineSpace(SA), DEs_sub);
RU_sub := CoordinateRing(U_sub); KU_sub := FunctionField(U_sub);
RA_sub := CoordinateRing(Ambient(U_sub)); KA_sub := FieldOfFractions(RA_sub);
gens_sub := [ KU_sub ! Evaluate(X`KA ! gen, [x_sub, y_sub]) : gen in X`RRgens ];
// TODO: This goes wrong in general because Magma took the easy way out.
// Gosh darn it.
gens_sub := [ RA_sub ! KA_sub ! gen : gen in gens_sub ];
den := LCM([ Denominator(gen) : gen in gens_sub ]);
gens_sub := [ RA_sub ! (den * gen) : gen in gens_sub ];
X`globgens := gens_sub; X`DEs_sub := DEs_sub;
return X`globgens, X`DEs_sub;

end function;


function GlobalBasis(X, d)

g := X`g; gens := X`globgens; RA := X`RA;
if d le g then
    return [ RA ! gens[1] ];
end if;
if d le (2*g + 1) then
    return [ RA ! gen : gen in gens[1..(d + 1 - g)] ];
end if;
B := [ RA ! gen : gen in gens ];
while #B lt (d + 1 - g) do
    Append(~B, B[2]*B[#B - g]);
end while;
return B;

end function;


function GlobalProductBasis(X, Y, d)

hX, hY := ExtractHomomorphismsRing(X, Y);
xs := [ X`RA ! b : b in GlobalBasis(X, DegreeBound(X, d)) ];
ys := [ Y`RA ! b : b in GlobalBasis(Y, DegreeBound(Y, Y`g)) ];
return [ hX(x)*hY(y) : x in xs, y in ys ];

end function;
