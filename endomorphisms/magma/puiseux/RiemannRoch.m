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


forward RRGenerators;
forward RRBasis;
forward RREvaluations;

forward ProductBasis;
forward ProductEvaluations;

forward GlobalGenerators;
forward GlobalBasis;
forward GlobalProductBasis;


function RRGenerators(X)
// Gives all Riemann-Roch functions needed to determine the others as monomials
// in them

if assigned X`RRgens then
    return X`RRgens;
end if;
g := X`g; KU := X`KU; P0 := X`P0;
V, phiV := RiemannRochSpace((2*g + 1)*Divisor(P0));
X`RRgens := [ KU ! phiV(V.i) : i in [1..(Dimension(V) - 1)] ];
return X`RRgens;

end function;


function RRBasis(X, d)
// Gives a basis for pole order up to d at a non-Weierstrass point

g := X`g; KU := X`KU; gens := X`RRgens;
if d le g then
    return [ KU ! 1 ];
end if;
if d le (2*g + 1) then
    return [ KU ! 1 ] cat gens[1..(d - g)];
end if;
B := [ KU ! 1 ] cat gens;
while #B lt (1 - g + d) do
    Append(~B, B[2]*B[#B - g]);
end while;
return B;

end function;


function RREvaluations(X, d, P)
// Gives a basis for pole order up to d at a non-Weierstrass point, evaluated

g := X`g; gens := X`RRgens; KA := X`KA; F := Parent(P[1]);
if d le g then
    return [ F ! 1 ];
end if;
if d le (2*g + 1) then
    return [ F ! 1 ] cat ([ Evaluate(KA ! gen, P) : gen in gens ])[1..(d - g)];
end if;
evs := [ F ! 1 ] cat [ Evaluate(KA ! gen, P) : gen in gens ];
while #evs lt (1 - g + d) do
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
xs := [ X`KA ! b : b in RRBasis(X, d + X`g + 1) ];
ys := [ Y`KA ! b : b in RRBasis(Y, 2*(Y`g) + 1) ];
return [ hX(x)*hY(y) : x in xs, y in ys ];

end function;


function ProductEvaluations(X, Y, d, P, Qs)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Evaluations for divisors of degree d coming from the ambient of X.
 */

xs := RREvaluations(X, d + X`g + 1, P);
yss := [ RREvaluations(Y, 2*(Y`g) + 1, Q) : Q in Qs ];
return [ [ x*y : x in xs, y in ys ] : ys in yss ];

end function;


function GlobalGenerators(X)
// Find equations in affine space from vector with elements of fraction field

if assigned X`globgens and assigned X`DEs_sub then
    return X`globgens, X`DEs_sub;
end if;

// First factor of product
if X`is_hyperelliptic then
    P := X`P0; g := X`g;
    RA<x,y> := X`RA; SA<u,v> := PolynomialRing(X`F, 2);
    u_sub := (1 / (x - P[1])); v_sub := y / (x - P[1])^(g + 1);
    x_sub := (1 / u) + P[1];   y_sub := v / u^(g + 1);
    DEs_sub := [ SA ! (u^(2*g + 2) * Evaluate(DE, [x_sub, y_sub])) : DE in X`DEs ];
    X_sub := Curve(AffineSpace(SA), DEs_sub);
    RX_sub := CoordinateRing(X_sub); KX_sub := FunctionField(X_sub);
    RA_sub := CoordinateRing(Ambient(X_sub)); KA_sub := FieldOfFractions(RA_sub);
    gens_sub := [ KX_sub ! Evaluate(X`KA ! gen, [x_sub, y_sub]) : gen in X`RRgens ];
    /* NOTE: We really need a two-step coercion */
    gens_sub := [ RA_sub ! KA_sub ! gen : gen in gens_sub ];
else
    P := X`P0; g := X`g;
    RA<x,y> := X`RA; SA<u,v> := PolynomialRing(X`F, 2);
    u_sub := (1 / (x - P[1])); v_sub := y / (x - P[1]);
    x_sub := (1 / u) + P[1];   y_sub := v / u;
    DEs_sub := [ SA ! (u^Degree(DE) * Evaluate(DE, [x_sub, y_sub])) : DE in X`DEs ];
    X_sub := Curve(AffineSpace(SA), DEs_sub);
    RX_sub := CoordinateRing(X_sub);
    KX_sub := FunctionField(X_sub);
    RA_sub := CoordinateRing(Ambient(X_sub));
    KA_sub := FieldOfFractions(RA_sub);
    gens_sub := [ KX_sub ! Evaluate(X`KA ! gen, [x_sub, y_sub]) : gen in X`RRgens ];
    /* NOTE: We really need a two-step coercion */
    //print P;
    //print [ KA_sub ! gen : gen in gens_sub ];
    gens_sub := [ RA_sub ! KA_sub ! gen : gen in gens_sub ];
end if;
X`globgens := gens_sub;
X`DEs_sub := DEs_sub;
return X`globgens, X`DEs_sub;

end function;


function GlobalBasis(X, d)

g := X`g; gens := X`globgens; RA := X`RA;
if d le g then
    return [ RA ! 1 ];
end if;
if d le (2*g + 1) then
    return [ RA ! 1 ] cat ([ RA ! gen : gen in gens ])[1..(d - g)];
end if;
B := [ RA ! 1 ] cat [ RA ! gen : gen in gens ];
while #B lt (1 - g + d) do
    Append(~B, B[2]*B[#B - g]);
end while;
return B;

end function;


function GlobalProductBasis(X, Y, d)

globX := X`globgens; globY := Y`globgens;
hX, hY := ExtractHomomorphismsRing(X, Y);
xs := [ X`RA ! b : b in GlobalBasis(X, d + X`g + 1) ];
ys := [ Y`RA ! b : b in GlobalBasis(Y, 2*(Y`g) + 1) ];
return [ hX(x)*hY(y) : x in xs, y in ys ];

end function;
