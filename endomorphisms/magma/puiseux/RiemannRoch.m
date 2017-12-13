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


forward RRGenerators;
forward RRBasis;
forward RREvaluate;


function RRGenerators(X)
// Gives all Riemann-Roch functions needed to determine the others as monomials
// in them

g := X`g; KU := X`KU; P0 := X`P0;
V, phiV := RiemannRochSpace((2*g + 1)*Divisor(P0));
return [ KU ! phiV(V.i) : i in [1..(Dimension(V) - 1)] ];

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


function RREvaluate(X, P, d)
// Gives a basis for pole order up to d at a non-Weierstrass point, evaluated

g := X`g; gens := X`RRgens; KA := X`KA; F := Parent(P[1]);
if d le g then
    return [ F ! 1 ];
end if;
if d le (2*g + 1) then
    return [ F ! 1 ] cat [ Evaluate(KA ! gen, P) : gen in gens ];
end if;
evs := [ F ! 1 ] cat [ Evaluate(KA ! gen, P) : gen in gens ];
while #evs lt (1 - g + d) do
    Append(~evs, evs[2]*evs[#evs - g]);
end while;
return evs;

end function;
