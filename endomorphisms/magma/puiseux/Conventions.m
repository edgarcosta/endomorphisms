/***
 *  Conventions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


forward ExtractHomomorphismsRing;
forward ExtractHomomorphismsField;
forward ExtractPoints;


/*
 * The order in which x(P), y(P), x(Q), y(Q) and hence x1, y1, x2, y2 are used
 * in the product space. Note that because of the lexicographical ordering
 * variables that occur later are eliminated for first.
 */

/* x(P) to 4th comp, y(P) to 2nd comp, etc */
VariableOrder := [4, 2, 3, 1];
// note, converting and uncoverting from this ordering is just the transposition (1,4)
VariableOrderPerm := SymmetricGroup(4)!VariableOrder;
assert Order(VariableOrderPerm) eq 2;


function toProduct(RX, RY, RXY)
    return hom< RX -> RXY | var[1..2] >, hom< RY -> RXY | var[3..4] >
           where var := [RXY.j : j in VariableOrder];
end function;


function ExtractHomomorphismsRing(X, Y)
/* Homomorphism from factor to product corresponding to variable order */

RAX := X`RA;
RAY := Y`RA;
varord := VariableOrder;
RAXY := PolynomialRing(X`F, 4, "grevlex");
hX, hY := toProduct(RAX, RAY, RAXY);


// this only applies if X and Y are hyperelliptic
// as we are trying to get the projection to the x coordinate
/* RAX is not the actual ring in which we are working, but we need some
 * polynomial ring with two generators, so this one will do */
/*
seqxs := [ RAX ! 0 : i in [1..4] ];
seqxs[varord[1]] := RAX.1; seqxs[varord[3]] := RAX.2;
*/
seqxs := PermuteSequence([ RAX.1, 0, RAX.2,  0], VariableOrderPerm); // [0, 0, RAX.2, RAX.4];
seqxsinv := PermuteSequence([RAXY.i : i in [1..4]], VariableOrderPerm)[1 .. 4 by 2]; // [RAXY.4, RAXY.3];
seqxsinv := [ RAXY.varord[i] : i in [1,3] ]; // RAXY.4, RAXY.3
hxs := hom< RAXY -> RAX | seqxs >; hxsinv := hom< RAX -> RAXY | seqxsinv >;

// we dot not use y in the end
/*
seqys := [ RAX ! 0 : i in [1..4] ];
seqys[varord[1]] := RAX.1; seqys[varord[4]] := RAX.2;
seqysinv := [ RAXY.varord[i] : i in [1,4] ];
hys := hom< RAXY -> RAX | seqys >; hysinv := hom< RAX -> RAXY | seqysinv >;
*/

return hX, hY, hxs, hxsinv; // hys, hysinv;

end function;


function ExtractHomomorphismsField(X, Y)
/* Homomorphism from factor to product corresponding to variable order */
    KAX := X`KA; KAY := Y`KA;
    varord := VariableOrder;
    RAXY := PolynomialRing(X`F, 4, "grevlex");
    KAXY := FieldOfFractions(RAXY);
    return toProduct(KAX, KAY, KAXY);
end function;


function ExtractPoints(P, Q)
/* Reorders coordinates of P and Q to match that of VariableOrder */
  return [ Q[2], P[2], Q[1], P[1] ];
  return PermuteSequence(Eltseq(P) cat Eltseq(Q), VariableOrderPerm);
seq := [ P[1], P[2], Q[1], Q[2] ];
varord := VariableOrder();
// [
return [ seq[Index(varord, i)] : i in [1..#varord] ];

end function;
