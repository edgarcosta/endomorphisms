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


forward VariableOrder;
forward ExtractHomomorphismsRing;
forward ExtractHomomorphismsField;
forward ExtractPoints;


function VariableOrder()
/*
 * The order in which x(P), y(P), x(Q), y(Q) and hence x1, y1, x2, y2 are used
 * in the product space. Note that because of the lexicographical ordering
 * variables that occur later are eliminated for first.
 */

/* x(P) to 4th comp, y(P) to 2nd comp, etc */
// TODO: Test other ones
return [4, 2, 3, 1];

end function;


function ExtractHomomorphismsRing(X, Y)

RAX := X`RA; RAY := Y`RA;
varord := VariableOrder();
// TODO: Test other orderings
RAXY := PolynomialRing(X`F, 4, "grevlex");
seqX := [ RAXY.varord[i] : i in [1..2] ];
seqY := [ RAXY.varord[i] : i in [3..4] ];
hX := hom< RAX -> RAXY | seqX >;
hY := hom< RAY -> RAXY | seqY >;
return hX, hY;

end function;


function ExtractHomomorphismsField(X, Y)

KAX := X`KA; KAY := Y`KA;
varord := VariableOrder();
// TODO: Test other orderings
RAXY := PolynomialRing(X`F, 4, "grevlex");
KAXY := FieldOfFractions(RAXY);
seqX := [ KAXY.varord[i] : i in [1..2] ];
seqY := [ KAXY.varord[i] : i in [3..4] ];
hX := hom< KAX -> KAXY | seqX >;
hY := hom< KAY -> KAXY | seqY >;
return hX, hY;

end function;


function ExtractPoints(X, Y, P, Q)
/* Reflects order in VariableOrder */

seq := [ P[1], P[2], Q[1], Q[2] ];
varord := VariableOrder();
return [ seq[Index(varord, i)] : i in [1..#varord] ];

end function;
