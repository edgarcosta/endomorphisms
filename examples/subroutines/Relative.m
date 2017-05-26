/***
 *  Examples of relative splitting fields
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

AttachSpec("../../spec");

print "Over QQ:";

F := Rationals();
R<x> := PolynomialRing(F);
f := x^3 + x + 1;
g := x^2 - 5;

K := NumberFieldExtra(f);
SetInfinitePlace(K, InfinitePlaces(K)[2]);

K := NumberField(f);
L := RelativeSplittingFieldExtra([f, g]);
print L;
print L`iota;
//print Roots(f, L);

print "Over a number field:";

F<r> := NumberField(x^2 - 2);
R<x> := PolynomialRing(F);
f := x^3 + x - r;
g := x^2 - 5;
g := x^2 - 5 + r;

K := NumberField(f);
L := RelativeSplittingFieldExtra([f, g]);
print L;
print L`iota;
//print Roots(f, L);

SetInfinitePlace(L, InfinitePlaces(L)[7]);
print L`iota;
print F`iota;

F := Rationals();
R<x> := PolynomialRing(F);
F<r> := NumberField(x^2 - 2);
R<x> := PolynomialRing(F);
f := x^3 + x + r;
time K := SplittingField(f);
print K;
time K := RelativeSplittingFieldExtra(f);
print K;
print K`iota;

RestrictInfinitePlace(K, K);
print K`iota;

exit;
