/***
 *  Verifies that endomorphism ring is saturated
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

// TODO: Generalize to isogenies

intrinsic VerifySaturated(GeoEndoRep::SeqEnum, P::.) -> BoolElt, .
{Returns a boolean that indicates whether the endomorphism ring in GeoEndoRep
is saturated in the corresponding algebra.}

genHoms := [ gen[2] : gen in GeoEndoRep ];
// Creation of relevant algebras
g := #Rows(genHoms[1]) div 2;
// Ambient matrix algebra, plus generators of the endomorphism ring
A := Algebra(MatrixRing(Rationals(), 2*g)); GensA := [ A ! Eltseq(genHom) : genHom in genHoms ];
// As a subalgebra
B := sub<A | GensA>; GensB := [ B ! gen : gen in GensA ];
// As an associative algebra
C := AssociativeAlgebra(B); GensC := [ C ! gen : gen in GensB ];

OC := Order(Integers(), GensC);
DOC := Discriminant(OC); DOM := Discriminant(MaximalOrder(C));
test, ind := IsSquare(DOC / DOM); ind := Integers() ! ind;

ps := [ tup[1] : tup in Factorization(ind) ];
for p in ps do
    test, genHom := VerifySaturatedAtPrime(GeoEndoRep, P, p);
    if not test then
        return false, genHom;
    end if;
end for;
return true, "";

end intrinsic;


intrinsic VerifySaturatedAtPrime(GeoEndoRep::SeqEnum, P::., p::RngIntElt) -> BoolElt, .
{Returns a boolean that indicates whether the endomorphism ring in GeoEndoRep
is saturated in the corresponding algebra at p.}

/* Uses an extremely naive algorithm by excluding every intermediate element
 * directly */
genHoms := [ gen[2] : gen in GeoEndoRep ];
CC := Parent(P[1,1]); RR := RealField(CC); JP := ComplexStructure(P);
I := [0..(p - 1)]; d := #genHoms; CP := CartesianPower(I, d);

for tup in [ tup : tup in CP | not &and[ c eq 0 : c in tup ] ] do
    genHom := &+[ tup[i]*genHoms[i] : i in [1..#genHoms] ]/p;
    if &and[ c in Integers() : c in Eltseq(genHom) ] then
        D := JP*genHom - genHom*JP;
        if &and[ (RR ! Abs(d)) lt RR`epscomp : d in Eltseq(D) ] then
            return false, genHom;
        end if;
    end if;
end for;
return true, "";

end intrinsic;
