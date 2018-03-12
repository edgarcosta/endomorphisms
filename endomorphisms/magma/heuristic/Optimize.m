/***
 *  Improve representation of endomorphisms
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic TransferInfinitePlace(h::Map, iotaK::.) -> .
{Return the image of the infinite place iotaK under the isomorphism h.}

/* TODO: This should not be global */
prec := 100;
K := Domain(h); L := Codomain(h);
if Type(Domain(h)) eq FldRat then
    return InfinitePlaces(L)[1];
end if;
for iotaL in InfinitePlaces(L) do
    if Abs(Evaluate(h(K.1), iotaL : Precision := prec) - Evaluate(K.1, iotaK : Precision := prec)) lt 10^(-prec + 10) then
        return iota;
    end if;
end for;

end intrinsic;


intrinsic TransferGenerator(gen::SeqEnum, K::Fld, L::Fld) -> SeqEnum
{Returns the generators in As transformed by some isomorphism from K to L.}

test, h := IsIsomorphic(K, L);
if assigned K`iota then
    L`iota := TransferInfinitePlace(h, K`iota);
end if;
genTan, genHom, genApp := Explode(gen);
genTanL := Matrix([ [ h(c) : c in Eltseq(row) ] : row in Rows(genTan) ]);
return [* genTanL, genHom, genApp *];

end intrinsic;


intrinsic TransferMatrices(As::SeqEnum, K::Fld, L::Fld) -> SeqEnum
{Returns the matrices in As transformed by some isomorphism from K to L.}

test, h := IsIsomorphic(K, L);
if assigned K`iota then
    L`iota := TransferInfinitePlace(h, K`iota);
end if;
return [ Matrix([ [ h(c) : c in Eltseq(row) ] : row in Rows(A) ]) : A in As ];

end intrinsic;
