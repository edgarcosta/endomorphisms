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

intrinsic TransferMatrices(As::SeqEnum, K::Fld, L::Fld) -> SeqEnum
{Transforms the matrices a under some isomorphism from K to L.}

L<s> := L;
test, h := IsIsomorphic(K, L);
return [ Matrix([ [ h(c) : c in Eltseq(row) ] : row in Rows(A) ]) : A in As ];

end intrinsic;
