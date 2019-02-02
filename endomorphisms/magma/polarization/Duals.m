/***
 *  Dualizing a morphism between PPAVs
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


/* TODO: Still very incomplete. For the right notion of dual, also outside the
 * PPAV case, we need to work directly with the definition via antilinear maps.
 * In the PPAV case this can be transferred via the polarization. */

intrinsic Dualization(R::., EL::AlgMatElt, EM::AlgMatElt) -> .
{Dualizes a homology representation with respect to the alternating pairings EL
and EM on the source and target.}

rowsEM := Rows(EM); rowsEL := Rows(EL); rowsR := Rows(R);
rowsDual := [ ];
for rowEM in rowsEM do
    /* Use values of EM when fixing basis element of M on the left: */
    v := Matrix(rowEM);
    /* Pair those with the images under R of the basis of L: */
    w := Matrix([[ (v*Transpose(Matrix(rowR)))[1,1] : rowR in rowsR ]]);
    /* Recognize the corresponding functional" */
    Append(~rowsDual, Eltseq(Solution(EL, w)));
end for;
return Matrix(rowsDual);

end intrinsic;
