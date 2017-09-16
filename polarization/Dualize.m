/***
 *  Dualizing a morphism between PPAVs
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Hanselman  (jeroen.hanselman@uni-ulm.de)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

function Dualization(R, EL, EM);

rowsEM := Rows(EM); rowsEL := Rows(EL); rowsR := Rows(R);
rowsDual := [ ];
for rowEM in rowsEM do
    /* Use values of EM when fixing basis element of M on the left: */
    v := Matrix(rowEM);
    /* Pair those with the images under R of the basis of L: */
    w := Matrix([[ (v*Transpose(Matrix(rowR)))[1,1] : rowR in rowsR ]]);
    /* Recognize the coresponding functional" */
    Append(~rowsDual, Eltseq(Solution(EL, w)));
end for;
return Matrix(rowsDual);

end function;
