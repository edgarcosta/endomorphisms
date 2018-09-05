/***
 *  Linear algebra functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic NumericalLeftSolve(A::., B::.) -> .
{Returns the solution X to the equation X * A = B.}

//return B * A^(-1);
R := BaseRing(A);
return NumericalSolution(A, B : Epsilon := R`epscomp);

end intrinsic;


intrinsic NumericalRightSolve(A::., B::.) -> .
{Returns the solution X to the equation A * X = B.}

return Transpose(NumericalLeftSolve(Transpose(A), Transpose(B)));

end intrinsic;


intrinsic SubmatrixOfRank(M::., rk::RngIntElt : ColumnsOrRows := "Columns") -> .
{Returns a submatrix of M of rank rk, together with the corresponding list of
rows or columns. Returns an error if such a matrix does not seem to exist.
ColumnsOrRows specifies whether columns of rows are used.}
/* TODO: Use an LU decomposition instead */

/* Reducing to the case of columns */
if ColumnsOrRows eq "Columns" then
    N, s0 := SubmatrixOfRank(Transpose(M), rk : ColumnsOrRows := "Rows");
    return Transpose(N), s0;
end if;

/* Elementary invariants */
CC := BaseRing(M); RM := Rows(M);
r := #RM; c := #Rows(Transpose(M));

/* Prefer obvious choice if possible */
s0 := [ 1..rk ];
M0 := Matrix([ RM[i] : i in s0 ]);
if NumericalRank(M0 : Epsilon := CC`epsinv) eq rk then
    return M0, s0;
end if;

for s in Subsets({1..r}, rk) do
    s0 := [ c : c in s ];
    M0 := Matrix([ RM[i] : i in s ]);
    if NumericalRank(M0 : Epsilon := CC`epsinv) eq rk then
        return M0, s0;
    end if;
end for;
error "Failed to find submatrix of the desired rank";

end intrinsic;


intrinsic InvertibleSubmatrix(M::. : IsPeriodMatrix := false) -> .
{Returns an invertible submatrix of M. We can indicate that M is a period
matrix with respect to a symplectic basis, which trivially speeds up the
calculation by extracting the first m columns.}

r := #Rows(M); c := #Rows(Transpose(M)); m := Min(r, c);
/* Speedup for period matrices with respect to a symplectic basis */
if IsPeriodMatrix then
    return Submatrix(M, [1..m], [1..m]), [1..m];
end if;
if r gt c then
    return SubmatrixOfRank(M, m : ColumnsOrRows := "Rows");
else
    return SubmatrixOfRank(M, m : ColumnsOrRows := "Columns");
end if;

end intrinsic;


intrinsic HorizontalSplitMatrix(M::.) -> .
{Returns the horizontal join of the real and imaginary part of M, so (Re | Im).}

CC := BaseRing(M); RR := RealField(CC);
MSplitRe := Matrix(RR, [ [ Real(c) : c in Eltseq(r)] : r in Rows(M) ]);
MSplitIm := Matrix(RR, [ [ Im(c) : c in Eltseq(r)] : r in Rows(M) ]);
return HorizontalJoin([ MSplitRe, MSplitIm ]);

end intrinsic;


intrinsic VerticalSplitMatrix(M::.) -> .
{Returns the vertical join of the real and imaginary part of M, so Re over Im.}

CC := BaseRing(M); RR := RealField(CC);
MSplitRe := Matrix(RR, [ [ Real(c) : c in Eltseq(r)] : r in Rows(M) ]);
MSplitIm := Matrix(RR, [ [ Im(c) : c in Eltseq(r)] : r in Rows(M) ]);
return VerticalJoin([ MSplitRe, MSplitIm ]);

end intrinsic;


intrinsic CombineVerticallySplitMatrix(MSplit::., CC::FldCom) -> .
{Returns the combination of the vertically split matrix MSplit back into a full
period matrix over CC. So Re over Im becomes Re + i*Im.}

r := #Rows(MSplit); c := #Rows(Transpose(MSplit));
MRe := Matrix(CC, Submatrix(MSplit, [1..(r div 2)],   [1..c]));
MIm := Matrix(CC, Submatrix(MSplit, [(r div 2)+1..r], [1..c]));
return MRe + CC.1*MIm;

end intrinsic;


intrinsic IntegralLeftKernel(M::. : OneRow := false, EndoRep := false) -> .
{Returns simultaneous integral cancellations of all the rows of M.}

RR := BaseRing(M);
MI := IdentityMatrix(RR, #Rows(M));
if EndoRep then
    eps := Minimum([ Abs(c) : c in Eltseq(M) | not Abs(c) lt RR`epscomp ]);
    MJ := HorizontalJoin(MI, (10^12 / eps) * M);
else
    MJ := HorizontalJoin(MI, (10^(Precision(RR) - 30)) * M);
end if;
MJ := Matrix(Integers(), [ [ Round(c) : c in Eltseq(row) ] : row in Rows(MJ) ]);

L, K := LLL(MJ);
rowsK := Rows(K); rowsK0 := [ ];
if OneRow then
    rowsK := [ rowsK[1] ];
end if;

for row in rowsK do
    test1 := &and[ Abs(c) lt RR`height_bound : c in Eltseq(row) ];
    if test1 then
        prod := Matrix(RR, [ Eltseq(row) ])*M;
        test2 := &and[ Abs(c) lt RR`epscomp : c in Eltseq(prod) ];
        if test2 then
            Append(~rowsK0, row);
        end if;
    end if;
end for;
if #rowsK0 eq 0 then
    return Matrix([ [ 0 : row in Rows(M) ] ]), false;
else
    return Matrix(rowsK0), true;
end if;

end intrinsic;


intrinsic IntegralRightKernel(M::.) -> .
{Returns simultaneous integral cancellations of all the columns of M.}

K, test := IntegralLeftKernel(Transpose(M));
return Transpose(K), test;

end intrinsic;


intrinsic ConjugateMatrix(sigma::Map, M::.) -> .
{Returns the transformation of the matrix M by the field automorphism sigma.}

return Matrix([ [ sigma(elt) : elt in Eltseq(row) ] : row in Rows(M) ]);

end intrinsic;


intrinsic MatrixInBasis(M::., Bs::SeqEnum) -> .
{Returns a vector that describes M as a rational combination of the elements in
Bs. Assume that base field of M is at most a double extension of the field of
rationals.}

MBs := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(b) ] : b in Eltseq(B) ] : B in Bs ]);
MM := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(m) ] : m in Eltseq(M) ] ]);
return Matrix(Solution(MBs, MM));

end intrinsic;


intrinsic SaturateLattice(L::., M::. : ColumnsOrRows := "Columns") -> .
{Given a basis of a lattice L and a generating set of a lattice M in which L is
of finite index, returns a basis of M along with matrices that give expressions
of the provided generating sets in this basis. The flag ColumnsOrRows specifies
whether column of row vectors are interpreted as generating the lattice.}

/* In the end we have L = T B, M = U B in case of Rows and L = B T, M = B U in
 * case of Columns. */

if ColumnsOrRows eq "Columns" then
    B, T, U := SaturateLattice(Transpose(L), Transpose(M) : ColumnsOrRows := "Rows");
    return Transpose(B), Transpose(T), Transpose(U);
end if;

CC := BaseRing(L);
subL, s0 := InvertibleSubmatrix(L);
subM := Submatrix(M, [1..#Rows(M)], s0);
S := NumericalLeftSolve(subL, subM);
S, test := FractionalApproximationMatrix(S);
if not test then
    error "No suitable fractional approximation found";
end if;
/* At this point we have S L = M, where S has rational entries and an integral
* inverse */

/* Now we write S = R S0, where R is integral and where S0 has an integral
 * inverse */
S0 := ChangeRing(Matrix(Basis(Lattice(S))), Rationals());
S0i := S0^(-1); R := S * S0i;
/* The result is that the following B, T, U can be used */
B := ChangeRing(S0, CC) * L; T := S0i; U := R;

/* Final sanity check */
test1 := Minimum([ Abs(c) : c in Eltseq(L - ChangeRing(T, CC)*B) ]) lt CC`epscomp;
test2 := Minimum([ Abs(c) : c in Eltseq(M - ChangeRing(U, CC)*B) ]) lt CC`epscomp;
if not (test1 and test2) then
    error "Error in determining saturated lattice";
end if;
return B, T, U;

end intrinsic;
