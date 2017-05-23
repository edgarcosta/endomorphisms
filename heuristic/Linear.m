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
{Solves X * A = B.}

// TODO: NumericalKernel should be used.
return B * A^(-1);

end intrinsic;


intrinsic SubmatrixOfRank(M::., rk::RngIntElt : ColumnsOrRows := "Rows") -> .
{Finds a submatrix of M of rank rk, returning an error if such a matrix does not seem to exist. ColumnsOrRows specifies which are used.}

// Reducing to columns
if ColumnsOrRows eq "Columns" then
    M := Transpose(M);
end if;

// Elementary invariants
CC := BaseRing(M);
r := #Rows(M); c := #Rows(Transpose(M));

RM := Rows(M);
for s in Subsets({1..r}, rk) do
    s0 := s;
    M0 := Matrix([ RM[i] : i in s ]);
    if NumericalRank(M0 : Epsilon := CC`epsinv) eq rk then
        if ColumnsOrRows eq "Rows" then
            return M0, s0;
        else
            return Transpose(M0), s0;
        end if;
    end if;
end for;
error Error("Failed to find submatrix of the desired rank");

end intrinsic;


intrinsic InvertibleSubmatrix(M::.) -> .
{Finds an invertible submatrix of M of maximal possible dimension, returning an error if such a matrix does not seem to exist.}

return SubmatrixOfRank(M, Minimum(#Rows(M), #Rows(Transpose(M))));

end intrinsic;


intrinsic SplitMatrix(M::.) -> .
{Splits M into real and imaginary part, stacked on top of each other.}

CC := BaseRing(M); RR := RealField(CC);
MSplitRe := Matrix(RR, [ [ Real(c) : c in Eltseq(r)] : r in Rows(M) ]);
MSplitIm := Matrix(RR, [ [ Im(c) : c in Eltseq(r)] : r in Rows(M) ]);
return HorizontalJoin([MSplitRe, MSplitIm]);

end intrinsic;


intrinsic CombineMatrix(MSplit::., CC::FldCom) -> .
{Combines MSplit back into a full period matrix with coefficients in CC.}

// Basic invariants
r := #Rows(MSplit); c := #Rows(Transpose(MSplit));
// Combining parts
MRe := Matrix(CC, Submatrix(MSplit, [1..r], [1..(c div 2)]));
MIm := Matrix(CC, Submatrix(MSplit, [1..r], [(c div 2)+1..c]));
return MRe + CC.1*MIm;

end intrinsic;


intrinsic IntegralLeftKernel(M::.) -> .
{Returns simultaneous integral cancellations of all the rows of M.}
// TODO: Integral splitting here?

RR := BaseRing(M);
MI := IdentityMatrix(RR, #Rows(M));
MJ := HorizontalJoin(MI, (1 / RR`epsLLL) * M);
L, K := LLL(MJ);
return K;

end intrinsic;


intrinsic ConjugateMatrix(sigma::Map, M::.) -> .
{Transforms the matrix M by the field automorphism sigma.}

return Matrix([ [ sigma(elt) : elt in Eltseq(row) ] : row in Rows(M) ]);

end intrinsic;


intrinsic MatrixInBasis(M::., Bs::SeqEnum) -> .
{Writes M as a combination of the elements of B over the rationals, provided that we are only two levels above.}

MBs := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(b) ] : b in Eltseq(B) ] : B in Bs ]);
MM := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(m) ] : m in Eltseq(M) ] ]);
return Matrix(Solution(MBs, MM));

end intrinsic;


intrinsic SaturateLattice(FullLattice::., SubLattice::.) -> .
{Input: Matrices FullLattice, SubLattice whose rows generate lattices such that the latter has finite index in the former. Output: A matrix of the same size as SubLattice whose rows generate the same lattice as FullLattice.}

MatrixSize := Min(#Rows(SubLattice), #Rows(Transpose(FullLattice)));

SaturationMatrix, SatColumns := SubmatrixOfRank(SubLattice, MatrixSize : ColumnsOrRows := "Columns");
LMat := Matrix(Rationals(), 0,MatrixSize, []);
for row in Rows(FullLattice) do
    linearDependence := NumericalLeftSolve(SaturationMatrix, Matrix([[row[j] : j in SatColumns]]) );
    linearDependence := Matrix([ [ FractionalApproximation(c) : c in Eltseq(row) ] : row in Rows(linearDependence) ]);
    LMat := VerticalJoin(LMat, linearDependence);
end for;

M := Matrix(Basis(Lattice(LMat)));
M := Matrix(BaseRing(FullLattice),M);
return M*SubLattice, M;

end intrinsic;
