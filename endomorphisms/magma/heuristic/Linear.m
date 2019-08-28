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
{Returns the numerical solution X to the equation X * A = B.}

//return B * A^(-1);
R := BaseRing(A);
return NumericalSolution(A, B : Epsilon := R`epscomp);

end intrinsic;


intrinsic NumericalRightSolve(A::., B::.) -> .
{Returns the numerical solution X to the equation A * X = B.}

return Transpose(NumericalLeftSolve(Transpose(A), Transpose(B)));

end intrinsic;


intrinsic SubmatrixOfRank(M::., rk::RngIntElt : ColumnsOrRows := "Columns") -> .
{Returns a submatrix of M of rank rk, together with the corresponding list of
rows or columns. Returns an error if such a matrix does not seem to exist.
ColumnsOrRows specifies which of the two are culled down.}
/* TODO: Use an LU decomposition instead. However, no numerical version of that
 * seems to work for now. */

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


intrinsic IntegralLeftKernel(M::. : OneRow := false, CalcEndoRep := false) -> .
{Returns simultaneous integral cancellations of all the rows of M.}

RR := BaseRing(M);
eps := Minimum([ Abs(c) : c in Eltseq(M) | not Abs(c) lt RR`epscomp ]);
//eps := 1;
if CalcEndoRep then
    B := 10^12 / eps;
else
    B := 10^(Precision(RR) + 10) / eps;
end if;
MJ := Matrix(Integers(), [ [ Round(B * c) : c in Eltseq(row) ] : row in Rows(M) ]);
MI := IdentityMatrix(Integers(), #Rows(MJ));
MJ := HorizontalJoin(MI, MJ);

L, K := LLL(MJ);
rowsK := Rows(K);
if OneRow then
    rowsK := [ rowsK[1] ];
end if;

CCSmall := ComplexField(5);
rowsK0 := [ ];
for row in rowsK do
    print "";
    print row;
    ht := Max([ Abs(c) : c in Eltseq(row) ]);
    test1 := ht lt RR`height_bound;
    //test1 := true;

    if test1 then
        prod := Matrix(RR, [ Eltseq(row) ])*M;
        abs := Max([ Abs(c) : c in Eltseq(prod) ]);
        test2 := abs lt RR`epscomp;
        //test2 := abs lt 10^50*RR`epscomp;

        if test2 then
            vprint EndoFind, 3 : "";
            vprint EndoFind, 3 : "Height of row:", Round(Log(ht));
            vprint EndoFind, 3 : "Precision reached:", CCSmall ! abs;
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


intrinsic MatrixInBasis(M::., Bs::SeqEnum) -> .
{Returns a vector that describes M as a rational combination of the elements in Bs.}

if #Bs eq 0 then
    return false, 0;
end if;
MBs := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(b) ] : b in Eltseq(B) ] : B in Bs ]);
MM := Matrix(Rationals(), [ &cat[ &cat[ Eltseq(c) : c in Eltseq(m) ] : m in Eltseq(M) ] ]);
return IsConsistent(MBs, MM);

end intrinsic;


intrinsic IsRationalMultiple(v::., v0::.) -> .
{Returns whether v is a rational multiple of v0 or not, along with a corresponding scalar if this is the case.}

M := Matrix([ Eltseq(v) ]);
M0 := Matrix([ Eltseq(v0) ]);
test, tup := MatrixInBasis(M, [ M0 ]);
if not test then
    return false, 0;
end if;
return test, Eltseq(tup)[1];

end intrinsic;


intrinsic IsMultiplePolynomial(f::., f0::.) -> .
{Returns whether f is a multiple of f0 or not, along with a corresponding scalar if this is the case.}

if not Monomials(f) eq Monomials(f0) then
    return false, 0;
end if;
coeffs := Coefficients(f); coeffs0 := Coefficients(f0);
M := Matrix([ coeffs ]); M0 := Matrix([ coeffs0 ]);
return IsConsistent(M, M0);

end intrinsic;
