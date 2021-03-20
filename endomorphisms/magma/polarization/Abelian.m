/***
 *  Abelian category functionality
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

import "../heuristic/Saturate.m": SaturateLattice;


intrinsic KerModKer0(h::., P::., Q::.) -> .
{Returns the group of connected components of the kernel of the morphism h = (A, R) from P to Q.}

/* Check input */
A := h[1]; R := h[2]; CC := BaseRing(A);
assert Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]) le 10^20*CC`epscomp;
assert #Rows(R) eq #Rows(Transpose(Q)); assert #Rows(Transpose(R)) eq #Rows(Transpose(P));
test := &and[ IsIntegral(c) : c in Eltseq(R) ];
if not test then
    error "Homology representation not integral";
end if;

L := Lattice(Transpose(R));
return quo< PureLattice(L) | L >;

end intrinsic;


intrinsic Ker0(h::., P::., Q::.) -> .
{Returns the connected component of the kernel of the morphism h = (A, R) from P to Q.}

/* Check input */
A := h[1]; R := h[2]; CC := BaseRing(A);
assert Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]) le 10^20*CC`epscomp;
assert #Rows(R) eq #Rows(Transpose(Q)); assert #Rows(Transpose(R)) eq #Rows(Transpose(P));

/* Linear algebra on homology */
RMS := RMatrixSpace(Integers(), #Rows(R), #Rows(Transpose(R)));
R0 := RMS ! (Denominator(R)*R); B := Basis(Kernel(Transpose(R0)));
if #B eq 0 then
    return 0, [0, 0];
end if;
B := [ Eltseq(b) : b in B ];

/* Build corresponding matrix */
colsP := Rows(Transpose(P)); colsK := [ &+[ b[i]*colsP[i] : i in [1..#b] ] : b in B ];
K := Transpose(Matrix(colsK)); K, s := SubmatrixOfRank(K, #B div 2 : ColumnsOrRows := "Rows");
ihom := Transpose(Matrix(Rationals(), [ [ c : c in Eltseq(row) ] : row in B ]));
/* TODO: Think of a better way */
itan := TangentRepresentation(ihom, K, P);

/* Check output */
assert IsZero(R*ihom);
L := Lattice(Transpose(ihom));
assert L eq PureLattice(L);
test := Maximum([ Abs(c) : c in Eltseq(itan*K - P*ChangeRing(ihom, CC)) ]);
if test gt 10^20*CC`epscomp then
    error "Error in determining tangent representation:", ComplexField(5) ! test;
end if;

return K, [* itan, ihom *];

end intrinsic;


function QuotientByColumns(A)

CC := BaseRing(A);
r := #Rows(A); c := #Rows(Transpose(A));
V := VectorSpace(CC, r); VmodW := VectorSpace(CC, r - c);

/* Standardize */
A0, s0 := SubmatrixOfRank(A, c : ColumnsOrRows := "Rows");
A := A*A0^(-1);
s := s0 cat [ i : i in [1..r] | not i in s0 ];
P := PermutationMatrix(CC, s);
A := P*A;

/* Check */
I := Submatrix(A, 1,1,   c,c);
B := Submatrix(A, c+1,1, r - c,c);
test := Maximum([ Abs(entry) : entry in Eltseq(I - IdentityMatrix(CC, c)) ]);
if test gt 10^20*CC`epscomp then
    error "QuotientByColumns does not give identity matrix";
end if;

/* Projection morphism */
Q := HorizontalJoin(-B, IdentityMatrix(CC, r - c));
Q := Q*P^(-1);
q := hom< V -> VmodW | Transpose(Q) >;
return V, VmodW, q;

end function;


intrinsic Coker(h::., P::., Q::.) -> .
{Returns the cokernel of the morphism h = (A, R) from P to Q.}

/* Check input */
A := h[1]; R := h[2]; CC := BaseRing(A);
assert Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]) le 10^20*CC`epscomp;
assert #Rows(R) eq #Rows(Transpose(Q)); assert #Rows(Transpose(R)) eq #Rows(Transpose(P));

/* Build lattice and complement */
L := PureLattice(Lattice(Transpose(R))); Zd := StandardLattice(Degree(L));
CL, projL := quo< Zd | L >; gensCL := [ gen : gen in Generators(CL) ];
if #gensCL eq 0 then
    return 0, [0, 0];
end if;

/* Build new tangent matrix surjecting onto the appropriate number of columns */
B := [ Zd ! b : b in Basis(L) ]; BC := [ Zd ! (g @@ projL) : g in gensCL ];
colsQ := Rows(Transpose(Q)); colsC := [ &+[ b[i]*colsQ[i] : i in [1..#Eltseq(b)] ] : b in BC ];
A0, s := SubmatrixOfRank(A, #B div 2 : ColumnsOrRows := "Columns");

/* Build period matrix */
V, VmodW, q := QuotientByColumns(A0);
C := Transpose(Matrix([ Eltseq(q(V ! col)) : col in colsC ]));

/* Correction of homology representation: */
M := Matrix(Rationals(), [ Eltseq(b) : b in B cat BC ]);
phom := VerticalJoin(ZeroMatrix(Rationals(), #B, #BC), IdentityMatrix(Rationals(), #BC));
phom := Transpose(M^(-1)*phom);
ptan := TangentRepresentation(phom, Q, C);

/* Check output */
assert IsZero(phom*R);
test := Maximum([ Abs(c) : c in Eltseq(ptan*Q - C*ChangeRing(phom, CC)) ]);
if test gt 10^20*CC`epscomp then
    error "Error in determining tangent representation:", ComplexField(5) ! test;
end if;

return C, [* ptan, phom *];

end intrinsic;


intrinsic ImgInc(h::., P::., Q::.) -> .
{Returns the image of the morphism h = (A, R) from P to Q. This comes with an inclusion map into Q.}

/* Check input */
A := h[1]; R := h[2]; CC := BaseRing(A);
assert Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]) le 10^20*CC`epscomp;

/* Compose previous algorithms */
C, proj := Coker(h, P, Q);
if C eq 0 then
    A := IdentityMatrix(BaseRing(Q), #Rows(Q));
    R := IdentityMatrix(Rationals(), #Rows(Transpose(Q)));
    return Q, [* A, R *];
end if;

K, i := Ker0(proj, Q, C); itan, ihom := Explode(i);
/* Check output */
test := Maximum([ Abs(c) : c in Eltseq(itan*K - Q*ChangeRing(ihom, CC)) ]);
if test gt 10^20*CC`epscomp then
    error "Error in determining tangent representation:", ComplexField(5) ! test;
end if;
return K, i;

end intrinsic;


intrinsic ImgProj(h::., P::., Q::.) -> .
{Returns the image of the morphism h = (A, R) from P to Q. This comes with a projection map from P.}

/* Check input */
A := h[1]; R := h[2]; CC := BaseRing(A);
assert Maximum([ Abs(c) : c in Eltseq(A*P - Q*ChangeRing(R, CC)) ]) le 10^20*CC`epscomp;

/* Compose previous algorithms */
K, i := Ker0(h, P, Q);
if K eq 0 then
    A := IdentityMatrix(BaseRing(P), #Rows(P));
    R := IdentityMatrix(Rationals(), #Rows(Transpose(P)));
    return P, [* A, R *];
end if;

C, p := Coker(i, K, P); ptan, phom := Explode(p);
/* Check output */
test := Maximum([ Abs(c) : c in Eltseq(ptan*P - C*ChangeRing(phom, CC)) ]);
if test gt 10^20*CC`epscomp then
    error "Error in determining tangent representation:", ComplexField(5) ! test;
end if;
return C, p;

end intrinsic;


intrinsic ImgIdemp(idem::List, P::.) -> List
{Given an idempotent idem for the period matrix P, returns a corresponding
lattice and an analytic representation of the projection to it.}

/* g is the dimension of the new abelian variety */
A, R := Explode(idem); g := Rank(R) div 2;
ACC := TangentRepresentation(R, P);
CC := BaseRing(ACC); RR := RealField(CC);
/* At this point A P = P R */

/* Determine rows (corresponds to taking full-rank projection on components): */
QLarge, indices := SubmatrixOfRank(ACC * P, g : ColumnsOrRows := "Rows");
rowsA := Rows(A); rowsACC := Rows(ACC);
B := Matrix([ [ c : c in Eltseq(rowsA[i]) ] : i in indices ]);
BCC := Matrix([ [ c : c in Eltseq(rowsACC[i]) ] : i in indices ]);
/* At this point B P = QLarge R */

/* Determine columns (corresponds to finding generating elements of the lattice */
QLargeSplit := VerticalSplitMatrix(QLarge);
QSubSplit, indices := SubmatrixOfRank(QLargeSplit, 2*g : ColumnsOrRows := "Columns");
QSplit, T, U := SaturateLattice(QSubSplit, QLargeSplit : ColumnsOrRows := "Columns");
Q := CombineVerticallySplitMatrix(QSplit, CC);
/* We now have QLarge = Q U, so B P = Q U R. We return B and U R. */

S := U * R;
proj := [* B, S, BCC *];
/* Check output */
if Maximum([ Abs(c) : c in Eltseq(BCC*P - Q*ChangeRing(S, CC)) ]) gt 10^20*CC`epscomp then
    error "Error in determining projection";
end if;
proj := [* B, S *];
return Q, proj;

end intrinsic;
