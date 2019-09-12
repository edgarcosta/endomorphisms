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

/* Build new matrix */
B := [ Zd ! b : b in Basis(L) ]; BC := [ Zd ! (g @@ projL) : g in gensCL ];
colsQ := Rows(Transpose(Q)); colsC := [ &+[ b[i]*colsQ[i] : i in [1..#Eltseq(b)] ] : b in BC ];
A0, s := SubmatrixOfRank(A, #B div 2 : ColumnsOrRows := "Columns");

/* Modify invertible matrix A0 for numerical stability when taking quotients later */
CC := BaseRing(A0);
g := #Rows(A0); c := #Rows(Transpose(A0));
cols := [ Minimum([ i : i in [1..#Eltseq(row)] | Abs(row[i]) ge 10^20*CC`epscomp ]) : row in Rows(Transpose(A0)) ];
cols cat:= [ i : i in [1..g] | not i in cols ];
Perm := PermutationMatrix(CC, cols)^(-1);

/* Build period matrix */
V := VectorSpace(CC, g);
gensW := sub<V | [ c*Perm : c in Rows(Transpose(A0)) ] >;
W := sub< V | gensW >;
VmodW, q := quo< V | W >;
C := Transpose(Matrix([ Eltseq(q((V ! col)*Perm)) : col in colsC ]));

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
