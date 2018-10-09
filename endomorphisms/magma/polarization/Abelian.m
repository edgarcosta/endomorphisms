/***
 *  Abelian category functionality
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic KerModKer0(h::., P::., Q::.) -> .
{Returns the group of connected components of the kernel of the morphism h = (A, R) from P to Q.}
    A := h[1]; R := h[2];
    assert #Rows(R) eq #Rows(Transpose(Q)); assert #Rows(Transpose(R)) eq #Rows(Transpose(P));
    L := Lattice(Transpose(R));
    return quo< PureLattice(L) | L >;
end intrinsic;


intrinsic Ker0(h::., P::., Q::.) -> .
{Returns the connected component of the kernel of the morphism h = (A, R) from P to Q.}
    A := h[1]; R := h[2];
    assert #Rows(R) eq #Rows(Transpose(Q)); assert #Rows(Transpose(R)) eq #Rows(Transpose(P));
    RMS := RMatrixSpace(Integers(), #Rows(R), #Rows(Transpose(R)));
    /* Correct for right action taken by Magma */
    R0 := RMS ! R; B := Basis(Kernel(Transpose(R0)));
    if #B eq 0 then
        return 0, [0, 0];
    end if;
    B := [ Eltseq(b) : b in B ];
    colsP := Rows(Transpose(P)); colsK := [ &+[ b[i]*colsP[i] : i in [1..#b] ] : b in B ];
    K := Transpose(Matrix(colsK)); K, s := SubmatrixOfRank(K, #B div 2 : ColumnsOrRows := "Rows");
    ihom := Transpose(Matrix(Integers(), [ [ c : c in Eltseq(row) ] : row in B ]));
    itan := TangentRepresentation(ihom, K, P);
    return K, [* itan, ihom *];
end intrinsic;


intrinsic Coker(h::., P::., Q::.) -> .
{Returns the cokernel of the morphism h = (A, R) from P to Q.}
    A := h[1]; R := h[2];
    assert #Rows(R) eq #Rows(Transpose(Q)); assert #Rows(Transpose(R)) eq #Rows(Transpose(P));
    L := PureLattice(Lattice(Transpose(R))); Zd := StandardLattice(Degree(L));
    CL, projL := quo< Zd | L >; gensCL := [ gen : gen in Generators(CL) ];
    if #gensCL eq 0 then
        return 0, [0, 0];
    end if;
    B := [ Zd ! b : b in Basis(L) ]; BC := [ Zd ! g @@ projL : g in gensCL ];
    colsQ := Rows(Transpose(Q)); colsC := [ &+[ b[i]*colsQ[i] : i in [1..#Eltseq(b)] ] : b in BC ];
    /* Taking image and checking on which columns it gives a full subspace: */
    A0, s := SubmatrixOfRank(A, #B div 2 : ColumnsOrRows := "Columns");
    CC := BaseRing(A0);
    V := VectorSpace(CC, #Rows(A0));
    W := sub< V | [ Eltseq(c) : c in Rows(Transpose(A0)) ] >;
    VmodW, q := quo< V | W >;
    //C := Transpose(Matrix([ [ col[i] : i in [1..#Eltseq(col)] | not i in s ] : col in colsC ]));
    C := Transpose(Matrix([ Eltseq(q(col)) : col in colsC ]));
    /* Correction of homology representation: */
    M := Matrix([ Eltseq(b) : b in B cat BC ]);
    phom := VerticalJoin(ZeroMatrix(Integers(), #B, #BC), IdentityMatrix(Integers(), #BC));
    phom := Transpose(M^(-1)*phom);
    ptan := TangentRepresentation(phom, Q, C);
    return C, [* ptan, phom *];
end intrinsic;


intrinsic Img(h::., P::., Q::.) -> .
{Returns the image of the morphism h = (A, R) from P to Q.}

C, proj := Coker(h, P, Q);
if C eq 0 then
    A := IdentityMatrix(BaseRing(Q), #Rows(Q));
    R := IdentityMatrix(Integers(), #Rows(Transpose(Q)));
    return Q, [* A, R *];
end if;
return Ker0(proj, Q, C);

end intrinsic;
