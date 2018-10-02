/***
 *  Abelian category functionality
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

/* TODO: These have changed because of transposes. And add factorizations */


intrinsic KerModKer0(h::., P::., Q::.) -> .
{Returns the group of connected components of the kernel of the morphism h = (A, R) from P to Q.}
    A := h[1]; R := h[2];
    L := Lattice(R);
    return quo< PureLattice(L) | L >;
end intrinsic;


intrinsic Ker0(h::., P::., Q::.) -> .
{Returns the connected component of the kernel of the morphism h = (A, R) from P to Q.}
    A := h[1]; R := h[2];
    RMS := RMatrixSpace(Integers(), #Rows(R), #Rows(Transpose(R)));
    R0 := RMS ! R; B := Basis(Kernel(R0));
    if #B eq 0 then
        return 0, [0, 0];
    end if;
    B := [ Eltseq(b) : b in B ];
    rowsP := Rows(P); rowsK := [ &+[ b[i]*rowsP[i] : i in [1..#b] ] : b in B ];
    K := Matrix(rowsK); K, s := SubmatrixOfRank(K, #B div 2 : ColumnsOrRows := "Columns");
    ihom := Matrix(Integers(), [ [ c : c in Eltseq(row) ] : row in B ]);
    /* TODO: See if there is a better way to do this */
    itan := TangentRepresentation(ihom, K, P);
    return K, [* itan, ihom *];
end intrinsic;


intrinsic Coker(h::., P::., Q::.) -> .
{Returns the cokernel of the morphism h = (A, R) from P to Q.}
    A := h[1]; R := h[2];
    L := PureLattice(Lattice(R)); Zd := StandardLattice(Degree(L));
    CL, projL := quo< Zd | L >; gensCL := [ gen : gen in Generators(CL) ];
    if #gensCL eq 0 then
        return 0, [0, 0];
    end if;
    B := [ Zd ! b : b in Basis(L) ]; BC := [ Zd ! g @@ projL : g in gensCL ];
    rowsQ := Rows(Q); rowsC := [ &+[ b[i]*rowsQ[i] : i in [1..#Eltseq(b)] ] : b in BC ];
    /* Taking image and checking on which columns it gives a full subspace: */
    _, s := SubmatrixOfRank(A, #B div 2 : ColumnsOrRows := "Columns");
    C := Matrix([ [ row[i] : i in [1..#Eltseq(row)] | not i in s ] : row in rowsC]);
    /* Correction of homology representation: */
    M := Matrix([ Eltseq(b) : b in B cat BC ]);
    phom := VerticalJoin(ZeroMatrix(Integers(), #B, #BC), IdentityMatrix(Integers(), #BC));
    phom := M^(-1)*phom;
    /* TODO: See if there is a better way to do this */
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
