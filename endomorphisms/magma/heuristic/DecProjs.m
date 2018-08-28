/***
 *  Elliptic curves from decompositions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic IdempotentsFromLattice(Lat::List) -> .
{Given a lattice Lat, returns idempotents over a small field extension.}
/* TODO: Change field to be even smaller */

entries := Lat[2];
entries := Reverse(entries);
GeoEndoStruct := entries[#entries][2];
L := entries[#entries][1][2];
num_idemsgeo := NumberOfIdempotentsFromStructure(GeoEndoStruct);
if num_idemsgeo eq 1 then
    /* Stop if there are no geometric idempotents and give some dummy output */
    K := entries[1][1][2];
    RestrictInfinitePlace(L, K);
    return [ ], K;
end if;

/* Now go up until we meet the same number of idempotents */
n := #entries; i := 1;
while i le n do
    EndoStruct := entries[i][2];
    num_idems := NumberOfIdempotentsFromStructure(EndoStruct);
    if num_idems eq num_idemsgeo then
        /* Get elements, their base field, and inherit embedding */
        idems := IdempotentsFromStructure(EndoStruct);
        K := BaseRing(idems[1][1]);
        RestrictInfinitePlace(L, K);
        return idems, K;
    end if;
    i +:= 1;
end while;

end intrinsic;


intrinsic NumberOfIdempotentsFromStructure(EndoStruct::List) -> RngIntElt
{Returns the number of idempotents for the endomorphism structure EndoStruct.}

EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
factors_QQ := EndoDesc[1];
num_idems := 0;
for factor_QQ in factors_QQ do
    albert, _, dim_sqrt, disc := Explode(factor_QQ);
    // TODO: the usual nastiness with powers of a quaternion algebra is again
    // not covered, so watch out in that case.
    if (albert in ["I", "IV"]) and (disc eq 1) then
        num_idems +:= dim_sqrt;
    else
        num_idems +:= 1;
    end if;
end for;
return num_idems;

end intrinsic;


intrinsic IdempotentsFromStructure(EndoStruct::List) -> List
{Returns idempotents for the endomorphism structure EndoStruct.}

g := #Rows(EndoStruct[1][1][1]);
EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
C, GensC := Explode(EndoAlg);
Ds := DirectSumDecomposition(C);
idemsC := &cat[ IdempotentsFromFactor(D, C, g) : D in Ds ];
idemsrep := MatricesFromIdempotents(idemsC, EndoStruct);
return idemsrep;

end intrinsic;


intrinsic IdempotentsFromFactor(D::., C::., g::RngIntElt) -> .
{Returns idempotents for the factor D in the direct sum decomposition of the
endomorphism algebra C. The genus of the curve involved in the background needs
to be passed as g.}

if g le 3 then
    return IdempotentsFromFactorG3(D, C);
else
    error "Higher genus not implemented yet";
end if;

end intrinsic;


intrinsic IdempotentsFromFactorG3(D::., C::.) -> .
{Returns idempotents for the factor D in the direct sum decomposition of the
endomorphism algebra C. The genus of the curve involved in the background needs
to be passed as g. We assume the genus to be at most 3.}

E1, f1 := AlgebraOverCenter(D);
//F := ClearFieldDenominator(BaseRing(E1));
//if Type(F) eq FldNum then
//    F := ClearFieldDenominator(F);
//    F := Polredabs(F);
//end if;
//E2, f2 := ChangeRing(E1, F);
E2 := E1;
if not IsCommutative(E2) then
    test_dim, d := IsSquare(Dimension(E2));
    if d eq 2 then
        test_quat, Q, f3 := IsQuaternionAlgebra(E2);
        if test_quat then
            test_mat, M, f4 := IsMatrixRing(Q : Isomorphism := true);
            //f := f1 * f2 * f3 * f4;
            f := f1 * f3 * f4;
            invf := Inverse(f);
            return [ C ! invf(M ! [1,0,0,0]), C ! invf(M ! [0,0,0,1]) ];
        end if;
    elif d eq 3 then
        //idems_E2 := IdempotentsInMatrixAlgebra(E2);
        //invf1 := Inverse(f1);
        //return [ C ! invf1(idem_E2) : idem_E2 in idems_E2 ];
        error "All cases in IdempotentsFromFactorG3 fell through";
    else
        error "All cases in IdempotentsFromFactorG3 fell through";
    end if;
end if;
return [ C ! D ! 1 ];

end intrinsic;


intrinsic MatricesFromIdempotents(idems::SeqEnum, EndoStruct::List) -> SeqEnum
{Returns the matrix representations corresponding to the idempotents in idems,
using the endomorphism structure EndoStruct.}

EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
GensTan := [ gen[1] : gen in EndoRep ];
GensHom := [ gen[2] : gen in EndoRep ];
GensApp := [ gen[3] : gen in EndoRep ];
C, GensC := Explode(EndoAlg);

idems := [ [ Rationals() ! c : c in Eltseq(idem) ] : idem in idems ];
GensC := [ [ Rationals() ! c : c in Eltseq(gen) ] : gen in GensC ];
idems := [ Eltseq(MatrixInBasis(idem, GensC)) : idem in idems ];

idemsRep := [ ];
for idem in idems do
    idemTan := &+[ idem[i] * GensTan[i] : i in [1..#GensTan] ];
    idemHom := &+[ idem[i] * GensHom[i] : i in [1..#GensHom] ];
    idemApp := &+[ idem[i] * GensApp[i] : i in [1..#GensApp] ];
    Append(~idemsRep, [* idemTan, idemHom, idemApp *]);
end for;
return idemsRep;

end intrinsic;


intrinsic IdempotentDenominator(idem::.) -> RngIntElt
{Returns the degree of the morphism corresponding to the idempotent idem.}
/* TODO: Use polarizations to get at the degree instead */

idemTan, idemHom, idemApp := Explode(idem);
return LCM([ Denominator(c) : c in Eltseq(idemHom) ]);

end intrinsic;


intrinsic ProjectionFromIdempotent(P::., idem::List) -> List
{Given an idempotent idem for the period matrix P, returns a corresponding
lattice and an analytic representation of the projection to it.}

/* g is the dimension of the new abelian variety */
A, R, ACC := Explode(idem); g := Rank(R) div 2;
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
Q := CombineMatrix(QSplit, CC);
/* We now have QLarge = Q U, so B P = Q U R. We return B and U R. */

S := U * R;
proj := [* B, S, BCC *];
if Maximum([ Abs(c) : c in Eltseq(BCC*P - Q*ChangeRing(S, CC)) ]) gt CC`epscomp then
    error "Error in determining projection";
end if;
return Q, proj;

end intrinsic;
