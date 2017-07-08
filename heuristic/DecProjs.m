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
{Finds idempotents over the smallest possible field by running over Lat.}

entries := Lat[2];
entries := Reverse(entries);
GeoEndoStruct := entries[#entries][2];
L := entries[#entries][1][2];
num_idemsgeo := NumberOfIdempotentsFromStructure(GeoEndoStruct);
if num_idemsgeo eq 1 then
    K := entries[1][1][2];
    RestrictInfinitePlace(L, K);
    return [ ], K;
end if;

n := #entries; i := 1;
while i le n do
    K := entries[i][1][2];
    EndoStruct := entries[i][2];
    num_idems := NumberOfIdempotentsFromStructure(EndoStruct);
    if num_idems eq num_idemsgeo then
        idems := IdempotentsFromStructure(EndoStruct);
        RestrictInfinitePlace(L, K);
        return idems, K;
    end if;
    i +:= 1;
end while;

end intrinsic;


intrinsic NumberOfIdempotentsFromStructure(EndoStruct::List) -> RngIntElt
{Finds number of idempotents.}

EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
factors_QQ := EndoDesc[1];
num_idems := 0;
for factor_QQ in factors_QQ do
    albert, _, dim_sqrt, disc := Explode(factor_QQ);
    // TODO: the usual nastyness with powers of a quaternion algebra is again
    // not covered.
    if (albert in ["II", "IV"]) and (disc eq 1) then
        num_idems +:= dim_sqrt;
    else
        num_idems +:= 1;
    end if;
end for;
return num_idems;

end intrinsic;


intrinsic IdempotentsFromStructure(EndoStruct::List) -> List
{Finds idempotents.}

g := #Rows(EndoStruct[1][1][1]);
EndoRep, EndoAlg, EndoDesc := Explode(EndoStruct);
C, GensC := Explode(EndoAlg);
Ds := DirectSumDecomposition(C);
idemsC := &cat[ IdempotentsFromFactor(D, C, g) : D in Ds ];
idemsrep := MatricesFromIdempotents(idemsC, EndoStruct);
return idemsrep;

end intrinsic;


intrinsic IdempotentsFromFactor(D::., C::., g::RngIntElt) -> .
{Idempotents originating from the factor D.}

if g le 3 then
    return IdempotentsFromFactorG3(D, C);
else
    error "Higher genus not implemented yet";
end if;

end intrinsic;


intrinsic IdempotentsFromFactorG3(D::., C::.) -> .
{Idempotents originating from the factor D.}

// TODO: This 
E1, f1 := AlgebraOverCenter(D);
//F := ClearFieldDenominator(BaseRing(E1));
//if Type(F) eq FldNum then
//    F := OptimizedRepresentation(F);
//    F := ClearFieldDenominator(F);
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
{Recovers matrices corresponding to idems.}

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


intrinsic ProjectionFromIdempotent(P::., idem::List) -> List
{From an idempotent, extracts corresponding lattice and projection.}

idemTan, idemHom, idemApp := Explode(idem);
// Extract the complex field
CC := BaseRing(P); RR := RealField(CC);

// Create analytic idempotent and project:
PEllHuge := P * Transpose(idemApp);

// Compute rank of projection
gQuotient := NumericalRank(PEllHuge : Epsilon := RR`epsinv);

PEllHugeSplit := SplitMatrix(PEllHuge);

PEllBigSplit, row_numbers := SubmatrixOfRank(PEllHugeSplit, 2*gQuotient : ColumnsOrRows := "Rows"); // extract 2g independent rows
PEllBigSplit, M := SaturateLattice(PEllHugeSplit, PEllBigSplit); // ensure that these rows generate the full lattice

PEllBig := CombineMatrix(PEllBigSplit, CC); // go back to the complex representation

LatticeMatrix, s0 := SubmatrixOfRank(PEllBig, gQuotient : ColumnsOrRows := "Columns"); // extract g columns (i.e. decide which projection to use)

/* Do not see the need for the next step:
PreliminaryLatticeMatrix := Transpose(PreliminaryLatticeMatrix); // necessary before calling SaturateLattice
PEllBig := Transpose(PEllBig);
LatticeMatrix := Transpose(SaturateLattice(PEllBig, PreliminaryLatticeMatrix));
*/

rowsTan := Rows(idemTan); projTan := Matrix([ rowsTan[i] : i in s0 ]);
rowsHom := Rows(idemHom); projHom := Matrix([ rowsHom[i] : i in s0 ]);
rowsApp := Rows(idemApp); projApp := Matrix([ rowsApp[i] : i in s0 ]);
d := IdempotentDenominator(idem);
proj := [* d*projTan, d*projHom, d*projApp *];

return LatticeMatrix, proj;
end intrinsic;


intrinsic IdempotentDenominator(idem::.) -> RngIntElt
{Degree of morphism from given idempotent.}

idemTan, idemHom, idemApp := Explode(idem);
return LCM([ Denominator(c) : c in Eltseq(idemHom) ]);

end intrinsic;
