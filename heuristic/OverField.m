/***
 *  Subrings of the geometric endomorphism ring over subfields
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic EndomorphismRepresentation(GeoEndoRep::SeqEnum, GalK::List) -> SeqEnum
{Extracts basis over subfield determines by a list of automorphisms.}

gensTan := [ gen[1] : gen in GeoEndoRep ];
gensHom := [ gen[2] : gen in GeoEndoRep ];
gensApp := [ gen[3] : gen in GeoEndoRep ];
L := BaseRing(gensTan[1]);

// Done if we are dealing with the rationals
if Type(L) eq FldRat then
    return GeoEndoRep;
end if;

gensH, Gphi := Explode(GalK);
K := FixedField(L, [ Gphi(genH) : genH in gensH ]);
// Done if we are dealing with the full extension
if #gensH eq 0 then
    return GeoEndoRep;
end if;

// The vector space representing the full endomorphism algebra
n := #gensTan;
Ker := VectorSpace(Rationals(), n);
// Successively filter by writing down the conditions for a matrix to be fixed
// under a given generator
for genH in gensH do
    sigma := Gphi(genH);
    Msigma := Matrix([ MatrixInBasis(ConjugateMatrix(sigma, genTan), gensTan) : genTan in gensTan ]);
    Msigma -:= IdentityMatrix(Rationals(), n);
    Ker meet:= Kernel(Msigma);
end for;

// Retrieve representations for a basis of of the endomorphism ring by taking a
// pure lattice / saturation
Lat := PureLattice(Lattice(Matrix(Basis(Ker))));
B := Basis(Lat);

// Constructing said basis
gens := [ ];
for b in B do
    genTan := &+[ b[i] * gensTan[i] : i in [1..n] ];
    // Coercion to subfield
    genTan := Matrix(K, genTan);
    genHom := &+[ b[i] * gensHom[i] : i in [1..n] ];
    genApp := &+[ b[i] * gensApp[i] : i in [1..n] ];
    Append(~gens, [* genTan, genHom, genApp *]);
end for;
return gens;

end intrinsic;


intrinsic EndomorphismRepresentation(GeoEndoRep::SeqEnum, K::Fld) -> SeqEnum
{Extracts basis over a general field.}

L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K);
return EndomorphismRepresentation(GeoEndoRep, GalK);

end intrinsic;


function CompareFields(K1, K2);
// Input:   Two subfields, fields, or polynomials.
// Output:  A comparison function: field with smaller degrees are smaller.

if Degree(K1) lt Degree(K2) then
    return -1;
elif Degree(K1) eq Degree(K2) then
    return 0;
else
    return 1;
end if;

end function;


intrinsic SubgroupGeneratorsUpToConjugacy(L::Fld, K::Fld) -> List
{Finds the subgroup generators up to conjugacy that correspond to the intersection of L and K, where L is Galois.}

if (Degree(K) eq 1) or (Type(L) eq FldRat) then
    Gp, Gf, Gphi := AutomorphismGroup(L);
    return [* Generators(Gp), Gphi *];
end if;

// Take smallest subfield of K that fits inside L and take corresponding fixed group
S := Subfields(K); S := [ tup[1] : tup in S ];
Sort(~S, CompareFields); Reverse(~S);
for M in S do
    test, f := IsSubfield(M, L);
    if test then
        KL := sub<L | f(M.1)>;
        break;
    end if;
end for;

Gp, Gf, Gphi := AutomorphismGroup(L);
if test then
    // Geometric case first:
    if Degree(K) eq Degree(L) then
        return [* [ ], Gphi *];
    end if;
    return [* Generators(FixedGroup(L, KL)), Gphi *];
else
    // Case where we need to go down all the way
    return [* Generators(Gp), Gphi *];
end if;

end intrinsic;
