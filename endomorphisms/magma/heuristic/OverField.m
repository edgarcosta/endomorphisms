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


intrinsic EndomorphismRepresentation(GeoEndoRep::SeqEnum, GalK::List, F::Fld) -> SeqEnum
{Given a geometric representation GeoEndoRep, a list of automorphisms GalK of
their field of definition, and a field F, returns a basis of the endomorphisms
defined over the subfield determined by GalK.}

gensTan := [ gen[1] : gen in GeoEndoRep ];
gensHom := [ gen[2] : gen in GeoEndoRep ];
gensApp := [ gen[3] : gen in GeoEndoRep ];

/* Boundary cases */
L := BaseRing(gensTan[1]);
gensH, Gphi := Explode(GalK);
/* Case where no extension is needed to find the geometric endomorphism ring */
if not IsRelativeExtension(L, F) then
    return GeoEndoRep;
/* Case where we ask for the geometric endomorphism ring */
elif #gensH eq 0 then
    return GeoEndoRep;
else
    H := sub< Domain(Gphi) | gensH >;
    /* Case where we ask for the geometric endomorphism ring (again) */
    if #H eq 1 then
        return GeoEndoRep;
    end if;
end if;

/* The vector space representing the full endomorphism algebra */
n := #gensTan;
Ker := VectorSpace(Rationals(), n);
/* Successively filter by writing down the conditions for a matrix to be fixed
 * under a given generator */
for genH in gensH do
    sigma := Gphi(genH);
    Msigma := Matrix([ MatrixInBasis(ConjugateMatrix(sigma, genTan), gensTan) : genTan in gensTan ]);
    Msigma -:= IdentityMatrix(Rationals(), n);
    Ker meet:= Kernel(Msigma);
end for;

/* Retrieve representations for a basis of of the endomorphism ring by taking a
 * pure lattice / saturation */
Lat := PureLattice(Lattice(Matrix(Basis(Ker))));
B := Basis(Lat);

/* Constructing said basis */
gens := [ ];
K := GeneralFixedField(L, [ Gphi(genH) : genH in gensH ]);
for b in B do
    genTan := &+[ b[i] * gensTan[i] : i in [1..n] ];
    /* Coercion to subfield */
    genTan := Matrix(K, genTan);
    genHom := &+[ b[i] * gensHom[i] : i in [1..n] ];
    genApp := &+[ b[i] * gensApp[i] : i in [1..n] ];
    Append(~gens, [* genTan, genHom, genApp *]);
end for;
return gens;

end intrinsic;


intrinsic EndomorphismRepresentation(GeoEndoRep::SeqEnum, K::Fld, F::Fld) -> SeqEnum
{Given a geometric representation GeoEndoRep, a subfield K of their field of
definition, and a field F, returns a basis of the endomorphisms defined over
the subfield determined by GalK.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, F);
return EndomorphismRepresentation(GeoEndoRep, GalK, F);

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


function CompareGroups(G1, G2);
// Input:   Two subgroups or groups.
// Output:  A comparison function: groups with smaller cardinality are smaller.

if #G1 lt #G2 then
    return -1;
elif #G1 eq #G2 then
    return 0;
else
    return 1;
end if;

end function;


intrinsic SubgroupGeneratorsUpToConjugacy(L::Fld, K::Fld, F::Fld) -> List
{Finds the subgroup generators up to conjugacy that correspond to the
intersection of the extensions L and K of F. It is assumed that L is Galois
and, for now, that F is the field of rationals.}

/* Case where L and K coincide */
if L eq K then
    return [* [ ], [ ] *];
end if;

/* Case where either K or L is small */
if (Degree(K) eq 1) or (Degree(L) eq 1) then
    Gp, Gf, Gphi := AutomorphismGroup(L);
    return [* Generators(Gp), Gphi *];
end if;

/* General case: take group corresponding to largest subfield of L that fits
 * inside K */
/* TODO: This is not very elegant, but the reason for this is that FixedGroup
 * fails for relative extensions */
Gp, Gf, Gphi := AutomorphismGroup(L);
Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
Sort(~Hs, CompareGroups);
for H in Hs do
    M := FixedField(L, [ Gphi(h) : h in H ]);
    test, f := IsSubfield(M, K);
    if test then
        return [* Generators(H), Gphi *];
    end if;
end for;

end intrinsic;
