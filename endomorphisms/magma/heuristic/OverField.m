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


intrinsic EndomorphismRepresentationOverSubfield(GeoEndoRep::SeqEnum, GalK::List) -> SeqEnum
{Given a geometric representation GeoEndoRep and a list of automorphisms GalK of
their field of definition, returns a basis of the endomorphisms defined over
the subfield determined by GalK.}

As := [ gen[1] : gen in GeoEndoRep ];
Rs := [ gen[2] : gen in GeoEndoRep ];

/* Boundary cases */
L := BaseRing(As[1]);
gensH, Gphi := Explode(GalK);
/* Case where no extension is needed to find the geometric endomorphism ring */
if Degree(L) eq 1 then
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
n := #As;
Ker := VectorSpace(Rationals(), n);
/* Successively filter by writing down the conditions for a matrix to be fixed
 * under a given generator */
for genH in gensH do
    sigma := Gphi(genH);
    Msigma := Matrix([ MatrixInBasis(ConjugateMatrix(sigma, A), As) : A in As ]);
    Msigma -:= IdentityMatrix(Rationals(), n);
    Ker meet:= Kernel(Msigma);
end for;

/* Retrieve representations for a basis of of the endomorphism ring by taking a
 * pure lattice / saturation */
Lat := PureLattice(Lattice(Matrix(Basis(Ker))));
B := Basis(Lat);

/* Constructing said basis */
gens := [ ];
K, res := FixedFieldExtra(L, [ Gphi(genH) : genH in gensH ]);
K0, h := ImproveFieldExtra(K);
for b in B do
    A := &+[ b[i] * As[i] : i in [1..n] ];
    R := &+[ b[i] * Rs[i] : i in [1..n] ];
    /* Coercion to subfield */
    A := CoerceToSubfieldMatrix(A, L, K, res);
    A := Matrix([ [ h(c) : c in Eltseq(row) ] : row in Rows(A) ]);
    Append(~gens, [* A, R *]);
end for;
return gens;

end intrinsic;


intrinsic EndomorphismRepresentationOverSubfield(GeoEndoRep::SeqEnum, K::Fld, h::Map) -> SeqEnum
{Given a geometric representation GeoEndoRep and a subfield K of their field of
definition, returns a basis of the endomorphisms defined over the subfield
determined by GalK.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, h);
return EndomorphismRepresentationOverSubfield(GeoEndoRep, GalK);

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


intrinsic SubgroupGeneratorsUpToConjugacy(L::Fld, K::Fld, h::Map) -> List
{Finds the subgroup generators up to conjugacy that correspond to the
intersection of the extensions L and K of their common base field.}

/* Case where L and K coincide */
if L eq K then
    return [* [ ], [ ] *];
end if;

/* Case where either K or L is small */
if (Degree(K) eq 1) or (Degree(L) eq 1) then
    Gp, Gf, Gphi := AutomorphismGroupPari(L);
    return [* Generators(Gp), Gphi *];
end if;

/* General case: take group corresponding to largest subfield of L that fits
 * inside K */
/* TODO: Deal with case where K is not a subfield of L, but not as cruddily as
 * Magma */
Gp, Gf, Gphi := AutomorphismGroupPari(L);
Helts := [ h : h in Gp | Gphi(h)(h(K.1)) eq h(K.1) ];
H := sub< Gp | Helts >;
return [* Generators(H), Gphi *];

end intrinsic;
