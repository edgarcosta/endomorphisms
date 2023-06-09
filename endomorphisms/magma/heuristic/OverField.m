/***
 *  Subrings of the geometric endomorphism ring over subfields
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

forward SubgroupGeneratorsUpToConjugacy;


// FIXME: Sachi and Edgar: this should be called Homomorphism....
intrinsic EndomorphismRepresentation(GeoEndoRep::SeqEnum, GalK::List) -> SeqEnum
{Given a geometric representation GeoEndoRep and a list of automorphisms GalK of
their field of definition, returns a basis of the endomorphisms defined over
the subfield determined by GalK.}

As := [ gen[1] : gen in GeoEndoRep ];
Rs := [ gen[2] : gen in GeoEndoRep ];
L := BaseRing(As[1]); gensH, Gphi := Explode(GalK);

/* Case where no extension is needed to find the geometric endomorphism ring */
if Degree(L) eq 1 then
    return GeoEndoRep, CanonicalInclusionMap(L, L);
/* Case where we ask for the geometric endomorphism ring */
elif #gensH eq 0 then
    return GeoEndoRep, CanonicalInclusionMap(L, L);
else
    H := sub< Domain(Gphi) | gensH >;
    /* Case where we ask for the geometric endomorphism ring (again) */
    if #H eq 1 then
        return GeoEndoRep, CanonicalInclusionMap(L, L);
    end if;
end if;

/* The vector space representing the full endomorphism algebra */
n := #As;
Ker := VectorSpace(Rationals(), n);
/* Successively filter by writing down the conditions for a matrix to be fixed
 * under a given generator */
for genH in gensH do
    sigma := Gphi(genH);
    rows := [ ];
    for A in As do
        test, row := MatrixInBasis(ConjugateMatrix(sigma, A), As);
        assert test;
        Append(~rows, row);
    end for;
    Msigma := Matrix(rows);
    Msigma -:= IdentityMatrix(Rationals(), n);
    Ker meet:= Kernel(Msigma);
end for;

/* Retrieve representations for a basis of of the endomorphism ring by taking a
 * pure lattice / saturation */
Lat := PureLattice(Lattice(Matrix(Basis(Ker))));
B := Basis(Lat);

/* Constructing said basis */
gens := [ ];
K, hKL := FixedFieldExtra(L, [ Gphi(genH) : genH in gensH ]);
for b in B do
    A := &+[ b[i] * As[i] : i in [1..n] ];
    R := &+[ b[i] * Rs[i] : i in [1..n] ];
    /* Coercion to subfield */
    A := CoerceToSubfieldMatrix(A, L, K, hKL);
    Append(~gens, [* A, R *]);
end for;
return gens, hKL;

end intrinsic;


intrinsic EndomorphismRepresentation(GeoEndoRep::SeqEnum, K::Fld, h::Map) -> SeqEnum
{Given a geometric representation GeoEndoRep and a subfield K of their field of
definition, returns a basis of the endomorphisms defined over the subfield
determined by GalK.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(GeoEndoRep[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, h);
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


function SubgroupGeneratorsUpToConjugacy(L, K, hKL)
// Finds the subgroup generators up to conjugacy that correspond to the
// subfield K of L.

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
Gp, Gf, Gphi := AutomorphismGroupPari(L);
Helts := [ h : h in Gp | Gphi(h)(hKL(K.1)) eq hKL(K.1) ];
H := sub< Gp | Helts >;
return [* Generators(H), Gphi *];

end function;
