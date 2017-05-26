/**
 *  Lattices of subfield data
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


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


intrinsic EndomorphismLattice(GeoEndoRep::SeqEnum : Optimize := false) -> SeqEnum
{Returns the lattice of endomorphisms by (conjugacy class of) subfield.}

gensTan := [ gen[1] : gen in GeoEndoRep ];
gensHom := [ gen[2] : gen in GeoEndoRep ];
gensApp := [ gen[3] : gen in GeoEndoRep ];
L := BaseRing(gensTan[1]);
Gp, Gf, Gphi := AutomorphismGroup(L);

Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
Sort(~Hs, CompareGroups);

Lat := [ ];
// The code of this first (geometric) step is a copy of that below except for
// shorthand extraction. Of course it could be simpler, but there is no time
// loss as a result.
OverK := [* *];
gensH := [ ]; GalK := [* gensH, Gphi *];
K := L;

// TODO: Indicate class group and treat the relative case (scaffolding in place).
K_seq := [ Integers() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
K_desc := [* K_seq, K *];
EndoStruct := EndomorphismStructure(GeoEndoRep, GalK);
Append(~OverK, K_desc);
Append(~OverK, EndoStruct);
Append(~Lat, OverK);
Shorthand := SatoTateShorthand(EndoStruct);

for H in Hs[2..#Hs] do
    OverK := [* *];
    gensH := Generators(H); GalK := [* gensH, Gphi *];
    K := FixedField(L, [ Gphi(genH) : genH in gensH ]);

    K := ClearFieldDenominator(K);
    if (Type(K) eq FldNum and Optimize) then
        K := OptimizedRepresentation(K);
        K := ClearFieldDenominator(K);
    end if;

    // TODO: Indicate class group and treat the relative case (scaffolding in place).
    K_seq := [ Integers() ! c : c in Eltseq(MinimalPolynomial(K.1)) ];
    K_desc := [* K_seq, K *];
    EndoStruct := EndomorphismStructure(GeoEndoRep, GalK : Shorthand := Shorthand);
    Append(~OverK, K_desc);
    Append(~OverK, EndoStruct);
    Append(~Lat, OverK);
end for;
return Lat;

end intrinsic;
