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


intrinsic EndomorphismLattice(GeoEndoRep::SeqEnum, F::Fld) -> List
{Returns the lattice of endomorphisms by (conjugacy class of) subfield.}

F_seq := FieldDescription(F, BaseRing(F));
base := [* F_seq, F *];

L := BaseRing(GeoEndoRep[1][1][1]);
if (not IsRelativeExtension(L, F)) then
    entry, Shorthand, sthash := EndomorphismLatticeGeometricStep(GeoEndoRep, F);
    entries := [ entry ]; sthashes := [ sthash ]; realstrs := [ Sort(entry[2][2][3]) ];
    Gp := Sym(1); Hs := [ Gp ];
    return [* base, entries *], [* Gp, Hs, sthashes, realstrs *];
end if;

Gp, Gf, Gphi := AutomorphismGroup(L);
Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
Sort(~Hs, CompareGroups);

entry, Shorthand, sthash := EndomorphismLatticeGeometricStep(GeoEndoRep, F);
entries := [ entry ]; sthashes := [ sthash ]; realstrs := [ entry[2][2][3] ];
for H in Hs[2..#Hs] do
    gensH := Generators(H); GalK := [* gensH, Gphi *];
    entry, sthash := EndomorphismLatticeGeneralStep(GeoEndoRep, GalK, Shorthand, F);
    Append(~entries, entry); Append(~sthashes, sthash); Append(~realstrs, entry[2][2][3]);
end for;
return [* base, entries *], [* Gp, Hs, sthashes, realstrs *];

end intrinsic;


intrinsic EndomorphismLatticeGeometricStep(GeoEndoRep::SeqEnum, F::Fld) -> List
{Returns the geometric entry of the endomorphism lattice.}

entry := [* *];

L := BaseRing(GeoEndoRep[1][1][1]);
L_seq := FieldDescription(L, F);
L_desc := [* L_seq, L *];
Append(~entry, L_desc);

GalL := [* [ ], [ ] *];
EndoStruct := EndomorphismStructure(GeoEndoRep, GalL, F);
Append(~entry, EndoStruct);

//Append(~entry, ClassNumber(AbsoluteField(K)));

Shorthand := SatoTateShorthand(EndoStruct);
sthash := SatoTateHash(GeoEndoRep, GalL);
return entry, Shorthand, sthash;

end intrinsic;


intrinsic EndomorphismLatticeGeneralStep(GeoEndoRep::SeqEnum, GalK::List, Shorthand::MonStgElt, F::Fld) -> List
{Part of the above.}

entry := [* *];

L := BaseRing(GeoEndoRep[1][1][1]);
gensH, Gphi := Explode(GalK);
K := GeneralFixedField(L, [ Gphi(genH) : genH in gensH ]);
if HasBaseQQ(K) and not IsQQ(K) then
    K := Polredabs(K);
end if;
K_seq := FieldDescription(K, F);
K_desc := [* K_seq, K *];
Append(~entry, K_desc);

EndoStruct := EndomorphismStructure(GeoEndoRep, GalK, F : Shorthand := Shorthand);
Append(~entry, EndoStruct);

//Append(~entry, ClassNumber(AbsoluteField(K)));

sthash := SatoTateHash(GeoEndoRep, GalK);

return entry, sthash;

end intrinsic;
