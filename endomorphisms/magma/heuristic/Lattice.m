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


intrinsic EndomorphismLattice(GeoEndoRep::SeqEnum) -> List
{Returns the lattice of endomorphisms by (conjugacy class of) subfield.}

L := BaseRing(GeoEndoRep[1][1]);
F := BaseRing(L); F_seq := FieldDescription(F);
base := [* F_seq, F *];

if Degree(L) eq 1 then
    entry, stpart, Shorthand := EndomorphismLatticeGeometricStep(GeoEndoRep);
    entries := [ entry ]; stparts := [ stpart ]; realstrs := [ entry[2][2][3] ];
    Gp := Sym(1); Hs := [ Gp ];
    return [* base, entries *], CanonizeSatoTateHash([* Gp, Hs, stparts, realstrs *]);
end if;

Gp, Gf, Gphi := AutomorphismGroup(L);
Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
Sort(~Hs, CompareGroups);

entry, stpart, Shorthand := EndomorphismLatticeGeometricStep(GeoEndoRep);
entries := [ entry ]; stparts := [ stpart ]; realstrs := [ entry[2][2][3] ];
for H in Hs[2..#Hs] do
    gensH := Generators(H); GalK := [* gensH, Gphi *];
    entry, stpart := EndomorphismLatticeGeneralStep(GeoEndoRep, GalK, Shorthand);
    Append(~entries, entry); Append(~stparts, stpart); Append(~realstrs, entry[2][2][3]);
end for;
return [* base, entries *], CanonizeSatoTateHash([* Gp, Hs, stparts, realstrs *]);

end intrinsic;


intrinsic EndomorphismLatticeGeometricStep(GeoEndoRep::SeqEnum) -> List
{Returns the geometric entry of the endomorphism lattice.}

entry := [* *];

L := BaseRing(GeoEndoRep[1][1][1]);
L_seq := FieldDescription(L); L_desc := [* L_seq, L *];
Append(~entry, L_desc);

GalL := [* [ ], [ ] *];
EndoStruct := EndomorphismStructure(GeoEndoRep, GalL);
Append(~entry, EndoStruct);

Shorthand := SatoTateShorthand(EndoStruct);
stpart := SatoTateHashPart(GeoEndoRep, GalL);
return entry, stpart, Shorthand;

end intrinsic;


intrinsic EndomorphismLatticeGeneralStep(GeoEndoRep::SeqEnum, GalK::List, Shorthand::MonStgElt) -> List
{Part of the above.}

entry := [* *];

L := BaseRing(GeoEndoRep[1][1]);
gensH, Gphi := Explode(GalK);
vprint EndoFind: "";
vprint EndoFind: "Galois group:", gensH;
K := FixedFieldExtra(L, [ Gphi(genH) : genH in gensH ]);
vprint EndoFind: "Base field:", K;
K := ImproveField(K); K_seq := FieldDescription(K);
vprint EndoFind: "Improved base field:", K;
K_desc := [* K_seq, K *];
Append(~entry, K_desc);

EndoStruct := EndomorphismStructure(GeoEndoRep, GalK : Shorthand := Shorthand);
Append(~entry, EndoStruct);

stpart := SatoTateHashPart(GeoEndoRep, GalK);

return entry, stpart;

end intrinsic;
