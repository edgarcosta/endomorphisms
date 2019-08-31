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

import "SatoTate.m": SatoTateShorthand;
//forward EndomorphismLatticeGeometricStep;
//forward EndomorphismLatticeGeneralStep;
forward EndomorphismLatticeStep;


intrinsic EndomorphismLattice(GeoEndoRep::SeqEnum) -> List
{Returns the lattice of endomorphisms by (conjugacy class of) subfield.}

vprint EndoFind : "";
vprint EndoFind : "Determining endomorphism lattice...";
L := BaseRing(GeoEndoRep[1][1]); Gp, Gf, Gphi := AutomorphismGroupPari(L);
Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];

entries := [ ];
for H in Hs do
    gensH := Generators(H); GalK := [* gensH, Gphi *];
    entry := EndomorphismLatticeStep(GeoEndoRep, GalK);
    Append(~entries, entry);
end for;

GalL := IdentifyGroup(AutomorphismGroupPari(L));
EndRR := entries[1][3];

vprint EndoFind : "done determining endomorphism lattice.";
return < EndRR, GalL, entries >;

end intrinsic;


function EndomorphismLatticeStep(GeoEndoRep, GalK);
// Returns the entry of the endomorphism lattice over the field corresponding
// to GalK.

entry := < >;
gensH, Gphi := Explode(GalK); H := sub< Domain(Gphi) | gensH >;
Append(~entry, IdentifyGroup(H));

vprint EndoFind, 2 : "";
vprint EndoFind, 2 : "Generators of Galois group in lattice:", gensH;
L := BaseRing(GeoEndoRep[1][1]);
K := FixedFieldExtra(L, [ Gphi(genH) : genH in gensH ]);
K_seq := FieldDescriptionExtra(K);
vprint EndoFind, 2 : "Corresponding field:", K;
Append(~entry, < c : c in K_seq >);

EndoStruct := EndomorphismData(GeoEndoRep, GalK); EndoDesc := EndoStruct[3];
for x in EndoDesc do Append(~entry, x); end for;
return EndoDesc;

end function;


//function CompareGroups(G1, G2);
//// Input:   Two subgroups or groups.
//// Output:  A comparison function: groups with smaller cardinality are smaller.
//
//if #G1 lt #G2 then
//    return -1;
//elif #G1 eq #G2 then
//    return 0;
//else
//    return 1;
//end if;
//
//end function;


//intrinsic EndomorphismLattice(GeoEndoRep::SeqEnum) -> List
//{Returns the lattice of endomorphisms by (conjugacy class of) subfield.}
//
//vprint EndoFind : "";
//vprint EndoFind : "Determining endomorphism lattice...";
//L := BaseRing(GeoEndoRep[1][1]);
//F := L`base; F_seq := FieldDescriptionExtra(F);
//base := [* F_seq, F *];
//
//if Degree(L) eq 1 then
//    entry, stpart, Shorthand := EndomorphismLatticeGeometricStep(GeoEndoRep);
//    entries := [ entry ]; stparts := [ stpart ]; realstrs := [ entry[2][2][3] ];
//    Gp := Sym(1); Hs := [ Gp ];
//    vprint EndoFind : "done determining endomorphism lattice.";
//    return [* base, entries *], CanonizeSatoTateHash([* Gp, Hs, stparts, realstrs *]);
//end if;
//
//Gp, Gf, Gphi := AutomorphismGroupPari(L);
//Hs := Subgroups(Gp); Hs := [ H`subgroup : H in Hs ];
//Sort(~Hs, CompareGroups);
//
//entry, stpart, Shorthand := EndomorphismLatticeGeometricStep(GeoEndoRep);
//entries := [ entry ]; stparts := [ stpart ]; realstrs := [ entry[2][2][3] ];
//for H in Hs[2..#Hs] do
//    gensH := Generators(H); GalK := [* gensH, Gphi *];
//    entry, stpart := EndomorphismLatticeGeneralStep(GeoEndoRep, GalK, Shorthand);
//    Append(~entries, entry); Append(~stparts, stpart); Append(~realstrs, entry[2][2][3]);
//end for;
//vprint EndoFind : "done determining endomorphism lattice.";
//return [* base, entries *];
//
//end intrinsic;


//function EndomorphismLatticeGeometricStep(GeoEndoRep)
//// Returns the geometric entry of the endomorphism lattice.
//
//entry := [* *];
//
//L := BaseRing(GeoEndoRep[1][1][1]);
//L_seq := FieldDescriptionExtra(L); L_desc := [* L_seq, L *];
//Append(~entry, L_desc);
//
//GalL := [* [ ], [ ] *];
//vprint EndoFind, 2 : "";
//vprint EndoFind, 2 : "Generators of Galois group in lattice:", [ ];
//vprint EndoFind, 2 : "Corresponding field:", L;
//EndoStruct := EndomorphismDataWithSatoTate(GeoEndoRep, GalL);
//Append(~entry, EndoStruct);
//
//Shorthand := SatoTateShorthand(EndoStruct);
//return entry, Shorthand;
//
//end function;


//function EndomorphismLatticeGeneralStep(GeoEndoRep, GalK, Shorthand)
//// Returns the entry of the endomorphism lattice over the field corresponding
//// to GalK.
//
//entry := [* *];
//
//L := BaseRing(GeoEndoRep[1][1]);
//gensH, Gphi := Explode(GalK);
//vprint EndoFind, 2 : "";
//vprint EndoFind, 2 : "Generators of Galois group in lattice:", gensH;
//K := FixedFieldExtra(L, [ Gphi(genH) : genH in gensH ]);
//K_seq := FieldDescriptionExtra(K);
//vprint EndoFind, 2 : "Corresponding field:", K;
//K_desc := [* K_seq, K *];
//Append(~entry, K_desc);
//
//EndoStruct := EndomorphismDataWithSatoTate(GeoEndoRep, GalK : Shorthand := Shorthand);
//Append(~entry, EndoStruct);
//
//stpart := SatoTateHashPart(GeoEndoRep, GalK);
//
//return entry, stpart;
//
//end function;
