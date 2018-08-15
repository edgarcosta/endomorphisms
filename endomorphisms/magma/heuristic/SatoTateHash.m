/***
 *  Hash for the Sato-Tate action
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic AutomorphismOrbits(G :: .) -> .
{To be added}
S:=ConjugacyClasses(G);
A:=AutomorphismGroup(G);
pi:=ClassMap(G);
T:=[i:i in [1..#S]];
for i in [1..#T] do
    for a in Generators(A) do
        T[pi(a(S[i][3]))] := T[i];
    end for;
end for;
X:=[[]:i in [1..#S]];
for i in [1..#S] do Append(~X[T[i]],<i,S[i][3]>); end for;
return [x:x in X|#x ne 0];
end intrinsic;


intrinsic AutomorphismOrbitsSubgroups(G :: .) -> .
{To be added}
S := [ rec`subgroup : rec in Subgroups(G) ];
A := AutomorphismGroup(G);
function pi(H, S)
for i in [1..#S] do
    if IsConjugate(G, H, S[i]) then
        return i;
    end if;
end for;
end function;
T := [i : i in [1..#S]];
for i in [1..#T] do
    for a in Generators(A) do
        T[pi(a(S[i]), S)] := T[i];
    end for;
end for;
X:=[[] : i in [1..#S]];
for i in [1..#S] do Append(~X[T[i]], <i,S[i]>); end for;
return [ x : x in X | #x ne 0 ];
end intrinsic;


intrinsic GaloisRepresentationOfElement(sigma :: ., GeoEndoRep :: .) -> .
{To be added}
rows := [ ];
As := [ galrep[1] : galrep in GeoEndoRep ]; L := BaseRing(As[1]);
for A in As do
    B := ConjugateMatrix(sigma, A);
    row := MatrixInBasis(B, As);
    Append(~rows, row);
end for;
return Matrix(rows);
end intrinsic;


intrinsic GaloisRepresentationOfConjugacyClasses(GeoEndoRep :: ., GalK :: .) -> .
{To be added}
As := [ galrep[1] : galrep in GeoEndoRep ]; L := BaseRing(As[1]);
gensGp, Gphi := Explode(GalK);
if #gensGp eq 0 then
    Gp := Sym(1);
    return Gp, [ < Id(Gp), IdentityMatrix(Rationals(), #As) > ];
end if;
Gp := sub< Domain(Gphi) | gensGp >;
CCs := ConjugacyClasses(Gp);
galreps := [ ];
for cc in CCs do
    M := GaloisRepresentationOfElement(Gphi(cc[3]), GeoEndoRep);
    Append(~galreps, < cc[3], M >);
end for;
return Gp, galreps;
end intrinsic;


intrinsic TracesOfConjugacyClasses(GeoEndoRep :: ., GalK :: .) -> .
{To be added}
Gp, galreps := GaloisRepresentationOfConjugacyClasses(GeoEndoRep, GalK);
return Gp, [ < rep[1], Trace(rep[2]) > : rep in galreps ];
end intrinsic;


intrinsic MapTraces(G :: ., trs :: .) -> .
{To be added}
D := SmallGroupDatabase();
id0 := IdentifyGroup(D, G);
o, n := Explode(id0);
G0 := SmallGroup(D, o, n);
test, iso := IsIsomorphic(G, G0);
trs0 := [ < iso(tr[1]), tr[2] > : tr in trs ];
return G0, trs0, id0;
end intrinsic;


intrinsic CanonizeTraces(G :: ., trs :: .) -> .
{To be added}
G0, trs0, id0 := MapTraces(G, trs);
f := ClassMap(G0);
trs0 := [ [ f(tr[1]), tr[2] ] : tr in trs0 ];
orbits := [ [ elt[1] : elt in orb ] : orb in AutomorphismOrbits(G0) ];
trs0tups := [ ];
for orb in orbits do
    trs0tup := [ ];
    for elt in orb do
        for tr in trs0 do
            if tr[1] eq elt then
                Append(~trs0tup, tr[2]);
            end if;
        end for;
    end for;
    Append(~trs0tups, Sort(trs0tup));
end for;
return [ [ c : c in id0 ], &cat(trs0tups) ];
end intrinsic;


intrinsic CanonizeRepresentation(GeoEndoRep :: .) -> .
{To be added}
As := [ galrep[1] : galrep in GeoEndoRep ]; L := BaseRing(As[1]);
Gp, Gf, Gphi := AutomorphismGroup(L);
D := SmallGroupDatabase();
id0 := IdentifyGroup(D, Gp);
o, n := Explode(id0);
G0 := SmallGroup(D, o, n);
test, iso := IsIsomorphic(G0, Gp);
galreps := [ GaloisRepresentationOfElement(Gphi(iso(g)), GeoEndoRep) : g in G0 ];
return < [ c : c in id0 ], [ Eltseq(M) : M in galreps ] >;
end intrinsic;


intrinsic SatoTateHash(GeoEndoRep :: ., GalK :: .) -> .
{To be added}
Gp, trs := TracesOfConjugacyClasses(GeoEndoRep, GalK);
return CanonizeTraces(Gp, trs);
end intrinsic;


function CompareHashes(hash1, hash2);
if hash1[1] lt hash2[1] then
    return -1;
elif hash1[1] gt hash2[1] then
    return 1;
end if;
if #hash1[2] lt #hash2[2] then
    return -1;
elif #hash1[2] gt #hash2[2] then
    return 1;
end if;
if hash1[2] lt hash2[2] then
    return -1;
elif hash1[2] gt hash2[2] then
    return 1;
end if;
return 0;
end function;


intrinsic CanonizeSatoTateHashes(sthashes :: .) -> .
{To be added}
Gp, Hs, hashes, realstrs := Explode(sthashes);
D := SmallGroupDatabase();
id0 := IdentifyGroup(D, Gp);
o, n := Explode(id0);
G0 := SmallGroup(D, o, n);
test, iso := IsIsomorphic(Gp, G0);
Hs0 := [ iso(H) : H in Hs ];
Kss0 := AutomorphismOrbitsSubgroups(G0);
hashes0 := [ ]; realstrs0 := [ ];
for Ks in Kss0 do
    hash0 := [ ]; realstr0 := [ ];
    for K in Ks do
        for i in [1..#Hs0] do
            H := Hs0[i];
            if IsConjugate(G0, H, K[2]) then
                Append(~hash0, hashes[i]);
                Append(~realstr0, realstrs[i]);
                break;
            end if;
        end for;
    end for;
    Append(~hashes0, Sort(hash0, CompareHashes));
    Append(~realstrs0, Sort(realstr0));
end for;
return [* [ c : c in id0 ], hashes0, realstrs0 *];
end intrinsic;
