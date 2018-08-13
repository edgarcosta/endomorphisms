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


intrinsic GaloisRepresentationOfConjugacyClasses(GeoEndoRep :: .) -> .
{To be added}
As := [ galrep[1] : galrep in GeoEndoRep ]; L := BaseRing(As[1]);
Gp, Gf, Gphi := AutomorphismGroup(L);
CCs := ConjugacyClasses(Gp);
galreps := [ ];
for cc in CCs do
    M := GaloisRepresentationOfElement(Gphi(cc[3]), GeoEndoRep);
    Append(~galreps, [* cc[3], M *]);
end for;
return Gp, galreps;
end intrinsic;


intrinsic TracesOfConjugacyClasses(GeoEndoRep :: .) -> .
{To be added}
Gp, galreps := GaloisRepresentationOfConjugacyClasses(GeoEndoRep);
return Gp, [ [* rep[1], Trace(rep[2]) *] : rep in galreps ];
end intrinsic;


intrinsic MapTraces(G :: ., trs :: .) -> .
{To be added}
D := SmallGroupDatabase();
id0 := IdentifyGroup(D, G);
o, n := Explode(id0);
G0 := SmallGroup(D, o, n);
test, iso := IsIsomorphic(G, G0);
trs0 := [ [* iso(tr[1]), tr[2] *] : tr in trs ];
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
return [* [ c : c in id0 ], [ Eltseq(M) : M in galreps ] *];
end intrinsic;


intrinsic SatoTateHash(GeoEndoRep :: .) -> .
{To be added}
G, trs := TracesOfConjugacyClasses(GeoEndoRep);
return CanonizeTraces(G, trs);
end intrinsic;
