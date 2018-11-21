/***
 *  Maximal isotropic subgroups and resulting lattices
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

forward IsogenousPPLatticesG2;


intrinsic InducedPolarization(E::., R::. : ProjToIdem := true) -> .
{Given a matrix E corresponding to a polarization, returns the pushforward (default) or pullback of E along R. The pullback is the pairing directly induced via the map represented by R, and the pushforward is its dual.}

if ProjToIdem then
    Q := R*E*Transpose(R);
else
    Q := Transpose(R)*E*R;
    Q := Q^(-1);
    Q := -Denominator(Q)*Q;
end if;
return Q;

end intrinsic;


intrinsic FrobeniusFormAlternatingAlt(E::.) -> .
{Returns a different standard form E0 and a matrix T with T E Transpose(T) = E0.}

/* Some coercions are involved to make everything defined over QQ */
E1, T1 := FrobeniusFormAlternating(Matrix(ChangeRing(E, Integers())));
g := #Rows(E) div 2; S := Sym(2*g);
sigma := S ! (&cat[ [ i, g + i ] : i in [1..g] ]);
P := PermutationMatrix(Integers(), sigma);
E2 := P*E1*Transpose(P);
T2 := P*T1;
return ChangeRing(E2, Rationals()), ChangeRing(T2, Rationals());

end intrinsic;


function SymplecticSubmodulesPrime(p, d)
// This is really stupid: right cosets are better. However, do not see how to
// do that now and can get by without.

assert d mod 2 eq 0;
FF := FiniteField(p);
V := VectorSpace(FF, d);
B0 := [ V.i : i in [1..(d div 2)] ];
G := SymplecticGroup(d, FF);
Ws := [ ];
for g in G do
    W := sub< V | [ b*g : b in B0 ] >;
    if not W in Ws then
        Append(~Ws, W);
    end if;
end for;
return Ws;

end function;


function SymplecticSubmodulesPrimePower(pf, d)
// Based on a suggestion of John Voight

assert d mod 2 eq 0;
test, p, f := IsPrimePower(pf);
if f eq 1 then
    return SymplecticSubmodulesPrime(p, d);
end if;

R := quo< Integers() | pf >;
M := RSpace(R, d);

bases := CartesianPower(M, d);
submods := [ ];
for tup in bases do
    basis := [ e : e in tup ];
    /* Check that elements generate */
    if sub< M | basis > eq M then
        factors := [ ];
        /* In the end we have to consider d/2 pairs */
        for i in [1..(d div 2)] do
            factor := [ ];
            /* The pairs per (i, i + 1) */
            for e in [0..(f div 2)] do
                Append(~factor, [ (p^e)*basis[2*i - 1], p^(f - e)*basis[2*i] ]);
            end for;
            Append(~factors, factor);
        end for;
        /* Now take cartesian product and keep new spaces */
        CP := CartesianProduct(factors);
        for tup in CP do
            newbasis := &cat[ pair : pair in tup ];
            N := sub< M | newbasis >;
            if not N in submods then
                Append(~submods, N);
            end if;
        end for;
    end if;
end for;
return submods;

end function;


intrinsic SymplecticSubmodules(n::RngIntElt, d::RngIntElt : ProjToPP := true) -> .
{All symplectic submodules of index n in rank 2*d, or alternatively the maximal
symplectic submodules of (ZZ / n ZZ)^(2*d) with the canonical form.}

assert d mod 2 eq 0;
Fac := Factorization(n);
pfs := [ tup[1]^tup[2] : tup in Fac ];
L0 := Lattice(IdentityMatrix(Rationals(), d));
Lats := [ L0 ];

for pf in pfs do
    submods := SymplecticSubmodulesPrimePower(pf, d);
    Latsnew := [ ];
    for Lat in Lats do
        for submod in submods do
            if ProjToPP then
                B := ChangeRing(Matrix(Basis(Lat)), Rationals());
                M := pf * B;
                Mnew := Matrix(Rationals(), [ [ Integers() ! c : c in Eltseq(gen) ] : gen in Generators(submod) ]) * B;
            else
                M := ChangeRing(Matrix(Basis(Lat)), Rationals());;
                Mnew := (1/pf) * Matrix(Rationals(), [ [ Rationals() ! Integers() ! c : c in Eltseq(gen) ] : gen in Generators(submod) ]);
            end if;
            Append(~Latsnew, Lattice(VerticalJoin(M, Mnew)));
        end for;
    end for;
    Lats := Latsnew;
end for;
return Lats;

end intrinsic;


intrinsic IsogenousPPLattices(E::. : ProjToPP := true) -> .
{Given an alternating form E, finds the sublattices to ZZ^2d of smallest possible index on which E induces a principal polarization. These are returned in matrix form, that is, as a span of a basis in the rows. This basis is symplectic in the usual sense.}

g := #Rows(E) div 2;
if g eq 1 then
    return [ E ];
elif g eq 2 then
    return IsogenousPPLatticesG2(E : ProjToPP := ProjToPP);
else
    error "No functionality for principally polarized lattices in genus strictly larger than 2 yet";
end if;

end intrinsic;


/* TODO: Generalize (note that description below is not quite correct) */
function IsogenousPPLatticesG2(E : ProjToPP := true)
// Given an alternating form E, finds the sublattices to ZZ^2d of smallest
// possible index on which E induces a principal polarization. These are
// returned in matrix form, that is, as a span of a basis in the rows. This
// basis is symplectic in the usual sense.
/* In general, we would isolate the blocks with given d and deal with those one at a time */

E0, T0 := FrobeniusFormAlternatingAlt(E);
d := GCD([ Integers() ! c : c in Eltseq(E0) ]);
E0 := (1/d)*E0;
n := Integers() ! Abs(E0[3,4]);

Ts := [ ];
for Lat in SymplecticSubmodules(n, 2 : ProjToPP := ProjToPP) do
    Ur := ChangeRing(Matrix(Basis(Lat)), Rationals());
    if ProjToPP then
        U := DiagonalJoin(Ur, IdentityMatrix(Rationals(), 2));
    else
        U := DiagonalJoin(IdentityMatrix(Rationals(), 2), Ur);
    end if;
    Append(~Ts, U*ChangeRing(T0, Rationals()));
end for;

/* Transform back */
sigma := Sym(4) ! [1, 3, 2, 4];
P := PermutationMatrix(Integers(), sigma);
Ts := [ P*T : T in Ts ];

/* Sign */
for i in [1..#Ts] do
    T := Ts[i];
    E0 := T*E*Transpose(T);
    if Sign(E0[1,3]) lt 0 then
        T := DiagonalMatrix(Rationals(), [1,1,-1,1]) * T;
    end if;
    if Sign(E0[2,4]) lt 0 then
        T := DiagonalMatrix(Rationals(), [1,1,1,-1]) * T;
    end if;
    Ts[i] := T;
end for;
return Ts;

end function;
