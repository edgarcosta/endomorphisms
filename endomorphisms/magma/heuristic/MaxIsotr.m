/***
 *  Maximal isotropic subgroups and resulting lattices
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@unipi.it)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

forward IsogenousPPLatticesG2;


intrinsic InducedPolarization(E::., R::. : ProjToIdem := true) -> .
{Given a matrix E corresponding to a polarization, returns the pushforward (default) or pullback of E along R. The pullback is the pairing directly induced via the map represented by R, and the pushforward is its dual.}

if ProjToIdem then
    Q := (R*E^(-1)*Transpose(R))^(-1);
    Q *:= LCM([ Denominator(c) : c in Eltseq(Q) ]);
else
    Q := Transpose(R)*E*R;
end if;
return Q;

end intrinsic;


intrinsic FrobeniusFormAlternatingRational(E::.) -> .
{Version that works with matrices over QQ.}

E := ChangeRing(E, Integers());
E0, T0 := FrobeniusFormAlternating(E);
return ChangeRing(E0, Rationals()), ChangeRing(T0, Rationals());

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
/* Rows give vectors */
// TODO: Right cosets.

FF := FiniteField(p);
V := VectorSpace(FF, 2*d);
B0 := [ V.i : i in [1..d] ];
G := SymplecticGroup(2*d, FF);
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
// Based on a suggestion of John Voight; the dumb implementation is mine (...)
// TODO: Right cosets.

test, p, f := IsPrimePower(pf);
if f eq 1 and d eq 1 then
    return SymplecticSubmodulesPrime(p, d);
end if;

R := quo< Integers() | pf >;
M := RSpace(R, 2*d);

/* Find bases of rank d by enumeration */
bases := CartesianPower(M, 2*d);
bases := [ [ e : e in basis ] : basis in bases ];
bases := [ basis : basis in bases | sub< M | basis > eq M ];

/* Check with preserve new standard form E */
bases_new := [ ]; E := DiagonalJoin([ Matrix(R, [[0,1],[-1,0]]) : i in [1..d] ]);
for basis in bases do
    A := Matrix(R, [ Eltseq(b) : b in basis ]);
    /* Note that rows of A are important, hence difference from usual */
    if A*E*Transpose(A) eq E then
        Append(~bases_new, basis);
    end if;
end for;
bases := bases_new;

submods := [ ];
for basis in bases do
    /* Check that elements generate */
    factors := [ ];
    /* In the end we have to consider d pairs */
    for i in [1..d] do
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
end for;
return submods;

end function;


intrinsic SymplecticSubmodules(n::RngIntElt, d::RngIntElt) -> .
{All symplectic submodules of index n in rank 2*d. Uses rows.}

Fac := Factorization(n);
pfs := [ tup[1]^tup[2] : tup in Fac ];
L0 := Lattice(IdentityMatrix(Rationals(), 2*d));
Lats := [ L0 ];

for pf in pfs do
    submods := SymplecticSubmodulesPrimePower(pf, d);
    Latsnew := [ ];
    for Lat in Lats do
        for submod in submods do
            B := ChangeRing(Matrix(Basis(Lat)), Rationals());
            M := pf * B;
            Mnew := Matrix(Rationals(), [ [ Rationals() ! Integers() ! c : c in Eltseq(gen) ] : gen in Generators(submod) ]) * B;
            Append(~Latsnew, Lattice(VerticalJoin(M, Mnew)));
        end for;
    end for;
    Lats := Latsnew;
end for;
return Lats;

end intrinsic;


intrinsic SymplecticOvermodules(n::RngIntElt, d::RngIntElt) -> .
{All symplectic overmodules of index n in rank 2*d. Uses rows.}
// TODO: I forgot why this is mathematically correct...

Fac := Factorization(n);
pfs := [ tup[1]^tup[2] : tup in Fac ];
L0 := Lattice(IdentityMatrix(Rationals(), 2*d));
Lats := [ L0 ];

for pf in pfs do
    submods := SymplecticSubmodulesPrimePower(pf, d);
    Latsnew := [ ];
    for Lat in Lats do
        for submod in submods do
            M := ChangeRing(Matrix(Basis(Lat)), Rationals());;
            Mnew := (1/pf) * Matrix(Rationals(), [ [ Rationals() ! Integers() ! c : c in Eltseq(gen) ] : gen in Generators(submod) ]);
            Append(~Latsnew, Lattice(VerticalJoin(M, Mnew)));
        end for;
    end for;
    Lats := Latsnew;
end for;
return Lats;

end intrinsic;


intrinsic IsogenousPPLattices(E::. : ProjToPP := true) -> .
{Given an alternating form E, finds the sublattices to ZZ^2d of smallest possible index on which E induces a principal polarization. These are returned in matrix form, that is, as a span of a basis in the rows. This basis is symplectic in the usual sense. More concretely, applying T*E*Transpose(T) sends E to normal form, and T (resp. T^(-1)) is integral if PropToPP is false (resp. true).}

g := #Rows(E) div 2;
if g eq 1 then
    // TODO: Arguably we should return something here as well, since currently
    // the standard polarization may not be positive
    return [ IdentityMatrix(Rationals(), 2) ];
elif g eq 2 then
    return IsogenousPPLatticesG2(E : ProjToPP := ProjToPP);
else
    error "No functionality for principally polarized lattices in genus strictly larger than 2 yet";
end if;

end intrinsic;


function IsogenousPPLatticesG2(E : ProjToPP := true)
// TODO: Generalize

/* Take new normal form and find quotient of d_i */
E0, T0 := FrobeniusFormAlternatingAlt(E);
d := GCD([ Integers() ! c : c in Eltseq(E0) ]);
E0 := (1/d)*E0;
n := Integers() ! Abs(E0[3,4]);

/* Take relevant projections or inclusions */
if not ProjToPP then
    Lats := SymplecticSubmodules(n, 1);
    Us := [ ChangeRing(Matrix(Basis(Lat)), Rationals()) : Lat in Lats ];
    Us := [ DiagonalJoin(U, IdentityMatrix(Rationals(), 2)) : U in Us ];
else
    Lats := SymplecticOvermodules(n, 1);
    Us := [ ChangeRing(Matrix(Basis(Lat)), Rationals()) : Lat in Lats ];
    Us := [ DiagonalJoin(IdentityMatrix(Rationals(), 2), U) : U in Us ];
end if;

/* Correct sign and compose all transformations */
Ts := [ ];
for U in Us do
    E := U*E0*Transpose(U);
    _, T := FrobeniusFormAlternatingRational(E);
    Append(~Ts, T*U*T0);
end for;

/* Check distinctness of lattices obtained explicitly */
Lats := [ Lattice(T) : T in Ts ];
CP := CartesianPower([1..#Lats], 2);
for tup in [ tup : tup in CP | tup[1] ne tup[2] ] do
    assert Lats[tup[1]] ne Lats[tup[2]];
end for;
return Ts;

end function;
