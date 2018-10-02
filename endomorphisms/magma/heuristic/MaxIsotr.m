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


intrinsic InducedPolarization(E::., R::.) -> .
{bla}

Q := R*E*Transpose(R);
d := GCD([ Integers() ! c : c in Eltseq(Q) ]);
return Matrix(ChangeRing(Q/d, Integers()));

end intrinsic;


intrinsic FrobeniusFormAlternatingAlt(E::.) -> .
{T E Transpose(T) = E0.}

E1, T1 := FrobeniusFormAlternating(ChangeRing(E, Integers()));
g := #Rows(E) div 2; S := Sym(2*g);
sigma := S ! (&cat[ [ i, g + i ] : i in [1..g] ]);
P := PermutationMatrix(Integers(), sigma);
E2 := P*E1*Transpose(P);
T2 := P*T1;
return E2, T2;

end intrinsic;


function SymplecticSubmodulesPrimePower(pf, d)

assert d mod 2 eq 0;
test, p, f := IsPrimePower(pf);
assert test;

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


intrinsic SymplecticSubmodules(n::RngIntElt, d::RngIntElt) -> .
{All symplectic submodules.}

Fac := Factorization(n);
pfs := [ tup[1]^tup[2] : tup in Fac ];
L0 := Lattice(IdentityMatrix(Rationals(), d));
Lats := [ L0 ];

for pf in pfs do
    submods := SymplecticSubmodulesPrimePower(pf, d);
    Latsnew := [ ];
    for Lat in Lats do
        for submod in submods do
            //M := ChangeRing(Matrix(Basis(Lat)), Rationals());;
            //Mnew := (1/pf) * Matrix(Rationals(), [ [ Rationals() ! Integers() ! c : c in Eltseq(gen) ] : gen in Generators(submod) ]);
            B := ChangeRing(Matrix(Basis(Lat)), Rationals());
            M := pf * B;
            Mnew := Matrix(Rationals(), [ [ Integers() ! c : c in Eltseq(gen) ] : gen in Generators(submod) ]) * B;
            Append(~Latsnew, Lattice(VerticalJoin(M, Mnew)));
        end for;
    end for;
    Lats := Latsnew;
end for;
return Lats;

end intrinsic;


intrinsic IsogenousPPLatticesG2(E::.) -> .
{Finds isogenous lattices on which the polarization is principal.}

E0, T0 := FrobeniusFormAlternatingAlt(E);
n := Abs(E0[3,4]);

Ts := [ ];
for Lat in SymplecticSubmodules(n, 2) do
    Ur := ChangeRing(Matrix(Basis(Lat)), Rationals());
    //U := DiagonalJoin(IdentityMatrix(Rationals(), 2), Ur);
    U := DiagonalJoin(Ur, IdentityMatrix(Rationals(), 2));
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

end intrinsic;
