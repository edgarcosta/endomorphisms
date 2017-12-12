/***
 *  Riemann-Roch functionality
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


forward ExtractHomomorphismsField;
forward CandidateDivisorsNew;
forward IrreducibleComponentsFromBranchesNew;
forward DivisorFromMatrixByDegreeNew;


function ExtractHomomorphismsField(X, Y)

KX := X`K; KY := Y`K;
varord := VariableOrder();
// TODO: Test other orderings
Kprod := RationalFunctionField(X`F, 4);
seqX := [ Kprod.Index(varord, i) : i in [1..2] ];
seqY := [ Kprod.Index(varord, i) : i in [3..4] ];
hX := hom<KX -> Kprod | seqX >;
hY := hom<KY -> Kprod | seqY >;
return [ hX, hY ];

end function;


function CandidateDivisorsNew(X, Y, d)
/*
 * Input:   Two curves X and Y and a degree d.
 * Output:  Equations for divisors of degree d coming from the ambient of X.
 */

gX := X`g; fX := X`DEs[1]; RX := X`R; KX := X`K; P := X`P0;
gY := Y`g; fY := Y`DEs[1]; RY := Y`R; KY := Y`K; Q := Y`P0;
V, phiV := RiemannRochSpace((d + 2*gX)*Divisor(P));
V, phiV := RiemannRochSpace((d + gX)*Divisor(P));
divsX := [ KX ! phiV(v) : v in Basis(V) ];
W, phiW := RiemannRochSpace(3*gY*Divisor(Q));
W, phiW := RiemannRochSpace(2*gY*Divisor(Q));
divsY := [ KY ! phiW(w) : w in Basis(W) ];
divsX := [ divsX[#divsX] ] cat divsX[1..(#divsX - 1)];
divsY := [ divsY[#divsY] ] cat divsY[1..(#divsY - 1)];

hs := ExtractHomomorphismsField(X, Y);
CP := [ [* divX, divY *] : divX in divsX, divY in divsY ];
divs := [ &*[ hs[i](tup[i]) : i in [1..2] ] : tup in CP ];
return divs;

end function;


function IrreducibleComponentsFromBranchesNew(X, Y, fs, P, Qs : DivPP1 := false)
/*
 * Input:   Two curves X and Y,
 *          a basis of divisor equations fs,
 *          the precision n used when determining these,
 *          and branch expansions P and Qs.
 * Output:  The irreducible components that fit the given data.
 */

/* Recovering a linear system */
e := Maximum([ Maximum([ Denominator(Valuation(c - Coefficient(c, 0))) : c in Q ]) : Q in Qs ]);
prec := Precision(Parent(P[1]));
evss := [ [ Evaluate(f, ExtractPoints(X, Y, P, Q)) : Q in Qs ] : f in fs ];
min := Minimum([ Valuation(ev) : ev in &cat(evss) ]);
max := Minimum([ AbsolutePrecision(ev) : ev in &cat(evss) ]);
M := Matrix([ &cat[ [ Coefficient(ev, i/e) : i in [(e*min)..(e*max - 10)] ] : ev in evs ] : evs in evss ]);
min := -9;
max := 300/2;
M := Matrix([ &cat[ [ Coefficient(ev, i/e) : i in [(e*min)..(e*max)] ] : ev in evs ] : evs in evss ]);
B := Basis(Kernel(M));
/* Coerce back to ground field (possible because of echelon form) */
B := [ [ X`F ! c : c in Eltseq(b) ] : b in B ];
print "Dimension:";
print #B;
return 0;

/* Corresponding equations */
hX, hY := Explode(ExtractHomomorphismsRing(X, Y));
Rprod := Codomain(hX);
Aprod := AffineSpace(Rprod);
eqs := [ &+[ b[i] * fs[i] : i in [1..#fs] ] : b in B ];
eqs := [ Rprod ! Numerator(&+[ b[i] * fs[i] : i in [1..#fs] ]) : b in B ];
S := Scheme(Aprod, [ hX(DE) : DE in X`DEs ] cat [ hY(DE) : DE in Y`DEs ] cat eqs);

I := DefiningIdeal(S);
varord := VariableOrder();
J1 := ideal<Rprod | [ hX(DE) : DE in X`DEs ] cat [ Rprod.Index(varord, i + 2) - (X`P0)[i] : i in [1..2] ]>;
J2 := ideal<Rprod | [ hY(DE) : DE in Y`DEs ] cat [ Rprod.Index(varord, i) - (Y`P0)[i] : i in [1..2] ]>;
vprintf EndoCheck, 3 : "Calculating colon ideals... ";
repeat
    Iold := I;
    I := ColonIdeal(I, J1);
until I eq Iold;
repeat
    Iold := I;
    I := ColonIdeal(I, J2);
until I eq Iold;
vprintf EndoCheck, 3 : "done.\n";
S := Scheme(Aprod, I);
S := ReducedSubscheme(S);
vprint EndoCheck, 3 : "Dimension:";
vprint EndoCheck, 3 : Dimension(S);
vprint EndoCheck, 3 : "Degree and dimensions to factors:";
vprint EndoCheck, 3 : BiDimDeg(X, X, S);
/*
for I in IrreducibleComponents(S) do
    print BiDimDeg(X, X, I);
end for;
*/

/*
if DivPP1 then
*/

/* Corresponding scheme */
eqs := [ &+[ b[i] * fs[i] : i in [1..#fs] ] : b in B ];
Kprod := FieldOfFractions(Rprod);
Sprod := Scheme(Aprod, [ hX(DE) : DE in X`DEs ] cat [ hY(DE) : DE in Y`DEs ]);
Ss := [ ];
for f in eqs do
    Scl := ProjectiveClosure(Sprod);
    Rcl := CoordinateRing(Ambient(Scl));
    Ccl := FieldOfFractions(Rcl);
    h := hom< Kprod -> Ccl | [ Ccl.1, Ccl.2, Ccl.3, Ccl.4 ] >;
    h := hom< Rprod -> Rcl | [ Rcl.1, Rcl.2, Rcl.3, Rcl.4 ] >;
    S_num := Divisor(Scl, h(Numerator(f)));
    S_den := Divisor(Scl, h(Denominator(f)));
    S := S_num - S_den;
    Append(~Ss, S);
end for;
S := GCD(Ss);
vprint EndoCheck, 3 : Support(SignDecomposition(S));
vprint EndoCheck, 3 : "Dimension:";
vprint EndoCheck, 3 : Dimension(S);
vprint EndoCheck, 3 : BiDimDeg(X, X, S);

return 0;
return [ S ];

end function;


function DivisorFromMatrixByDegreeNew(X, Y, NormM, d : Margin := 2^4, DivPP1 := false, have_to_check := true)

vprintf EndoCheck, 2 : "Trying degree %o...\n", d;
fs := CandidateDivisorsNew(X, Y, d);
n := #fs + Margin;
vprintf EndoCheck, 2 : "Number of terms in expansion: %o.\n", n;

/* Take non-zero image branch */
vprintf EndoCheck, 2 : "Expanding... ";
P, Qs := ApproximationsFromTangentAction(X, Y, NormM, n);
vprintf EndoCheck, 4 : "Base point:\n";
_<t> := Parent(P[1]);
_<r> := BaseRing(Parent(P[1]));
vprint EndoCheck, 4 : P;
vprintf EndoCheck, 4 : "Resulting branches:\n";
vprint EndoCheck, 4 : Qs;
vprint EndoCheck, 4 : BaseRing(Parent(P[1]));
vprintf EndoCheck, 2 : "done.\n";

/* Fit a divisor to it */
vprintf EndoCheck, 2 : "Solving linear system... ";
ICs := IrreducibleComponentsFromBranchesNew(X, Y, fs, P, Qs : DivPP1 := DivPP1);
vprintf EndoCheck, 2 : "done.\n";

for S in ICs do
    DEs := DefiningEquations(S);
    vprintf EndoCheck, 2 : "Checking:\n";
    vprintf EndoCheck, 2 : "Step 1... ";
    test1 := CheckEquations(X, Y, P, Qs, DEs);
    vprintf EndoCheck, 2 : "done.\n";
    if test1 then
        vprintf EndoCheck, 2 : "Step 2... ";
        test2 := CheckIrreducibleComponent(X, Y, S);
        vprintf EndoCheck, 2 : "done.\n";
        if test2 then
            vprintf EndoCheck, 2 : "Divisor found!\n";
            return true, S;
        end if;
    end if;
end for;
return false, [ ];

end function;
