/***
 *  Polarizations and Rosati involutions
 *
 *  Copyright (C) 2016-2017
 *            Jeroen Hanselman (jeroen.hanselman@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

intrinsic FindPolarizationBasis(P::.) -> .
{docstring}
// Basic invariants

JP :=ComplexStructure(P);
RR := BaseRing(JP);

gP := #Rows(JP) div 2;
n := 4 * gP ^2;
// Building a formal matrix corresponding to all possible polarisations
R := PolynomialRing(RR, n);
vars := GeneratorsSequence(R);
M := Matrix(R, 2 * gP, 2 *gP, vars);

JP_R := ChangeRing(JP, R);
// Conditions that ensure that E(ix,iy) = E(x,y) and that E is anti-symmetric. 
Comm := Eltseq(JP_R* M * Transpose(JP_R)- M) cat Eltseq(M+Transpose(M));

// Splitting previous linear equations by formal variable

M :=  Matrix(RR, [ [MonomialCoefficient(c, var) : c in Comm] :var in vars ]);
//print M;

Ker := IntegralLeftKernel(M);

// Culling the correct polarizations using the conditions on E.

RR:=BaseRing(JP); Pols := [];
for r in Rows(Ker) do    
    alpha := Matrix(Rationals(), 2*gP, 2*gP, Eltseq(r));
// Culling the correct polarizations using the conditions on E.    
    Comm:=JP*alpha*Transpose(JP) - alpha;   
    Comm2 := alpha+Transpose(alpha);
    if &and([Abs(c) lt RR`epscomp : c in Eltseq(Comm)]) then    
        if &and([Abs(c) lt RR`epscomp : c in Eltseq(Comm2)]) then 
            Append(~Pols, alpha);
        end if;    
    end if; 
end for;  
return Pols;
end intrinsic;


intrinsic FindSymplecticBasis(M::.) -> .
{docstring}

V:=SymplecticSpace(Matrix(M));
S:=HyperbolicSplitting(V);

B:=[];
n:=NumberOfRows(M)/2;
for i in [1..n] do       
    for v in S[1][i] do
        Append(~B, Vector(v)); 
    end for;
end for;
return Matrix(B);
end intrinsic;


intrinsic FindPrinPolPeriodMatrix(P::.) -> .
{docstring}

M:=FindPolarizationBasis(P)[1];
N:= FindSymplecticBasis(M);

for i in [1..NumberOfRows(N)] do        
    for j in [1..NumberOfRows(Transpose(N))] do
        N := N*Denominator(N[i][j]);       
    end for;
end for;
return N;
end intrinsic;
