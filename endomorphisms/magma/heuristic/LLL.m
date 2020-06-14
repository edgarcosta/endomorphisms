/***
 *  Our version of Magma functionality
 */

function lindep(A,H,DELTA) n:=#A; m:=Max([Abs(a) : a in A]); H:=Ceiling(H);
if m eq 0 then return Polynomial([1]); end if; A:=[a/m : a in A];
M:=Matrix([[i eq j select 1 else 0 : i in [1..n]] cat
           [H^n*Real(A[j]),H^n*Imaginary(A[j])] : j in [1..n]]);
B,T:=LLL(M : Delta:=DELTA);
return Eltseq(T[1]);
end function;

function OurIntegerRelation(A, H : Delta:=0.75)
//{Given a sequence of complex numbers and a height bound,
// return the best linear relation up to that height, using LLL with Delta}
CC := Universe(A);
v:=lindep(A,H,Delta);
abs := Abs(&+[ A[i]*v[i] : i in [1..#A] ]);
/* Only a weak test, since root is evaluated later, which is the better test */
if abs gt 10^20*CC`epscomp then return false, 0; end if;
return true, v;
end function;

function OurAllLinearRelations(q, p);
//Given a sequence with entries from a real field, return the lattice of all
//(small) integer linear dependencies among the entries. The precision p is
//used for two purposes, namely a relation must be zero to approximately 10^(-p),
//and secondly the coefficients themselves must be bounded approximately by p.
U:=Universe(q);
if Type(U) eq FldCom then
    L1 := OurAllLinearRelations([ Real(x) : x in q ], p);
    L2 := OurAllLinearRelations([ Imaginary(x) : x in q ], p);
    return L1 meet L2;
end if;
I:=[1..#q]; hb := U`height_bound;
M:=Matrix([[q[j]] cat [i eq j select 10^(-p) else 0 : i in I] : j in I]);
A,B:=LLL(M);
E:=Rows(B)[[i : i in I | Norm(A[i]) lt hb*p^2*#q^2*10^(-2*p)]];
return Lattice(Matrix(E));
end function;

function OurAlgdep(q, p);
CC := Universe(q);
M := Transpose(Matrix([ [ Real(c) : c in q ], [ Im(c) : c in q ] ]));
B := 10^p;

MJ := Matrix(Integers(), [ [ Round(B * c) : c in Eltseq(row) ] : row in Rows(M) ]);
MI := IdentityMatrix(Integers(), #Rows(MJ)); MJ := HorizontalJoin(MI, MJ);
L, K := LLL(MJ); rowsK := Rows(K);
//print rowsK;

row1 := rowsK[1]; ht1 := Max([ Height(c) : c in Eltseq(row1) ]);
abs := Abs(&+[ row1[i]*q[i] : i in [1..#q] ]);
/* Only a weak test, since root is evaluated later, which is the better test */
test1 := ht1 lt CC`height_bound;
test2 := abs lt CC`epscomp;
if test1 and test2 then
    return true, Eltseq(row1);
end if;
return false, Eltseq(0*row1);
end function;

function OurAlgdepNew(q, p);
CC := Universe(q); assert Type(CC) eq FldCom;
L := OurAllLinearRelations(q, p);
if Rank(L) eq 0 then return false, 0; end if;
s := Eltseq(Basis(L)[1]);
abs := Abs(&+[ s[i]*q[i] : i in [1..#q] ]);
/* Only a weak test, since root is evaluated later, which is the better test */
if abs gt CC`epscomp then return false, 0; end if;
return true, s;
end function;
