printf "Testing  HeuristicEndomorphismAlgebra/Representation/Lattice over Qsqrt5/-3...";
prec := 50;

QQ := RationalsExtra(prec);
_<t> := PolynomialRing(QQ);
F<r> := NumberField(t^2 - t - 1);
R<x> := PolynomialRing(F);

fs := [*
  x^6 + Sqrt(F!5),
  x^6 + r, x^5 + x + 1 *];

F<r> := BaseNumberFieldExtra(t^2 - t + 1, prec);
R<x> := PolynomialRing(F);

fs cat:= [*
  x^6 + r,
  x^5 + r*x^3 + x,
  R ! [ -30*r + 42, -156*r + 312, -66*r + 186, -1456*r + 1040, -90*r + 126, 156*r - 312, -22*r + 62 ]
*];

// curves := [* HyperellipticCurve(f) : f in fs *];
// [Dimension(A) where _, A := HeuristicEndomorphismAlgebra(C : CC:=true) : C in curves];
dimensions := [ 8, 8, 1, 8, 4, 4 ];
// [HeuristicEndomorphismLattice(C) : C in curves]
// [<Eltseq(DefiningPolynomial(elt[1])), <Eltseq(DefiningPolynomial(elt[2,1])),elt[2,2]>, elt[3,1], elt[4,1], elt[4,3]> where elt := HeuristicDecomposition(C) : C in curves]

for i->f in fs do
  printf ".";
  C := HyperellipticCurve(f);
  _, A := HeuristicEndomorphismAlgebra(X : CC := true);
  assert Dimension(A) eq dimensions[i];
  assert HeuristicEndomorphismLattice(C) eq lattices[i];
  elt := HeuristicDecomposition(C);
  assert <Eltseq(DefiningPolynomial(elt[1])), <Eltseq(DefiningPolynomial(elt[2,1])),elt[2,2]>, elt[3,1], elt[4,1], elt[4,3]> eq decompositions[i];
end for;
