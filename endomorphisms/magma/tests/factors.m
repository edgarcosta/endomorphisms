prec := 2000;
R<x> := PolynomialRing(RationalsExtra(prec));

fs := [
  10*x^10 + 24*x^9 + 23*x^8 + 48*x^7 + 35*x^6 + 35*x^4 - 48*x^3 + 23*x^2 - 24*x + 10,
  x^5 - x,
  x^8 + x^6 + 5*x^4 - 3*x^2 + 17,
  x^9 + 2*x
];

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



