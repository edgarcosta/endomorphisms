printf "Testing the Genus 2 curve 262144.d.524288.2...";

prec := 50;

F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
f := -25*x^6 + 12*x^5 + 27*x^4 - 16*x^3 - 3*x^2 + 4*x + 1;
X := HyperellipticCurve(f);


P := PeriodMatrix(X);
desc := HeuristicEndomorphismAlgebra(X : CC := true);
assert desc eq <[ <2, 1> ], [<1, 4, [ -1, 1 ], 6, 2>], <1, 1>, 3>;
rep := HeuristicEndomorphismRepresentation(X);
assert #rep eq 1;
L := HeuristicEndomorphismFieldOfDefinition(X);
assert Coefficients(DefiningPolynomial(L)) eq [ 1, 0, 0, 0, 1 ];
Lat := HeuristicEndomorphismLattice(X);
assert #Lat[3] eq 5;
test_gl2 := HeuristicIsGL2(X);
assert not test_gl2;
dec := HeuristicDecomposition(X);
assert dec eq [* Rationals(), [* L, true *], [* [ [ 2, 1 ] ], [] *], [* [ [ 2, 1 ] ], [], [] *] *];

assert [] eq HeuristicDecompositionFactors(X);
