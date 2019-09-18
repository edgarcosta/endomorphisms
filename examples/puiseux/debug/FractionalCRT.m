import "../../../endomorphisms/magma/puiseux/FractionalCRT.m": RandomSplitPrime, FractionalCRTQQ, FractionalCRTSplit, ReduceConstantSplit, ReduceConstantSplitAlt, ReduceCurveSplit;
import "../../../endomorphisms/magma/puiseux/Initialize.m": InitializeCurve;

prec := 1000;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
f := x^12 - 6*x^11 + 21*x^10 - 50*x^9 + 90*x^8 - 126*x^7 + 493*x^6 - 1182*x^5 - 2550*x^4 + 6990*x^3 - 2619*x^2 - 1062*x + 31329;
f := x^3 - 2;

/* This costs time if f has large degree */
B := 200;
time ps := [ RandomSplitPrime(f, B) : i in [1..10] ];

num := Random(Integers(), 2^B);
repeat
    den := Random(Integers(), 2^B);
until den ne 0;
x := num/den;

rs := [* FiniteField(Norm(p)) ! x : p in ps *];
/* This is very fast indeed */
time print x eq FractionalCRTQQ(rs, ps);


prec := 1000;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
F<r> := NumberField(x^3 - x - 1);
R<x> := PolynomialRing(F);
f := x^2 - 3;
K<r> := NumberField(f);

B := 100;
time ps := [ RandomSplitPrime(f, B) : i in [1..10] ];

num1 := Random(Integers(), 2^B);
num2 := Random(Integers(), 2^B);
num3 := Random(Integers(), 2^B);
repeat
    den := Random(Integers(), 2^B);
until den ne 0;
x := (num1 + num2*F.1 + num3*F.1^2)/den;

rs := [* *];
for p in ps do
    FF, h := ResidueClassField(p);
    Append(~rs, h(x));
end for;
time print x eq FractionalCRTSplit(rs, ps);


FF, h := ResidueClassField(ps[1]);
time bla := ReduceConstantSplit(x, h);
time bla := ReduceConstantSplit(x / 101, h);
time bla := ReduceConstantSplit(x / 123, h);


prec := 1000;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);
f := x^8 + x + 1; h := R ! 0;
X := HyperellipticCurve(f, h);
P0 := X ! [0, 1];
time InitializeCurve(X, P0);

FF, h := ResidueClassField(ps[1]);
time X0 := ReduceCurveSplit(X, h);


prec := 1000;
F := RationalsExtra(prec);
R<x,y,z> := PolynomialRing(F, 3);
D := [-5..5];
f := &+[ Random(D)*mon : mon in MonomialsOfDegree(R, 4) ];
f := f - Evaluate(f, [0, 0, 1])*z^4;
X := PlaneCurve(f);
P0 := X ! [0, 0, 1];
time InitializeCurve(X, P0);

FF, h := ResidueClassField(ps[1]);
time X0 := ReduceCurveSplit(X, h);

