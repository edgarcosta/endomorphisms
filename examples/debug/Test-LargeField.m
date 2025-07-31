SetVerbose("EndoFind", 0);

F := RationalsExtra(1000);
R<x,y,z> := PolynomialRing(Rationals(), 3);

fs := [
-18*x^4+27*x^3*z-27*x^2*z^2-9*x*z^3+y^4-9*z^4,
-24*x^4+16*x^3*z-24*x^2*z^2-24*x*z^3+y^4-14*z^4,
-14*x^4-28*x^3*z-21*x^2*z^2-7*x*z^3+y^4-29*z^4,
14*x^4+28*x^3*z+21*x^2*z^2+7*x*z^3+y^4+29*z^4,
9*x^4-24*x^3*z-18*x^2*z^2+8*x*z^3+y^4+z^4,
-4*x^4-8*x^3*z+12*x^2*z^2+4*x*z^3+y^4-z^4,
-28*x^3*z+27*x*z^3+y^4,
-3*x^4-9*x^3*z-18*x^2*z^2-30*x*z^3+y^4-24*z^4,
-25*x^4-25*x^3*z-15*x^2*z^2+5*x*z^3+y^4-4*z^4,
-4*x^4+5*x^3*z-15*x^2*z^2-25*x*z^3+y^4-25*z^4,
25*x^4+25*x^3*z+15*x^2*z^2-5*x*z^3+y^4+4*z^4,
4*x^4-5*x^3*z+15*x^2*z^2+25*x*z^3+y^4+25*z^4,
-2*x^4-16*x^3*z+24*x^2*z^2+8*x*z^3+y^4-5*z^4,
2*x^4+16*x^3*z-24*x^2*z^2-8*x*z^3+y^4+5*z^4,
16*x^4+24*x^3*z+24*x^2*z^2+20*x*z^3+y^4+8*z^4,
-2*x^4+3*x^3*z-3*x^2*z^2-x*z^3+y^4-z^4,
-2*x^4+11*x^3*z-24*x^2*z^2+22*x*z^3+y^4-8*z^4,
x^4+3*x^3*z+6*x^2*z^2+10*x*z^3+y^4+8*z^4
];

for f in fs do
    print f;
    X := _PlaneCurve(f);
    L := HeuristicEndomorphismFieldOfDefinition(X);
    print L;
end for;

