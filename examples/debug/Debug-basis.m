SetVerbose("EndoFind", 3);
SetVerbose("CurveRec", 0);

function TransformForm(f, T : co := true, contra := false)
    R := Parent(f);
    vars := Matrix([ [ mon ] : mon in MonomialsOfDegree(R, 1) ]);
    if (not co) or contra then
        return Evaluate(f, Eltseq(ChangeRing(Transpose(T)^(-1), R) * vars));
    end if;
    return Evaluate(f, Eltseq(ChangeRing(T, R) * vars));
end function;

function RandomInvertibleMatrix(g, B)
D := [ -B..B ];
repeat
    T := Matrix(Rationals(), g,g, [ Random(D) : i in [1..g^2] ]);
until Determinant(T) eq 1;
return T;
end function;

function PolHom(f)
R := Parent(f);
S := PolynomialRing(BaseRing(R), 3);
g := (S.3^Degree(f))*Evaluate(f, [S.1/S.3, S.2/S.3]);
return S ! g;
end function;

function PeriodMatrixRetryQQ(f, CC);
while true do
    try
        RS := RiemannSurface(f : Precision := Precision(CC));
        P := ChangeRing(BigPeriodMatrix(RS), CC);
        return P, RS;
    catch e
        F := PolHom(f); F := TransformForm(F, RandomInvertibleMatrix(3, 2));
        X := PlaneCurve(F);
        f := DefiningPolynomial(AffinePatch(X, 1));
    end try;
end while;
end function;

prec := 300;
K := RationalsExtra(300);
R<x,y> := PolynomialRing(K, 2);
z := 1;
f := x^3*z+2*x^2*y^2+x^2*y*z+4*x^2*z^2+4*x*y^3+3*x*y^2*z-3*x*y*z^2+y^3*z-y^2*z^2-3*y*z^3-z^4;
RS := RiemannSurface(f : Precision := prec);
P := BigPeriodMatrix(RS);

R<x,y> := PolynomialRing(K, 2);
z := 1;
f := x^3*z+2*x^2*y^2+x^2*y*z+4*x^2*z^2+4*x*y^3+3*x*y^2*z-3*x*y*z^2+y^3*z-y^2*z^2-3*y*z^3-z^4;
RS := RiemannSurface(f : Precision := prec);
P := BigPeriodMatrix(RS);

T := RandomInvertibleMatrix(3, 2);
F := PolHom(f); G := TransformForm(F, T);
Y := PlaneCurve(G);
g := DefiningPolynomial(AffinePatch(Y, 1));
RS := RiemannSurface(g : Precision := prec);
Q := BigPeriodMatrix(RS);

TCC := EmbedMatrixExtra(T);
for U in [ TCC, TCC^(-1), Transpose(TCC), Transpose(TCC^(-1)) ] do
    //Q := ChangeRing(Q, ComplexFieldExtra(100));
    UCC := ChangeRing(U, BaseRing(Q));
    P0 := UCC*Q;
    I3CC := IdentityMatrix(BaseRing(Q), 3);
    try
        //print HomologyRepresentation(I3CC, P, P0);
        print SymplecticIsomorphismsCC(P, P0);
    catch e
        print "oink";
    end try;
end for;

exit;


prec := 300;
F := RationalsExtra(prec);
CC := F`CC;
R<x,y,z> := PolynomialRing(F, 3);
f := x^3*z+3*x^2*y^2+x^2*y*z+4*x^2*z^2+4*x*y^3+3*x*y^2*z-3*x*y*z^2+y^3*z-y^2*z^2-3*y*z^3-z^4;
X := PlaneCurve(f);

print "";
print "Curve:";
print X;

time P := PeriodMatrix(X);
time desc := HeuristicEndomorphismAlgebra(X : CC := true);
time rep := HeuristicEndomorphismRepresentation(X);
time L := HeuristicEndomorphismFieldOfDefinition(X);
time Lat := HeuristicEndomorphismLattice(X);
time test_gl2 := HeuristicIsGL2(X);
time dec := HeuristicDecomposition(X);

print "";
print "Heuristic endomorphism lattice:";
print Lat;

/*
facinfo := HeuristicJacobianFactors(X);
print "";
print "Heuristic Jacobian factors:";
print facinfo;

/*
exps, test, degs := IsogenyInformation(X : facinfo := facinfo);
print "";
print "Isogeny exponents:";
print exps;
print "";
print "Compatible with polarizations:";
print test;
print "";
print "Degrees (if applicable):";
print degs;
*/

