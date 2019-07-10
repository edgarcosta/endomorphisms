//SetVerbose("EndoFind", 1);
//SetVerbose("CurveRec", 2);

prec := 300;

F := RationalsExtra(prec);
P3<x,y,z,w> :=ProjectiveSpace(F, 3);

f1 := -y*z - 12*z^2 + x*w - 32*w^2;
f2 := y^3 + 108*x^2*z + 36*y^2*z + 8208*x*z^2 - 6480*y*z^2 + 74304*z^3 + 96*y^2*w
+ 2304*y*z*w - 248832*z^2*w + 2928*y*w^2 - 75456*z*w^2 + 27584*w^3;
X := Curve(P3, [f1, f2]);


print "";
print "Curve:";
print X;

Lat := HeuristicEndomorphismLattice(X);
print "";
print "Heuristic endomorphism lattice:";
print Lat;

test_gl2 := HeuristicIsGL2(X);
print "";
print "Heuristic GL_2-determination:";
print test_gl2;
