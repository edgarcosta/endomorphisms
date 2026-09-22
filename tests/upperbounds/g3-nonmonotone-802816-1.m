// Frozen genus-3 coverage from the 4,341-quartic escalating corpus at 54f95f3.
// Equations are copied verbatim from artifacts/corpus/prepared/g3.tsv.

AttachSpec("../../endomorphisms/magma/spec");
AttachSpec(GetEnv("POLRED_SPEC"));

SetVerbose("EndoFind", 0);

if #Pipe("command -v gp || true", "") eq 0 then
    print "g3-cells: gp (PARI/gp) not on PATH; skipping.";
    exit;
end if;

P2<x, y, z> := ProjectiveSpace(Rationals(), 2);
P2CR := CoordinateRing(P2);

function BuildCurve(data)
    return Curve(P2, P2CR ! eval("return " cat data cat ";"));
end function;


nonmonotone := [*
    <"802816.1", "x^3*z - x*y^2*z - x*z^3 + y^4 + y^2*z^2">,
*];

for entry in nonmonotone do
    label, data := Explode(entry);
    C := BuildCurve(data);
    error if Sort(&cat RealRepresentationBound(C, 50)) ne ["M_2(RR)", "RR"],
        Sprintf("%o: B = 50 did not produce [M_2(RR), RR]", label);
    error if Sort(&cat RealRepresentationBound(C, 200)) ne ["CC", "CC", "M_2(RR)"],
        Sprintf("%o: B = 200 did not produce [CC, CC, M_2(RR)]", label);
    error if Sort(&cat RealRepresentationBound(C, 400)) ne ["RR", "RR", "RR"],
        Sprintf("%o: B = 400 did not recover [RR, RR, RR]", label);
end for;
