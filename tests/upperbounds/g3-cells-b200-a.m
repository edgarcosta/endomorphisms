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


cells := [*
    <"3538944.3", 200, ["CC","CC","RR"], "x^3*y + x^3*z + x^2*y^2 + 2*x^2*y*z + x^2*z^2 + 2*x*y^3 + 2*x*z^3 + y^4 + z^4">,
    <"41600.1", 200, ["RR","RR","RR"], "x^3*z + x^2*y*z + x*y^3 - 2*x*y^2*z - x*y*z^2 - 2*x*z^3 - y^4 - y^3*z - 2*y^2*z^2 - y*z^3 - z^4">
*];

for entry in cells do
    label, maxB, expected, data := Explode(entry);
    C := BuildCurve(data);
    bound, B, status := RealRepresentationBound(C, expected : Bmax := maxB);
    error if status ne "sharp" or B gt maxB,
        Sprintf("%o: expected sharp at B <= %o; got %o at B = %o with bound %o",
                label, maxB, status, B, bound);
end for;
