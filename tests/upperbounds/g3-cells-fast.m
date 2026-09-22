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
    <"19652.1", 10, ["M_2(RR)","RR"], "x^3*z + x^2*y*z + x^2*z^2 + x*y^3 + x*y^2*z + x*y*z^2 + x*z^3 + y^4 + y^3*z">,
    <"823543.1", 10, ["M_3(CC)"], "x^3*z + x*y^3 + y*z^3">,
    <"771147.1", 10, ["M_3(RR)"], "x^4 + 2*x^3*z + x^2*y^2 + x^2*y*z + 3*x^2*z^2 + x*y^2*z + x*y*z^2 + 2*x*z^3 + y^4 + 2*y^3*z + 3*y^2*z^2 + 2*y*z^3 + z^4">,
    <"27651.1", 10, ["RR","RR"], "x^4 + 2*x^3*y + 5*x^3*z - 4*x^2*y^2 - 3*x^2*y*z - 2*x^2*z^2 - x*y^3 - 4*x*y^2*z - 3*x*y*z^2 - 6*x*z^3 + 2*y^4 + 2*y^3*z + 4*y^2*z^2 + 2*y*z^3 + 3*z^4">,
    <"480491.1", 10, ["RR","RR","RR"], "x^3*z + x^2*y^2 + x^2*y*z + x^2*z^2 - x*y^3 + x*y^2*z + x*z^3 - 2*y^3*z + y*z^3">,
    <"492075.1", 20, ["CC"], "x^3*z + x^2*z^2 + x*y^3 - x*z^3 + y^3*z">,
    <"2519424.1", 20, ["CC","M_2(CC)"], "x^3*z + x*y^3 + x*z^3 + y^3*z">,
    <"177147.1", 20, ["CC","M_2(RR)"], "x^3*z + x^2*z^2 + x*y^3 + x*z^3 + y^3*z">,
    <"90112.1", 20, ["CC","RR"], "x^3*y + x^3*z + x^2*y*z + 2*x^2*z^2 + x*y^3 + x*y^2*z + y^3*z + 2*y^2*z^2 - z^4">,
    <"843750.1", 20, ["CC","RR","RR"], "x^3*y + x^3*z + x^2*y^2 - 2*x^2*z^2 - x*y^3 - x*y^2*z - 2*x*y*z^2 - y^3*z - y^2*z^2 - y*z^3 + z^4">,
    <"16000.1", 20, ["M_2(CC)","RR"], "x^3*z + x^2*y^2 + x*y^3 - 2*x*y^2*z + x*z^3 - y^3*z + y^2*z^2">,
    <"2940.1", 20, ["M_2(RR)","RR"], "x^3*y + x^3*z + x^2*y^2 + 3*x^2*y*z + x^2*z^2 - 4*x*y^3 - 3*x*y^2*z - 3*x*y*z^2 - 4*x*z^3 + 2*y^4 + 3*y^2*z^2 + 2*z^4">,
    <"1830125.1", 20, ["M_3(RR)"], "x^4 + 2*x^3*z + 3*x^2*y^2 + 3*x^2*y*z + 3*x^2*z^2 + 3*x*y^2*z + 3*x*y*z^2 + 2*x*z^3 + y^4 + 2*y^3*z + 3*y^2*z^2 + 2*y*z^3 + z^4">,
    <"5978.1", 20, ["RR","RR"], "x^3*z + x^2*y^2 + x^2*y*z + x*y^3 + x*y^2*z + x*y*z^2 + x*z^3 + y^3*z + y^2*z^2">,
    <"48778.1", 20, ["RR","RR","RR"], "x^3*z + x^2*y*z + x*y^3 + x*y^2*z + x*y*z^2 + x*z^3 + y^4 + y^3*z">,
    <"2834352.1", 50, ["CC"], "x^3*z + y^4 + y^3*z + 2*y^2*z^2 + y*z^3 + z^4">,
    <"16000.2", 50, ["M_2(CC)","RR"], "x^3*z + x*y^3 + x*z^3 + y^4 + y^3*z">,
    <"26325.1", 50, ["M_2(RR)","RR"], "x^4 + 2*x^3*z + 3*x^2*y^2 + 3*x^2*y*z + 4*x^2*z^2 + 3*x*y^2*z + 3*x*y*z^2 + 3*x*z^3 + y^4 + 2*y^3*z + 4*y^2*z^2 + 3*y*z^3 + 2*z^4">,
    <"5835.1", 50, ["RR","RR"], "x^4 + 2*x^3*y + 2*x^3*z - 4*x^2*y^2 + 2*x^2*y*z - 4*x^2*z^2 - x*y^3 - x*z^3 + 2*y^4 - 3*y^3*z + 5*y^2*z^2 - 3*y*z^3 + 2*z^4">,
    <"27951.1", 50, ["RR","RR","RR"], "x^3*z + x^2*y*z + 3*x^2*z^2 + x*y^3 - 3*x*y^2*z + x*y*z^2 - x*z^3 + 2*y^4 + y^2*z^2">,
    <"73728.1", 50, ["CC","M_2(RR)"], "x^3*z + x^2*y^2 + 2*x^2*y*z + 2*x*y^3 - 2*x*y*z^2 + x*z^3 + y^4 - 2*y^3*z + y^2*z^2">,
    <"27702.1", 50, ["CC","RR"], "x^3*z + x^2*y^2 + x^2*y*z + x^2*z^2 + x*y^2*z + x*y*z^2 + x*z^3 + y^3*z + y^2*z^2 + y*z^3">,
    <"139264.1", 50, ["CC","RR","RR"], "x^4 - x^2*y^2 - x^2*z^2 + x*y^3 + x*y^2*z + x*y*z^2 + x*z^3 + y^4 + y^3*z + 2*y^2*z^2 + y*z^3 + z^4">
*];

for entry in cells do
    label, maxB, expected, data := Explode(entry);
    C := BuildCurve(data);
    bound, B, status := RealRepresentationBound(C, expected : Bmax := maxB);
    error if status ne "sharp" or B gt maxB,
        Sprintf("%o: expected sharp at B <= %o; got %o at B = %o with bound %o",
                label, maxB, status, B, bound);
end for;
