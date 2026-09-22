// Frozen genus-3 coverage from the 4,341-quartic escalating corpus at 54f95f3.
// The default remains below 60 s: all B <= 20 cells, five B = 50 cells, and
// the three non-monotone witnesses at B = 50 and B = 200. Set
// ENDO_G3_FULL_TESTS=1 to run all 35 smallest-conductor cells and the B = 400
// recovery. UpperBounds.m:68-78 explains why B200 can be looser than B50.
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
full := GetEnv("ENDO_G3_FULL_TESTS") eq "1";

function BuildCurve(data)
    return Curve(P2, P2CR ! eval("return " cat data cat ";"));
end function;

// <label, recorded B, expected tokens, plane-quartic equation>.
default_cells := [*
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
    <"27951.1", 50, ["RR","RR","RR"], "x^3*z + x^2*y*z + 3*x^2*z^2 + x*y^3 - 3*x*y^2*z + x*y*z^2 - x*z^3 + 2*y^4 + y^2*z^2">
*];

full_only_cells := [*
    <"73728.1", 50, ["CC","M_2(RR)"], "x^3*z + x^2*y^2 + 2*x^2*y*z + 2*x*y^3 - 2*x*y*z^2 + x*z^3 + y^4 - 2*y^3*z + y^2*z^2">,
    <"27702.1", 50, ["CC","RR"], "x^3*z + x^2*y^2 + x^2*y*z + x^2*z^2 + x*y^2*z + x*y*z^2 + x*z^3 + y^3*z + y^2*z^2 + y*z^3">,
    <"139264.1", 50, ["CC","RR","RR"], "x^4 - x^2*y^2 - x^2*z^2 + x*y^3 + x*y^2*z + x*y*z^2 + x*z^3 + y^4 + y^3*z + 2*y^2*z^2 + y*z^3 + z^4">,
    <"3326427.1", 100, ["CC"], "x^3*z + x*y^3 + y^4 + y^3*z + z^4">,
    <"131072.1", 100, ["CC","M_2(RR)"], "x^3*z + 2*x^2*y*z + 2*x^2*z^2 + x*y^3 - 2*x*y^2*z - 2*x*y*z^2 + x*z^3 - y^4 - y^3*z">,
    <"45056.1", 100, ["CC","RR"], "x^3*z + x^2*z^2 + 2*x*y^2*z + 2*x*y*z^2 + x*z^3 + y^4 + 2*y^3*z + y^2*z^2">,
    <"46080.1", 100, ["CC","RR","RR"], "x^4 + 2*x^3*y + 6*x^3*z - 4*x^2*y^2 - 2*x^2*y*z - 4*x^2*z^2 - x*y^3 - 3*x*y^2*z + 2*x*y*z^2 - 4*x*z^3 + 2*y^4 + 3*y^2*z^2 + 2*z^4">,
    <"655360.2", 100, ["M_2(CC)","RR"], "x^3*z + x^2*y^2 - 2*x*y^2*z - 2*x*z^3 - 2*y^4 - 2*y^2*z^2 - z^4">,
    <"6050.1", 100, ["RR","RR"], "x^3*z + x^2*y^2 + x*y^3 - x*y^2*z - 2*x*z^3 - y^2*z^2 - z^4">,
    <"35937.1", 100, ["RR","RR","RR"], "x^3*z + x^2*y^2 + x^2*y*z - x^2*z^2 + x*y^3 + x*y*z^2 + x*z^3 + y^3*z + y^2*z^2">,
    <"3538944.3", 200, ["CC","CC","RR"], "x^3*y + x^3*z + x^2*y^2 + 2*x^2*y*z + x^2*z^2 + 2*x*y^3 + 2*x*z^3 + y^4 + z^4">,
    <"79488.1", 200, ["CC","RR"], "x^3*z + x^2*y^2 + 4*x^2*z^2 + x*y^3 + 4*x*y^2*z + x*z^3 + 2*y^4 - y^3*z + y^2*z^2">,
    <"1053696.1", 200, ["CC","RR","RR"], "x^4 + x^3*z + 5*x^2*y^2 + 2*x^2*z^2 + 3*x*y^2*z + x*z^3 + 7*y^4 + 5*y^2*z^2 + z^4">,
    <"41600.1", 200, ["RR","RR","RR"], "x^3*z + x^2*y*z + x*y^3 - 2*x*y^2*z - x*y*z^2 - 2*x*z^3 - y^4 - y^3*z - 2*y^2*z^2 - y*z^3 - z^4">,
    <"111537.1", 400, ["RR","RR","RR"], "x^3*z + x^2*y*z + x^2*z^2 - x*y^3 - x*y*z^2 - x*z^3 + y^2*z^2">
*];

cells := full select default_cells cat full_only_cells else default_cells;
for entry in cells do
    label, maxB, expected, data := Explode(entry);
    C := BuildCurve(data);
    bound, B, status := RealRepresentationBound(C, expected : Bmax := maxB);
    error if status ne "sharp" or B gt maxB,
        Sprintf("%o: expected sharp at B <= %o; got %o at B = %o with bound %o",
                label, maxB, status, B, bound);
end for;

// These three models become less precise between the B50 and B200 calls.
nonmonotone := [*
    <"802816.1", "x^3*z - x*y^2*z - x*z^3 + y^4 + y^2*z^2">,
    <"802816.2", "x^3*z + x*y^2*z - x*z^3 + y^4 - y^2*z^2">,
    <"3781323.1", "x^4 + 4*x^3*y + 4*x^3*z + 6*x^2*y^2 + 3*x^2*y*z + 6*x^2*z^2 + 4*x*y^3 - 5*x*y^2*z - 5*x*y*z^2 + 4*x*z^3 + y^4 - 5*y^3*z + 7*y^2*z^2 - 5*y*z^3 + z^4">
*];

for entry in nonmonotone do
    label, data := Explode(entry);
    C := BuildCurve(data);
    error if Sort(&cat RealRepresentationBound(C, 50)) ne ["M_2(RR)", "RR"],
        Sprintf("%o: B = 50 did not produce [M_2(RR), RR]", label);
    error if Sort(&cat RealRepresentationBound(C, 200)) ne ["CC", "CC", "M_2(RR)"],
        Sprintf("%o: B = 200 did not produce [CC, CC, M_2(RR)]", label);
    if full then
        error if Sort(&cat RealRepresentationBound(C, 400)) ne ["RR", "RR", "RR"],
            Sprintf("%o: B = 400 did not recover [RR, RR, RR]", label);
    end if;
end for;
