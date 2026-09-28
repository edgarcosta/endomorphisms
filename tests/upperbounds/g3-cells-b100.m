// Frozen genus-3 coverage from the 4,341-quartic escalating corpus at 54f95f3.
// Equations are copied verbatim from artifacts/corpus/prepared/g3.tsv.

AttachSpec("../../endomorphisms/magma/spec");
AttachSpec(GetEnv("POLRED_SPEC"));

SetVerbose("EndoFind", 0);


P2<x, y, z> := ProjectiveSpace(Rationals(), 2);
P2CR := CoordinateRing(P2);

function BuildCurve(data)
    return Curve(P2, P2CR ! eval("return " cat data cat ";"));
end function;


cells := [*
    <"3326427.1", 128, ["CC"], "x^3*z + x*y^3 + y^4 + y^3*z + z^4">,
    <"131072.1", 64, ["CC","M_2(RR)"], "x^3*z + 2*x^2*y*z + 2*x^2*z^2 + x*y^3 - 2*x*y^2*z - 2*x*y*z^2 + x*z^3 - y^4 - y^3*z">,
    <"45056.1", 128, ["CC","RR"], "x^3*z + x^2*z^2 + 2*x*y^2*z + 2*x*y*z^2 + x*z^3 + y^4 + 2*y^3*z + y^2*z^2">,
    <"46080.1", 128, ["CC","RR","RR"], "x^4 + 2*x^3*y + 6*x^3*z - 4*x^2*y^2 - 2*x^2*y*z - 4*x^2*z^2 - x*y^3 - 3*x*y^2*z + 2*x*y*z^2 - 4*x*z^3 + 2*y^4 + 3*y^2*z^2 + 2*z^4">,
    <"655360.2", 64, ["M_2(CC)","RR"], "x^3*z + x^2*y^2 - 2*x*y^2*z - 2*x*z^3 - 2*y^4 - 2*y^2*z^2 - z^4">,
    <"6050.1", 64, ["RR","RR"], "x^3*z + x^2*y^2 + x*y^3 - x*y^2*z - 2*x*z^3 - y^2*z^2 - z^4">,
    <"35937.1", 128, ["RR","RR","RR"], "x^3*z + x^2*y^2 + x^2*y*z - x^2*z^2 + x*y^3 + x*y*z^2 + x*z^3 + y^3*z + y^2*z^2">
*];

for entry in cells do
    label, maxB, expected, data := Explode(entry);
    C := BuildCurve(data);
    bound, B, status := RealRepresentationBound(C, expected : Bmax := maxB);
    error if status ne "sharp" or B ne maxB,
        Sprintf("%o: expected sharp at B = %o; got %o at B = %o with bound %o",
                label, maxB, status, B, bound);
end for;
