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
    <"79488.1", 128, ["CC","RR"], "x^3*z + x^2*y^2 + 4*x^2*z^2 + x*y^3 + 4*x*y^2*z + x*z^3 + 2*y^4 - y^3*z + y^2*z^2">,
    <"1053696.1", 128, ["CC","RR","RR"], "x^4 + x^3*z + 5*x^2*y^2 + 2*x^2*z^2 + 3*x*y^2*z + x*z^3 + 7*y^4 + 5*y^2*z^2 + z^4">
*];

for entry in cells do
    label, maxB, expected, data := Explode(entry);
    C := BuildCurve(data);
    bound, B, status := RealRepresentationBound(C, expected : Bmax := maxB);
    error if status ne "sharp" or B ne maxB,
        Sprintf("%o: expected sharp at B = %o; got %o at B = %o with bound %o",
                label, maxB, status, B, bound);
end for;
