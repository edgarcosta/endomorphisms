// Blind curve-level upper-bound API coverage. Each row fixes B and pins the
// returned real-algebra token multiset without handing the API a target. The
// source curves and truth column come from corpus-cells.m; 3125.a1 documents
// that a small blind bound can be coarser than the geometric endomorphism algebra.
//
// Requires PARI/gp on PATH for center-field detection via Polredabs.
// Silent on success. Run from the root via: ./tests/run.sh blind

AttachSpec("../../endomorphisms/magma/spec");
AttachSpec(GetEnv("POLRED_SPEC"));

SetVerbose("EndoFind", 0);

R<x> := PolynomialRing(Rationals());

if #Pipe("command -v gp || true", "") eq 0 then
    print "blind: gp (PARI/gp) not on PATH; skipping.";
    exit;
end if;

// <label, B, blind bound, geometric truth, f-coeffs, h-coeffs (low to high)>.
cases := [*
    <"28561.c1", 10, ["CC", "CC"], ["CC", "CC"],
        [0, 0, -2, 2, 2, -3, -2], [1]>,
    <"3125.a1", 10, ["M_4(RR)"], ["CC", "CC"],
        [0, 0, 0, 0, 0, -1], [1]>,
    <"20736.ce1", 10, ["M_2(RR)"], ["M_2(RR)"],
        [0, -2, 0, 5, -3, -5, 1], [1, 1, 1, 1]>,
    <"8100.o1", 10, ["M_2(CC)"], ["M_2(RR)"],
        [1, 3, -42, 43, 21, -60, -28], [1]>,
    <"8192.d1", 10, ["M_2(RR)", "RR"], ["CC", "CC"],
        [-2, 0, 8, 0, -9, 0, 2], []>,
    <"729.a1", 10, ["M_2(CC)"], ["M_2(CC)"],
        [16, 0, 0, -5], [0, 0, 0, 1]>,
    <"529.a1", 10, ["RR", "RR"], ["RR", "RR"],
        [0, 1, 0, -1, 0, -1], [1, 0, 1, 1]>,
    <"121.a1", 10, ["M_2(RR)"], ["M_2(RR)"],
        [-2, 4, 2, 5, 2, 1], [1, 1, 1, 1]>,
    <"6859.a1", 10, ["CC", "RR"], ["CC", "RR"],
        [0, -1, -1, -11, -4, -16, 16], [1]>,
    <"363.a1", 10, ["RR", "RR"], ["RR", "RR"],
        [0, -1, 0, -2, 4, -2], [0, 1, 0, 1]>
*];

for entry in cases do
    label, B, expected, truth, fc, hc := Explode(entry);
    C := HyperellipticCurve(R ! fc, R ! hc);
    got := Sort(&cat RealRepresentationBound(C, B));
    error if got ne Sort(expected),
        Sprintf("%o at B = %o: expected blind bound %o, got %o",
                label, B, Sort(expected), got);
    if label eq "3125.a1" then
        error if got eq Sort(truth) or not RealRepresentationEmbeds(truth, got),
            Sprintf("%o at B = %o: expected a strictly coarser sound bound than %o, got %o",
                    label, B, Sort(truth), got);
    end if;
end for;

procedure test_tuple_28561()
    C := HyperellipticCurve(R ! [0, 0, -2, 2, 2, -3, -2], R ! [1]);
    ok, msg, eta, t, output, total_dim := EndomorphismAlgebraUpperBound(C, 10);
    error if not ok or msg ne "We have putatively computed eta and t. Under this assumption, we bounded the corresponding centers." or eta ne 2 or t ne 1 or total_dim ne 4,
        Sprintf("28561.c1 at B = 10: expected <true, message, 2, 1, output, 4>, got <%o, %o, %o, %o, %o, %o>",
                ok, msg, eta, t, output, total_dim);
    error if #output ne 1 or #output[1] ne 4 or output[1][1] ne 1 or output[1][2] ne 2 or
             #output[1][3] ne 2 or Degree(output[1][3][1]) ne 4 or #output[1][3][2] ne 3 or
             output[1][4] ne ["CC", "CC"],
        Sprintf("28561.c1 at B = 10: unexpected output component %o", output);
end procedure;

procedure test_tuple_3125()
    C := HyperellipticCurve(R ! [0, 0, 0, 0, 0, -1], R ! [1]);
    ok, msg, eta, t, output, total_dim := EndomorphismAlgebraUpperBound(C, 10);
    error if not ok or msg ne "We have putatively computed eta and t. Under this assumption, we bounded the corresponding centers." or eta ne 8 or t ne 1 or total_dim ne 16,
        Sprintf("3125.a1 at B = 10: expected <true, message, 8, 1, output, 16>, got <%o, %o, %o, %o, %o, %o>",
                ok, msg, eta, t, output, total_dim);
    error if #output ne 1 or #output[1] ne 4 or output[1][1] ne 4 or output[1][2] ne 2 or
             #output[1][3] ne 2 or Degree(output[1][3][1]) ne 1 or #output[1][3][2] ne 1 or
             output[1][4] ne ["M_4(RR)"],
        Sprintf("3125.a1 at B = 10: unexpected output component %o", output);
end procedure;

test_tuple_28561();
test_tuple_3125();
