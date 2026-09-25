// Regression tests for ambiguous real-representation tokens.
// Silent on success. Run from the repository root via: ./tests/run.sh possibilities

AttachSpec("../../endomorphisms/magma/spec");
AttachSpec(GetEnv("POLRED_SPEC"));

Qx<x> := PolynomialRing(Rationals());
possibility := "oneof{CC|RR+RR+RR+RR}";

// Treating a multi-factor alternative as one opaque token makes both valid
// embeddings fail. The parser must splice the selected list into the product.
assert RealRepresentationEmbeds(["CC"], [possibility]);
assert RealRepresentationEmbeds(
    ["RR", "RR", "RR", "RR"], [possibility]);
assert not RealRepresentationEmbeds(["HH", "HH"], [possibility]);

// Exact spelling does not make an ambiguous answer sharp. This genus-4 curve
// genuinely returns the pre-existing type II/III disjunction at B = 8.
C := HyperellipticCurve(x^9 - 6*x^7 - x^5 - 6*x^3 + x);
bound, B, status := RealRepresentationBound(
    C, ["M_2(RR) or HH"] : Bmax := 8);
assert bound eq ["M_2(RR) or HH"];
assert B eq 8;
assert status eq "exhausted";
assert RealRepresentationEmbeds(["HH"], bound);

// Empty alternatives and factors are malformed. Magma's Split discards empty
// pieces, so each spelling must be rejected before its pieces are parsed.
for token in [
    "oneof{CC||RR}",
    "oneof{|CC|RR}",
    "oneof{CC|RR|}",
    "oneof{CC+|RR}",
    "oneof{CC++CC|RR}",
    "oneof{CC|+RR}"
] do
    rejected := false;
    try
        result := RealRepresentationEmbeds(["CC"], [token]);
    catch e
        rejected := true;
    end try;
    assert rejected;
end for;
