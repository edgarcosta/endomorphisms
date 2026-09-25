#!/usr/bin/env python3
"""Regression tests for report.classify, run as: python3 test_report.py

Silent on success, exit 0. On failure it prints each bad case and exits 1.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import report

# (expected, got, status, verdict, unknown tokens).
CASES = [
    # The bug: 4 > 3, but C x R has no unital injective map into M_2(R). The
    # two central idempotents would be complementary and rank 1, so the C
    # summand would have to land in P M_2(R) P = R. Units (2, 1) cannot sum
    # to 2 with both c_i positive.
    ("CC,RR", "M_2(RR)", "ok", "mismatch", []),

    # The five non-sharp shapes present in the g2 corpus; all stay 'over'.
    ("CC,RR", "CC,CC", "ok", "over", []),       # C -> C, R -> C
    ("CC,RR", "RR,M_2(RR)", "ok", "over", []),  # R -> R, C -> M_2(R)
    ("RR,RR", "CC,CC", "ok", "over", []),       # one R into each C
    ("CC,RR", "CC,M_2(RR)", "ok", "over", []),  # C -> C, R^2 -> M_2(R)
    ("RR,RR", "M_2(RR)", "ok", "over", []),     # R x R diagonal in M_2(R)

    # sharp: the multisets are compared after sorting and stripping.
    ("CC,RR", "CC,RR", "ok", "sharp", []),
    ("RR, CC", "CC,RR", "ok", "sharp", []),
    ("M_2(CC)", "M_2(CC)", "ok", "sharp", []),
    ("M_4(RR)", "M_4(RR)", "ok", "sharp", []),
    ("M_2(HH)", "M_2(HH)", "ok", "sharp", []),
    ("M_3(HH)", "M_3(HH)", "ok", "sharp", []),

    # under: the documented limitation, quartic CM invisible to prime-field
    # Frobenius. Units are 2 and 2, the targets 1 and 1, so nothing embeds.
    ("CC,CC", "RR,RR", "ok", "under", []),
    ("CC,RR", "", "ok", "under", []),

    # error wins over the algebras, whatever they are.
    ("CC,RR", "CC,RR", "fail", "error", []),
    ("RR", "M_2(RR)", "magma error", "error", []),

    # Unknown tokens keep the old behaviour on either side.
    ("RR", "M_5(RR)", "ok", "mismatch", ["M_5(RR)"]),
    ("QQ", "M_2(RR)", "ok", "mismatch", ["QQ"]),

    # u = 1, so c = 2 fills dim_R R^2 = 2.
    ("RR", "M_2(RR)", "ok", "over", []),
    # u = 1 max(2, 1) = 2 = dim_R R^2, c = 1.
    ("CC", "M_2(RR)", "ok", "over", []),
    # Three units of 1, each c_i positive, cannot sum to 2; dim 4 > 3.
    ("RR,RR,RR", "M_2(RR)", "ok", "mismatch", []),
    # u = 2 max(1, 2) = 4 = dim_R C^2, c = 1: M_2(R) -> M_2(C) by inclusion.
    ("M_2(RR)", "M_2(CC)", "ok", "over", []),
    # u = 1 max(4, 1) = 4 > 2: H needs M_4(R). Equal dimension, so not under.
    ("HH", "M_2(RR)", "ok", "mismatch", []),
    # Units 2 and 2, target dim_R C^2 = 4, c = (1, 1): the diagonal.
    ("CC,CC", "M_2(CC)", "ok", "over", []),

    # The ambiguous token, on the got side: only the M_2(RR) resolution works.
    ("RR,RR", "M_2(RR) or HH", "ok", "over", []),
    # Neither resolution: units (4, 4) into H, units (2, 1) into M_2(R).
    ("CC,RR", "M_2(RR) or HH", "ok", "mismatch", []),
    # On the expected side: M_2(R) and H both embed in M_2(C), u = 4 = target.
    ("M_2(RR) or HH", "M_2(CC)", "ok", "over", []),

    # A centre possibility token may resolve to a list of several factors.
    ("CC", "oneof{CC|RR+RR+RR+RR}", "ok", "over", []),
    ("RR,RR,RR,RR", "oneof{CC|RR+RR+RR+RR}", "ok", "over", []),
    # Even literal equality with an ambiguous bound is not evidence of sharpness.
    ("oneof{CC|RR+RR+RR+RR}", "oneof{CC|RR+RR+RR+RR}",
     "ok", "over", []),

    # Empty alternatives and factors are malformed, even when Split-like
    # parsers would discard the empty piece.
    ("CC", "oneof{CC||RR}", "ok", "mismatch", ["oneof{CC||RR}"]),
    ("CC", "oneof{|CC|RR}", "ok", "mismatch", ["oneof{|CC|RR}"]),
    ("CC", "oneof{CC|RR|}", "ok", "mismatch", ["oneof{CC|RR|}"]),
    ("CC", "oneof{CC+|RR}", "ok", "mismatch", ["oneof{CC+|RR}"]),
    ("CC", "oneof{CC++CC|RR}", "ok", "mismatch", ["oneof{CC++CC|RR}"]),
    ("CC", "oneof{CC|+RR}", "ok", "mismatch", ["oneof{CC|+RR}"]),
]


def main():
    failures = []
    for expected, got, status, verdict, unknown in CASES:
        result = report.classify(expected, got, status)
        if result != (verdict, unknown):
            failures.append("classify(%r, %r, %r) -> %r, want %r"
                            % (expected, got, status, result,
                               (verdict, unknown)))
    dim_result = report.total_dim(["oneof{CC|RR+RR+RR+RR}"])
    if dim_result != (4, []):
        failures.append("total_dim(disjunction) -> %r, want (4, [])"
                        % (dim_result,))
    for line in failures:
        print(line)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
