#!/usr/bin/env python3
"""Normalize a corpus source file into the TSV the Magma driver reads.

Output columns (tab separated, one header line):

    id       unique key for the curve
    kind     'g2' (hyperelliptic minimal model) or 'g3' (plane quartic)
    data     the model, verbatim, in the form run_corpus.m knows how to eval
    expected comma-joined sorted multiset of RR-representation strings
    meta     extra stratum label carried for the report (g2: geom_end_alg)

Sources
-------
g2  LMFDB g2c_curves CSV with columns
      label,eqn,real_geom_end_alg,geom_end_alg,is_simple_geom,cond
    'eqn' is "[[f coeffs],[h coeffs]]", low-to-high, for y^2 + h(x) y = f(x).
    'real_geom_end_alg' is LMFDB text and is translated to the intrinsic's
    convention by REAL_GEOM_END_ALG below.

g3  one line per curve, 'id:[quartic in x,y,z]:[python list of RR strings]'.
    The third field is already in the intrinsic's convention. Ids in that file
    are NOT unique (they are a discriminant-like key, 425 of them repeat), so
    every id gets a '.k' occurrence suffix, k counting from 1 in file order.
    Raw ids are all digits, so the suffixed ids cannot collide.

Usage
-----
    prepare_input.py SOURCE {g2,g3} [-o OUT.tsv]
"""

import argparse
import ast
import collections
import csv
import sys

# LMFDB real_geom_end_alg -> the strings RealRepresentationBound emits.
# These six values are exhaustive for g2c_curves.
REAL_GEOM_END_ALG = {
    "R": ["RR"],
    "R x R": ["RR", "RR"],
    "C x R": ["CC", "RR"],
    "C x C": ["CC", "CC"],
    "M_2(R)": ["M_2(RR)"],
    "M_2(C)": ["M_2(CC)"],
}

HEADER = ("id", "kind", "data", "expected", "meta")


def clean(s):
    """Collapse anything that would corrupt the TSV framing."""
    return " ".join(str(s).split())


def load_g2(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = {"label", "eqn", "real_geom_end_alg"} - set(reader.fieldnames or ())
        if missing:
            sys.exit("%s: missing column(s) %s" % (path, ", ".join(sorted(missing))))
        for lineno, row in enumerate(reader, start=2):
            text = row["real_geom_end_alg"].strip()
            if text not in REAL_GEOM_END_ALG:
                sys.exit("%s:%d: unknown real_geom_end_alg %r" % (path, lineno, text))
            yield (
                clean(row["label"]),
                "g2",
                clean(row["eqn"]),
                ",".join(sorted(REAL_GEOM_END_ALG[text])),
                clean(row.get("geom_end_alg", "")),
            )


def load_g3(path):
    seen = collections.Counter()
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            head = line.find(":")
            tail = line.rfind(":")
            if head < 0 or tail <= head:
                sys.exit("%s:%d: expected 'id:[quartic]:[strings]'" % (path, lineno))
            raw_id, body, rr_text = line[:head], line[head + 1:tail], line[tail + 1:]
            body = body.strip()
            if not (body.startswith("[") and body.endswith("]")):
                sys.exit("%s:%d: quartic field is not bracketed" % (path, lineno))
            quartic = body[1:-1].strip()
            try:
                rr = ast.literal_eval(rr_text.strip())
            except (SyntaxError, ValueError):
                sys.exit("%s:%d: cannot parse %r" % (path, lineno, rr_text))
            if not isinstance(rr, (list, tuple)) or not all(isinstance(s, str) for s in rr):
                sys.exit("%s:%d: third field is not a list of strings" % (path, lineno))
            raw_id = clean(raw_id)
            seen[raw_id] += 1
            yield (
                "%s.%d" % (raw_id, seen[raw_id]),
                "g3",
                clean(quartic),
                ",".join(sorted(rr)),
                "",
            )


LOADERS = {"g2": load_g2, "g3": load_g3}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("source", help="path to the corpus source file")
    ap.add_argument("kind", choices=sorted(LOADERS), help="how to parse SOURCE")
    ap.add_argument("-o", "--output", help="destination TSV (default: stdout)")
    args = ap.parse_args(argv)

    out = open(args.output, "w") if args.output else sys.stdout
    ids = set()
    n = 0
    try:
        out.write("\t".join(HEADER) + "\n")
        for row in LOADERS[args.kind](args.source):
            if row[0] in ids:
                sys.exit("duplicate id %r; ids must be unique" % (row[0],))
            ids.add(row[0])
            out.write("\t".join(row) + "\n")
            n += 1
    finally:
        if args.output:
            out.close()
    print("prepare_input: wrote %d rows (%s)" % (n, args.kind), file=sys.stderr)


if __name__ == "__main__":
    main()
