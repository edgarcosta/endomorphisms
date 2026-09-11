#!/usr/bin/env python3
"""Cut a reproducible stratified sample out of a prepared corpus TSV.

The stratum of a row is the pair (expected, meta), columns 4 and 5 of the
prepared TSV, which is exactly what report.py groups by.

Reproducibility is the point of this script. Selection does not use
random.sample: which rows a PRNG hands back depends on the order the pools were
built in and on the CPython version's sampling algorithm, so the same seed does
not survive a Python upgrade or a re-ordered input. Instead each row gets the
digest

    blake2b(seed || "\\0" || id)

and a stratum keeps the rows with the smallest digests. That depends on nothing
but the seed and the id, so the same (seed, spec, input) always yields the same
id set, and adding curves to the input only ever adds candidates rather than
reshuffling the whole sample. Rows are written in input order.

Usage
-----
    sample_corpus.py PREPARED.tsv -o SAMPLE.tsv [--seed S]
                     [--default N|all] [--take 'EXPECTED|META=N|all'] ...
    sample_corpus.py PREPARED.tsv --list-strata

A stratum spec is 'EXPECTED|META=COUNT', where EXPECTED and META are the two
column values verbatim (META is empty for g3, so the key ends in '|'), and
COUNT is a row count or the word 'all'. Strata not named take --default; a
stratum smaller than its count is taken whole. Example, the standard genus-2
sample:

    sample_corpus.py prepared/g2.tsv -o prepared/g2_sample.tsv --seed 0 \\
        --default 200 --take 'M_2(RR)|QM=all' \\
        --take 'CC,CC|CM=all' --take 'CC,CC|CM x CM=all'
"""

import argparse
import collections
import hashlib
import sys

HEADER_FIELDS = ("id", "kind", "data", "expected", "meta")
ALL = "all"


def stratum_key(fields):
    expected = fields[3] if len(fields) > 3 else ""
    meta = fields[4] if len(fields) > 4 else ""
    return "%s|%s" % (expected, meta)


def digest(seed, cid):
    return hashlib.blake2b(("%s\0%s" % (seed, cid)).encode("utf-8"),
                           digest_size=16).digest()


def parse_take(text):
    if "=" not in text:
        sys.exit("--take wants 'EXPECTED|META=COUNT', got %r" % text)
    key, count = text.rsplit("=", 1)
    return key, parse_count(count, "--take " + text)


def parse_count(text, where):
    text = text.strip()
    if text.lower() == ALL:
        return None
    try:
        n = int(text)
    except ValueError:
        sys.exit("%s: count must be an integer or 'all', got %r" % (where, text))
    if n < 0:
        sys.exit("%s: count must not be negative" % where)
    return n


def load(path):
    """Yield (index, key, id, line) for every data row, plus the header line."""
    header = None
    rows = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if lineno == 1 and fields[0] == "id":
                header = line
                continue
            if len(fields) < 4:
                sys.exit("%s:%d: fewer than 4 columns" % (path, lineno))
            rows.append((len(rows), stratum_key(fields), fields[0], line))
    if header is None:
        header = "\t".join(HEADER_FIELDS) + "\n"
    return header, rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("source", help="prepared TSV from prepare_input.py")
    ap.add_argument("-o", "--output", help="destination TSV (default: stdout)")
    ap.add_argument("--seed", default="0",
                    help="selection seed, any string (default: 0)")
    ap.add_argument("--default", dest="default_count", default="200",
                    metavar="N", help="rows per unnamed stratum, or 'all' "
                                      "(default: 200)")
    ap.add_argument("--take", action="append", default=[], metavar="SPEC",
                    help="'EXPECTED|META=N' or 'EXPECTED|META=all', repeatable")
    ap.add_argument("--list-strata", action="store_true",
                    help="print the strata and their sizes, then exit")
    args = ap.parse_args(argv)

    header, rows = load(args.source)

    pools = collections.OrderedDict()
    for row in rows:
        pools.setdefault(row[1], []).append(row)

    if args.list_strata:
        for key in sorted(pools):
            print("%8d  %s" % (len(pools[key]), key))
        print("%8d  TOTAL (%d strata)" % (len(rows), len(pools)))
        return 0

    default_count = parse_count(args.default_count, "--default")
    take = {}
    for spec in args.take:
        key, count = parse_take(spec)
        take[key] = count
    unknown = sorted(set(take) - set(pools))
    if unknown:
        sys.exit("--take names stratum/strata not present in %s: %s\n"
                 "(run --list-strata to see the keys)"
                 % (args.source, ", ".join(repr(k) for k in unknown)))

    chosen = set()
    summary = []
    for key in sorted(pools):
        pool = pools[key]
        count = take.get(key, default_count)
        if count is None or count >= len(pool):
            picked = pool
        else:
            picked = sorted(pool, key=lambda r: (digest(args.seed, r[2]), r[2]))[:count]
        chosen.update(r[0] for r in picked)
        summary.append((key, len(pool), len(picked)))

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        out.write(header)
        for row in rows:
            if row[0] in chosen:
                out.write(row[3])
    finally:
        if args.output:
            out.close()

    for key, have, took in summary:
        print("sample_corpus: %6d of %6d  %s" % (took, have, key), file=sys.stderr)
    print("sample_corpus: %d rows (seed %r) -> %s"
          % (len(chosen), args.seed, args.output or "stdout"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
