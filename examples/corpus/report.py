#!/usr/bin/env python3
"""Aggregate run_corpus.sh result TSVs into a per-stratum sharpness table.

The verdict is computed here rather than in Magma so that the comparison does
not depend on Magma and Python agreeing on how to sort strings: both sides are
re-sorted before they are compared.

    sharp     the multisets are equal
    over      got is a genuine upper bound: expected embeds unitally in got,
              but the two are not equal
    under     dim(got) < dim(expected)
    mismatch  equal or larger dimension but structurally incomparable
    error     the Magma driver raised on this curve

'under' is expected in one place. For geometrically simple genus-2 curves whose
quartic CM is only defined over an extension of Q, the extra endomorphisms are
invisible to prime-field Frobenius, so the bound returns [RR, RR] where the
truth is [CC, CC]. eta is a min over primes, so raising B does not help. That
is a documented limitation of the method, not a harness bug, which is why it is
counted rather than treated as failure.

When several passes have computed the same id at different B, the rows are
folded to one. The reported algebra is NOT monotone in B (see README.md), so
the fold keeps the row with the smallest total dimension of 'got' rather than
the largest B; ties go to the larger B. --all-bounds turns the fold off.

Usage:
    report.py RESULTS.tsv [RESULTS.tsv ...] [--meta PREPARED.tsv ...]
                          [--list-nonsharp IDS.txt] [--by-id] [--all-bounds]
"""

import argparse
import collections
import re
import sys

BASE_DIMENSIONS = {"RR": 1, "CC": 2, "HH": 4}
MATRIX_TOKEN = re.compile(r"M_([1-9][0-9]*)\((RR|CC|HH)\)\Z")

VERDICTS = ("sharp", "over", "under", "mismatch", "error")

RESULT_COLUMNS = 7  # id, kind, B, expected, got, cputime, status


def parse_multiset(text):
    return sorted(part for part in (p.strip() for p in text.split(",")) if part)


def _factor_data(tok):
    """A simple real algebra token, as (matrix degree, base dimension)."""
    if tok in BASE_DIMENSIONS:
        return 1, BASE_DIMENSIONS[tok]
    match = MATRIX_TOKEN.fullmatch(tok)
    if match is None:
        return None
    return int(match.group(1)), BASE_DIMENSIONS[match.group(2)]


def _token_candidates(tok):
    """Alternative token lists represented by one possibly ambiguous token."""
    if tok.startswith("oneof{") and tok.endswith("}"):
        body = tok[6:-1]
        alternatives = body.split("|")
        if len(alternatives) < 2 or any(not alternative for alternative in alternatives):
            return ()
        candidates = tuple(tuple(alternative.split("+"))
                           for alternative in alternatives)
        if any(any(not token for token in alternative)
               for alternative in candidates):
            return ()
        return candidates
    return ((tok,),)


def _is_ambiguous(tok):
    return _token_candidates(tok) != ((tok,),)


def _unknown_tokens(items):
    unknown = []
    for tok in items:
        candidates = _token_candidates(tok)
        if not candidates:
            unknown.append(tok)
        elif candidates == ((tok,),):
            if _factor_data(tok) is None:
                unknown.append(tok)
        else:
            for alternative in candidates:
                unknown.extend(_unknown_tokens(alternative))
    return unknown


def total_dim(items):
    """Maximum total real dimension over resolutions, plus unknown tokens."""
    unknown = _unknown_tokens(items)
    if unknown:
        return 0, unknown
    return max(sum(k * k * d for k, d in map(_factor_data, resolution))
               for resolution in _resolutions(items)), []


def _support_masks(units, target):
    """Supports, as bitmasks over i, of the solutions in non-negative integers
    of sum_i c_i * units[i] == target, excluding the all-zero c.
    """
    found = set()

    def walk(i, rest, mask):
        if i == len(units):
            if rest == 0 and mask:
                found.add(mask)
            return
        for c in range(rest // units[i] + 1):
            walk(i + 1, rest - c * units[i],
                 (mask | (1 << i)) if c else mask)

    walk(0, target, 0)
    return found


def _embeds_resolved(exp, obs):
    """Is there a unital injective R-algebra map prod(exp) -> prod(obs)?

    Per factor M_n(E) of obs: sum_i c_i k_i max(dim D_i, dim E) = n dim E,
    then no factor of obs left unhit and no factor of exp left unused.
    """
    shape = [_factor_data(tok) for tok in exp]
    reach = set([0])
    for tok in obs:
        n, e = _factor_data(tok)
        masks = _support_masks([k * max(d, e) for k, d in shape], n * e)
        if not masks:
            return False
        reach = set(r | m for r in reach for m in masks)
    return (1 << len(exp)) - 1 in reach


def _resolutions(items):
    """Every way of splicing each ambiguous token's alternative list."""
    if not items:
        return [[]]
    tok = items[0]
    tail = items[1:]
    candidates = _token_candidates(tok)
    if candidates == ((tok,),):
        return [[tok] + rest for rest in _resolutions(tail)]
    return [head + rest
            for alternative in candidates
            for head in _resolutions(alternative)
            for rest in _resolutions(tail)]


_EMBEDS = {}


def embeds(exp, obs):
    """Does prod(exp) embed unitally in prod(obs), for some resolution of the
    ambiguous tokens? Both arguments are sorted token lists.
    """
    key = (",".join(exp), ",".join(obs))
    hit = _EMBEDS.get(key)
    if hit is None:
        hit = any(_embeds_resolved(e, o)
                  for e in _resolutions(exp) for o in _resolutions(obs))
        _EMBEDS[key] = hit
    return hit


def classify(expected, got, status):
    if status != "ok":
        return "error", []
    exp = parse_multiset(expected)
    obs = parse_multiset(got)
    if exp == obs and not any(_is_ambiguous(tok) for tok in obs):
        return "sharp", []
    de, unk_e = total_dim(exp)
    dg, unk_g = total_dim(obs)
    unknown = unk_e + unk_g
    if unknown:
        return "mismatch", unknown
    if embeds(exp, obs):
        return "over", []
    if dg < de:
        return "under", []
    return "mismatch", []


def tightness(row):
    """Fold key: smaller is a better bound. Sorts (rank, dim, -B).

    Rank first, because two of the three outcomes carry no bound to compare:
    an empty 'got' (RealRepresentationBound returned []) and an 'error' row
    both have dimension 0 and would otherwise win every fold. Only within
    rank 0 does the dimension mean anything, and there the smallest total
    dimension is the tightest bound, whatever B produced it.
    """
    obs = parse_multiset(row["got"])
    dim, unknown = total_dim(obs)
    if row["status"] != "ok" or unknown:
        rank = 2
    elif not obs:
        rank = 1
    else:
        rank = 0
    return (rank, dim, -int(row["B"]))


def load_results(paths):
    for path in paths:
        with open(path) as fh:
            for lineno, line in enumerate(fh, start=1):
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if fields[0] == "id" and len(fields) > 1 and fields[1] == "kind":
                    continue
                if len(fields) < RESULT_COLUMNS:
                    print("%s:%d: skipping short line" % (path, lineno), file=sys.stderr)
                    continue
                cid, kind, bnd, expected, got, cpu, status = fields[:RESULT_COLUMNS]
                try:
                    cpu = float(cpu)
                except ValueError:
                    cpu = 0.0
                yield dict(id=cid, kind=kind, B=bnd, expected=expected,
                           got=got, cputime=cpu, status=status)


def load_meta(paths):
    meta = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 5 and fields[0] != "id":
                    meta[fields[0]] = fields[4]
    return meta


def render(title, rows, out):
    """rows: list of (label, counter, n). Prints counts with percentages."""
    label_w = max([len("stratum")] + [len(r[0]) for r in rows])
    cols = ["n"] + list(VERDICTS)
    widths = {"n": max(6, len("n"))}
    for v in VERDICTS:
        widths[v] = max(len(v), 13)

    header = "stratum".ljust(label_w) + "  " + "n".rjust(widths["n"])
    for v in VERDICTS:
        header += "  " + v.rjust(widths[v])
    print(title, file=out)
    print("=" * len(header), file=out)
    print(header, file=out)
    print("-" * len(header), file=out)
    for label, counter, n in rows:
        line = label.ljust(label_w) + "  " + str(n).rjust(widths["n"])
        for v in VERDICTS:
            c = counter[v]
            cell = "-" if c == 0 else "%d (%.1f%%)" % (c, 100.0 * c / n)
            line += "  " + cell.rjust(widths[v])
        print(line, file=out)
    print("", file=out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("results", nargs="+", help="result TSVs from run_corpus.sh")
    ap.add_argument("--meta", nargs="*", default=[],
                    help="prepared TSVs, to recover the extra stratum column "
                         "(genus 2: geom_end_alg)")
    ap.add_argument("--list-nonsharp", metavar="FILE",
                    help="write the ids worth re-running at a higher B")
    ap.add_argument("--nonsharp-verdicts", default="over,mismatch,error",
                    help="which verdicts --list-nonsharp collects "
                         "(default: over,mismatch,error; 'under' is excluded "
                         "because a larger B cannot recover a missed factor)")
    ap.add_argument("--by-id", action="store_true",
                    help="also dump every non-sharp row")
    ap.add_argument("--all-bounds", action="store_true",
                    help="report every (id, B) row separately instead of "
                         "folding each id to its tightest row")
    args = ap.parse_args(argv)

    meta = load_meta(args.meta)
    out = sys.stdout

    by_kind = collections.defaultdict(collections.Counter)
    by_stratum = collections.defaultdict(collections.Counter)
    stratum_order = []
    overall = collections.Counter()
    cputime = collections.defaultdict(float)
    duplicates = 0
    superseded = 0
    nobound = 0
    unknown_tokens = collections.Counter()
    bounds = collections.Counter()
    nonsharp_ids = []
    nonsharp_rows = []
    wanted = set(v.strip() for v in args.nonsharp_verdicts.split(",") if v.strip())

    # The reported algebra is NOT monotone in B (README.md, "The bound is not
    # monotone in B"), so a larger B is not a better answer. Fold on the
    # tightest bound instead: smallest total dimension of 'got', ties to the
    # larger B. Since only non-sharp ids are ever re-run, and 'under' ids are
    # not re-run at all, this can lose a sharp row to a smaller-dimensional
    # 'under' row but can never turn a non-sharp id into a false 'sharp'.
    best = {}
    largest_b = {}
    for row in load_results(args.results):
        key = (row["id"], row["kind"], row["B"]) if args.all_bounds \
            else (row["id"], row["kind"])
        previous = best.get(key)
        if previous is None:
            best[key] = row
        elif args.all_bounds or row["B"] == previous["B"]:
            duplicates += 1
        else:
            superseded += 1
            if tightness(row) < tightness(previous):
                best[key] = row
        if not args.all_bounds:
            largest_b[key] = max(largest_b.get(key, 0), int(row["B"]))

    kept_smaller_b = sum(1 for key, row in best.items()
                         if int(row["B"]) < largest_b.get(key, 0))

    for row in best.values():
        verdict, unknown = classify(row["expected"], row["got"], row["status"])
        for tok in unknown:
            unknown_tokens[tok] += 1
        if row["status"] == "ok" and not parse_multiset(row["got"]):
            nobound += 1

        kind = row["kind"]
        expected = ",".join(parse_multiset(row["expected"]))
        extra = meta.get(row["id"], "")
        label = "%s  %s" % (kind, expected)
        if extra:
            label += "  [%s]" % extra

        if label not in by_stratum:
            stratum_order.append(label)
        by_stratum[label][verdict] += 1
        by_kind[kind][verdict] += 1
        overall[verdict] += 1
        cputime[kind] += row["cputime"]
        bounds[(kind, row["B"])] += 1

        if verdict != "sharp":
            if verdict in wanted:
                nonsharp_ids.append(row["id"])
            nonsharp_rows.append((row["id"], kind, verdict, expected,
                                  ",".join(parse_multiset(row["got"])) or "(none)",
                                  row["status"]))

    n_total = sum(overall.values())
    if n_total == 0:
        print("report: no result rows found", file=sys.stderr)
        return 1

    rows = []
    for label in sorted(stratum_order):
        counter = by_stratum[label]
        rows.append((label, counter, sum(counter.values())))
    render("Sharpness by expected multiset", rows, out)

    rows = []
    for kind in sorted(by_kind):
        counter = by_kind[kind]
        rows.append((kind, counter, sum(counter.values())))
    rows.append(("ALL", overall, n_total))
    render("Totals", rows, out)

    print("Bounds used:", ", ".join(
        "%s B=%s: %d" % (k, b, n) for (k, b), n in sorted(bounds.items())), file=out)
    for kind in sorted(cputime):
        n = sum(by_kind[kind].values())
        print("Magma cputime, %s: %.0f s total, %.3f s/curve over %d curves"
              % (kind, cputime[kind], cputime[kind] / n, n), file=out)
    if nobound:
        print("Note: %d row(s) established no bound at all "
              "(RealRepresentationBound returned []); they are counted as 'under'."
              % nobound, file=out)
    if superseded:
        print("Note: %d row(s) folded away; one row kept per id, the one with "
              "the smallest dim(got)." % superseded, file=out)
    if kept_smaller_b:
        print("Note: for %d id(s) the kept row is not the largest B computed "
              "(the bound is not monotone in B)." % kept_smaller_b, file=out)
    if duplicates:
        print("Note: %d duplicate (id, kind, B) row(s) ignored." % duplicates, file=out)
    if unknown_tokens:
        print("Note: unrecognised algebra string(s), counted as 'mismatch': %s"
              % ", ".join("%s x%d" % kv for kv in unknown_tokens.most_common()), file=out)

    if args.by_id:
        print("", file=out)
        print("Non-sharp rows", file=out)
        for r in sorted(nonsharp_rows, key=lambda r: (r[2], r[1], r[0])):
            print("  %-16s %-3s %-9s expected=%-18s got=%-18s %s" % r, file=out)

    if args.list_nonsharp:
        with open(args.list_nonsharp, "w") as fh:
            for cid in nonsharp_ids:
                fh.write(cid + "\n")
        print("", file=out)
        print("Wrote %d id(s) [%s] to %s"
              % (len(nonsharp_ids), args.nonsharp_verdicts, args.list_nonsharp),
              file=out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
