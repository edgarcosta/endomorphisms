# Corpus harness for `RealRepresentationBound`

Thanks for picking this up. The harness answers one question at scale: how often
does the upper-bound search

```
RealRepresentationBound(C::Crv, target::SeqEnum : Bmax := 1024)
```

return the known geometric endomorphism data sharply? It runs a prepared list
of curves through Magma in parallel, one process per chunk, and passes the
reference multiset as the target. `-B` sets `Bmax`, the cap on the search.

Five pieces, in the order you use them:

| file | what it does |
| --- | --- |
| `prepare_input.py` | turns a corpus source file into the driver's TSV |
| `sample_corpus.py` | cuts a reproducible stratified sample out of that TSV |
| `run_corpus.m` | the Magma driver: one chunk in, one result line per curve out |
| `run_corpus.sh` | shards, parallelizes, resumes, concatenates |
| `report.py` | aggregates results into a per-stratum sharpness table |

## PARI/gp is mandatory

This is the one thing that will silently ruin a run, so it comes first.

`FieldIntersectionMatrix` computes each center by **set intersection on
polredabs-canonicalized defining polynomials**. Polredabs is what makes field
equality testable with `eq` at all, and `Polredabs` shells out to PARI/gp. When
gp is missing, `Polred` catches the failed `Pipe`, prints a warning, and returns
a non-canonical fallback polynomial. Nothing raises. The meet then loses
elements, centers collapse to Q, and results degrade quietly: `["RR","CC"]`
comes back as `["RR","RR"]`.

So both `run_corpus.m` and `run_corpus.sh` refuse to start unless

```
Polredabs(x^2 - 5) eq x^2 - x - 1
```

Round-trip, not `command -v gp`: a gp on PATH that Magma cannot actually pipe to
would pass the second test and fail the first. A failed guard exits **2** and
writes no results.

Any recent PARI/gp works. If gp is not on the default PATH, Sage ships one:

```sh
export PATH="$(dirname "$(sage -sh -c 'command -v gp')"):$PATH"
```

Nothing in the harness hardcodes that location; the caller prepends it.

## Configuration

Two `AttachSpec` paths, both environment variables:

| variable | default | meaning |
| --- | --- | --- |
| `ENDO_SPEC` | `<repo>/endomorphisms/magma/spec`, derived from `run_corpus.sh`'s own location | this repository's Magma spec |
| `POLRED_SPEC` | unset | `CHIMP/MagmaPolred/spec`; leave unset if `MAGMA_USER_SPEC` already attaches it |
| `MAGMA` | `magma` | the Magma binary |

On the machine this was developed on, `POLRED_SPEC` is
`.../CHIMP/CHIMP/MagmaPolred/spec`. Set it to wherever your CHIMP checkout is.

## Running it

```sh
export POLRED_SPEC=/path/to/CHIMP/MagmaPolred/spec
export PATH="$(dirname "$(sage -sh -c 'command -v gp')"):$PATH"

# 1. normalize the source data
./prepare_input.py g2_nongeneric.csv g2 -o prepared/g2.tsv
./prepare_input.py g3_quartics.txt    g3 -o prepared/g3.tsv

# 1b. optionally cut a stratified sample instead of running everything
#     (see "Stratified sample" below for the exact spec)
./sample_corpus.py prepared/g2.tsv -o prepared/g2_sample.tsv --seed 0

# 2. run it
./run_corpus.sh -i prepared/g2.tsv -o results/g2.tsv -B 1024 -j 200 -c 50

# 3. report
./report.py results/g2.tsv --meta prepared/g2.tsv
```

`run_corpus.sh --help` lists the flags: `-B` climb cap, `-j` jobs, `-c` chunk
size, `-t` per-chunk timeout, `-w` work directory, `--ids` to restrict a pass
to a list of ids. The default cap is 1024, which runs the full search.

Sizing `-c`: each chunk pays one Magma startup plus one gp round-trip, about a
second. At 0.35 s/curve for genus 2 a chunk of 200-400 keeps that under 1%, and
smaller chunks give the scheduler more to balance and lose less work to a
timeout. On a 200-core box, `-j 200 -c 400` is a reasonable starting point.

### Resuming

`run_corpus.m` writes each result with `PrintFile`, which opens and closes the
file per call and therefore flushes. (Printing would not work: Magma
block-buffers redirected stdout, so a killed process would lose everything it
had computed.) Each invocation of `run_corpus.sh`

1. folds any chunk output left in the work directory into the results file,
   keeping one row per `(id, kind, B)`;
2. computes the pending set as the input ids that have no row in the results
   file *at this `B`*;
3. runs only those.

So `Ctrl-C`, `kill -9`, a chunk hitting `timeout`, or a machine reboot all cost
at most the curves that were mid-flight. Re-run the identical command line to
continue. The work directory (`<output>.work` by default) is removed only when
the pass completes with nothing missing.

Because it is removed wholesale, the results file must not live inside it.
`-w run -o run/results.tsv` would have the completing run delete its own
output, so the script refuses that combination up front with exit **64**.

Resume is currently unreliable for these escalating runs: result column 3 is
the rung reached, while the planner compares it with the requested cap. To
restart or extend a run, supply `--ids` explicitly. Separate result files per
pass are safe, and `report.py` accepts any number of them.

### The bound is not monotone in `B`

Read this before designing a multi-pass run: raising `B` can make the reported
algebra *bigger*.

eta itself is monotone. It is a `min` over primes, so more primes can only
lower it. The algebra is not, and the reason is that lowering eta throws away
the evidence collected so far. `EndomorphismAlgebraEtaBound` keeps
`eta_lower`, the per-prime endomorphism factorizations from the primes that
attain the current minimum, and resets it to `[]` the moment a smaller eta
turns up (`endomorphisms/magma/upperbounds/UpperBounds.m:68-72`, and the same
in the Sage original at `endomorphisms/UpperBounds/upper_bounds.py:65-73`,
which resets on both the smaller-eta and the fewer-factors branch).
`FieldIntersectionMatrix` then intersects only over the primes still in
`eta_lower`. A larger `B` that finds a new eta-minimizing prime can therefore
leave *fewer* primes in the intersection, and an intersection over fewer primes
is larger, so the centers are larger and so is the total dimension.

This is not a bug in the Magma port. The Magma is a faithful translation of the
Sage, which does the same thing; do not "fix" the intrinsic on the strength of
this note.

Three genus-3 curves in this corpus show it directly, all three with the true
answer `RR,RR,RR` (dimension 3):

| id | `B=64` | `B=256` | `B=512` |
| --- | --- | --- | --- |
| `802816.1` | `RR,M_2(RR)` (5) | `CC,CC,CC` (6) | `RR,RR,RR` (3) |
| `802816.2` | `RR,M_2(RR)` (5) | `CC,CC,CC` (6) | `RR,RR,RR` (3) |
| `3781323.1` | `RR,M_2(RR)` (5) | `CC,CC,CC` (6) | `RR,RR,RR` (3) |

So "the largest `B` you ran" is not a synonym for "the tightest bound you
have", and `report.py` does not treat it as one; see below.

### Two-pass strategy

The multi-pass strategy survives the previous section intact, for two reasons:
only ids that came back non-sharp are re-run, and the fold keeps the *tightest*
row rather than the newest one. A re-run can therefore waste `B` on a curve, but
it can never make the combined table worse than the cheap pass alone, and it can
never manufacture a false `sharp`: a row is only reported `sharp` if some pass
actually computed that multiset. Run the whole corpus cheaply, then spend the
expensive `B` only where it can help:

```sh
./run_corpus.sh -i prepared/g2.tsv -o results/g2.B64.tsv -B 64 -j 200 -c 400
./report.py results/g2.B64.tsv --meta prepared/g2.tsv --list-nonsharp retry.txt

./run_corpus.sh -i prepared/g2.tsv -o results/g2.B1024.tsv -B 1024 -j 200 -c 50 \
                --ids retry.txt
./report.py results/g2.B64.tsv results/g2.B1024.tsv --meta prepared/g2.tsv
```

`--list-nonsharp` collects the `over`, `mismatch` and `error` ids by default and
deliberately leaves out `under`. An `under` is the prime-field Frobenius
limitation documented below, not a `B` that was too small; and since the fold
keeps the smallest-dimension row, a larger-`B` row could not displace an
`under` one anyway. Cost goes roughly as `B`: measured at genus 2, 0.35 s/curve
at `B=50`, 0.73 s at `B=100`, 2.0 s at `B=200`.

Hand `report.py` both files and it folds each id to one row: the one with the
**smallest total dimension** of `got`, ties going to the larger `B`. That is
the tightest bound any pass established, which is the state of knowledge after
both passes. Rows carrying no bound at all (empty `got`, or status `error`)
also have dimension 0 and are ranked last rather than first, so they never win
a fold against a real answer. `--all-bounds` turns the fold off and reports
every `(id, B)` separately, which is how you see what the extra `B` actually
bought. The report prints a note counting the ids whose kept row is not the
largest `B` computed, which is the non-monotonicity showing up in your own data.

## Input formats

`prepare_input.py` emits `id, kind, data, expected, meta`, tab separated, with a
header. `expected` is the comma-joined sorted multiset in the intrinsic's
convention. `meta` is an extra stratum label the report can group by.

**`g2`**: LMFDB CSV with columns
`label,eqn,real_geom_end_alg,geom_end_alg,is_simple_geom,cond`. `eqn` is
`"[[f coeffs],[h coeffs]]"`, low-to-high, for the minimal model
`y^2 + h(x) y = f(x)`; the driver builds `HyperellipticCurve(R!f, R!h)`.
`real_geom_end_alg` is LMFDB text, translated as

| LMFDB | intrinsic |
| --- | --- |
| `R` | `RR` |
| `R x R` | `RR,RR` |
| `C x R` | `CC,RR` |
| `C x C` | `CC,CC` |
| `M_2(R)` | `M_2(RR)` |
| `M_2(C)` | `M_2(CC)` |

Those six values are exhaustive for the table; an unknown one aborts.
`meta` carries `geom_end_alg`, which is what separates `QM` from `M_2(Q)` and
geometrically simple `CM` from `CM x CM` in the report.

**`g3`**: one line per curve, `id:[quartic in x,y,z]:['M_2(RR)', 'RR']`. The
third field is already in the intrinsic's convention. The driver builds
`Curve(P2, CR ! eval("return " cat body cat ";"))` with `P2<x,y,z>`; Magma's
`eval` needs both the `return` and the generators in scope, which is why
`run_corpus.m` names the univariate ring's generator `t`. Ids in that file are
*not* unique (425 repeat), so every id gets a `.k` occurrence suffix counting
from 1 in file order. Raw ids are all digits, so the suffixed ids cannot
collide.

## Pulling the full genus-2 table

The staged corpus covers only the 199,627 non-generic curves. The full table is
`g2c_curves_new`, 6,216,959 rows; the 6,017,332 generic ones are
`real_geom_end_alg IN ('R','R x R') AND geom_end_alg <> 'RM'`.

Everything, for the big machine:

```sql
\copy (
    SELECT label, eqn, real_geom_end_alg, geom_end_alg, is_simple_geom, cond
    FROM g2c_curves_new
    ORDER BY cond, label
) TO 'g2_all.csv' WITH (FORMAT csv, HEADER true)
```

Just the non-generic part, which reproduces the staged `g2_nongeneric.csv`:

```sql
\copy (
    SELECT label, eqn, real_geom_end_alg, geom_end_alg, is_simple_geom, cond
    FROM g2c_curves_new
    WHERE NOT (real_geom_end_alg IN ('R', 'R x R') AND geom_end_alg <> 'RM')
    ORDER BY cond, label
) TO 'g2_nongeneric.csv' WITH (FORMAT csv, HEADER true)
```

Strata, for sizing a run:

| `real_geom_end_alg` | `geom_end_alg` | simple | rows |
| --- | --- | --- | ---: |
| `R x R` | `Q x Q` | no | 3,221,210 |
| `R` | `Q` | yes | 2,796,122 |
| `C x R` | `CM x Q` | no | 172,773 |
| `M_2(R)` | `M_2(Q)` | no | 19,993 |
| `R x R` | `RM` | yes | 4,990 |
| `M_2(C)` | `M_2(CM)` | no | 1,135 |
| `C x C` | `CM x CM` | no | 458 |
| `M_2(R)` | `QM` | yes | 206 |
| `C x C` | `CM` | yes | 72 |

At 0.35 s/curve, the full 6.2M rows is about 600 core-hours at `B=50`: three
hours on 200 cores.

### Stratified sample

`(expected, meta)` determines the stratum, so a sample can be cut straight out
of the prepared TSV. `sample_corpus.py` does that, reproducibly:

```sh
./sample_corpus.py prepared/g2.tsv --list-strata

./sample_corpus.py prepared/g2.tsv -o prepared/g2_sample.tsv --seed 0 \
    --default 200 \
    --take 'M_2(RR)|QM=all' --take 'CC,CC|CM=all' --take 'CC,CC|CM x CM=all'
```

That is the exact command behind `prepared/g2_sample.tsv`: 1,536 curves, all
206 QM and all 530 `C x C`, 200 from each of the other four strata. Re-running
it reproduces the file byte for byte.

A stratum spec is `EXPECTED|META=COUNT`, the two column values verbatim
(`META` is empty for g3, so the key ends in `|`), and `COUNT` is a row count or
the word `all`. Strata you do not name take `--default`; a stratum smaller than
its count is taken whole; naming a stratum that does not exist is an error
rather than a silent no-op, so a typo in a key cannot quietly shrink a sample.

Selection is by `blake2b(seed || id)`, keeping the smallest digests, not by
`random.sample`. Which rows a PRNG hands back depends on the order the pools
were built in and on the CPython version's sampling algorithm, so a seeded
`random.sample` recipe does not in fact reproduce its own sample across
machines. Hashing the id depends on nothing but the seed and the id, so the
same `(seed, spec, input)` always yields the same id set, and adding curves to
the corpus only adds candidates instead of reshuffling everything.

## Reading the report

Verdicts, computed in Python (not in Magma, so the comparison does not depend on
Magma and Python sorting strings identically; both sides are re-sorted before
comparing):

| verdict | meaning |
| --- | --- |
| `sharp` | the multisets are equal |
| `over` | got is a genuine upper bound: expected embeds unitally in got, but the two are not equal |
| `under` | `dim(got) < dim(expected)` |
| `mismatch` | equal or larger dimension but structurally incomparable |
| `error` | the driver raised on this curve |

Dimensions over R: `RR` 1, `CC` 2, `HH` 4, `M_2(RR)` 4, `M_2(CC)` 8,
`M_3(RR)` 9, `M_3(CC)` 18, `M_2(RR) or HH` 4. Greater total dimension is not
enough for `over`: `CC,RR` does not embed in `M_2(RR)`, so that pair is a
`mismatch`.

`over` and `mismatch` are the ones a larger `B` can move. `under` is not: it
means the bound missed endomorphisms that are really there.

### The expected `under`

For geometrically simple genus-2 curves whose quartic CM is defined only over an
extension of Q, the extra endomorphisms are invisible to prime-field Frobenius:
the bound returns `[RR, RR]` where the truth is `[CC, CC]`. Because eta is a
`min` over primes, raising `B` does not help. This is an inherent limitation of
taking Frobenius data at primes of Q, not a harness bug, which is why the report
counts these rows rather than treating them as failures. Expect them in the
`C x C` strata (530 curves in the non-generic corpus).

A row that reports `under` outside those strata is worth looking at.

## What this corpus actually measured

Numbers from an 8-core run, kept here because they are what the two-pass advice
above is based on.

**Genus 3, all 4,341 quartics.** At `B=50`: 94.6% sharp, 236 `over`, no `under`,
no `mismatch`, no errors. Re-running those 236 at `B=200` fixed 231 of them;
the last 5 (all in the `RR,RR,RR` stratum) went sharp at `B=500`. Three of those
5 are the non-monotonicity examples above: `B=200` made them worse before
`B=500` made them sharp. Three passes reach **100% sharp**, and 94.6% of the
corpus never leaves the cheap pass. Cost
climbs steeply with `B` for plane quartics: 0.88 s/curve at `B=50`, 14.5 s at
`B=200`, 57 s at `B=500`, so paying it only for the survivors is the whole game.

**Genus 2, 1,536-curve stratified sample** (all 206 QM, all 530 `C x C`, 200
each from the other four strata). This sample predates `sample_corpus.py` and
was drawn by a recipe that did not reproduce itself, so it is kept alongside
the results as `prepared/g2_sample.legacy.tsv`; the seed-0 sample that
`sample_corpus.py` now produces has the same stratum shape and shares 778 of
its 1,536 ids. At `B=50`: 94.1% sharp, 19 `over`, 72 `under`.
All 19 `over` went sharp at `B=200`. The 72 `under` are exactly the
geometrically simple quartic-CM curves described above; they are the entire
residue, so the sample is sharp everywhere the method can be sharp.

Worth knowing: the `under` set is **only** the 72 geometrically simple CM
curves, not both `C x C` strata. All 458 `CM x CM` curves come back sharp at
`B=50`. A product of two CM elliptic curves is visible to prime-field Frobenius
even though a simple abelian surface with the same real algebra is not.

## Behaviour under failure

* Every curve is wrapped in `try`/`catch`. A curve that raises gets status
  `error` with the message in the `got` field, and the chunk carries on. No
  asserts, no `SetQuitOnError`.
* Each chunk runs under `timeout` (`-t`, default 3600s), so a pathological curve
  cannot wedge a worker. Because results are flushed per line, the chunk resumes
  from where it stopped.
* A chunk that exits non-zero is logged in the work directory and does not abort
  the run. Only the gp guard, which is checked once up front by `run_corpus.sh`
  before any worker starts, stops the whole thing.
* A dispatcher that could not start the chunks at all is different: the run
  still harvests and reports, keeps the work directory, and then exits **1**,
  so a pass that computed nothing cannot be mistaken for a pass with nothing
  to do.
* Magma ignores the return value of the `write(2)` behind `PrintFile`, so a
  full or over-quota disk can leave a chunk file ending mid-row. The fold reads
  chunk files record by record rather than concatenating them, so a partial row
  is dropped on its own and that curve is simply re-run; it can never fuse with
  the next chunk's first row into a well-formed row for another id.
* `RealRepresentationBound` returning `[]` means no bound was established at
  all. That is recorded as an empty `got` with status `ok`; the report counts it
  as `under` and prints a separate note saying how many there were.
