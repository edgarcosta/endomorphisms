# Upperbounds tests

Thanks for helping check the upperbounds implementation. These tests need a
licensed Magma installation, PARI/gp on PATH, and the MagmaPolred spec.

- `upperbounds/unit.m`: polynomial operations, field intersections, algebra
  bounds, and curve-level examples.
- `upperbounds/lmfdb.m`: LMFDB examples across real endomorphism algebra types.
- `upperbounds/regressions.m`: previous bugs, ladder exhaustion, and hard cases.
- `upperbounds/corpus-cells.m`: the smallest-conductor curve in each of the 27
  populated (geom_end_alg, B-rung) cells from the 3,420,837-curve validation run.

From the repository root, run all tests or name a subset:

```sh
export PATH="/home/sage/sage-10.8/local/bin:$PATH"
export POLRED_SPEC=/home/edgarcosta/projects/CHIMP/CHIMP/MagmaPolred/spec
./tests/run.sh
./tests/run.sh unit
./tests/run.sh lmfdb corpus-cells
```

The runner also works from any directory when invoked by its path. Names may
include `.m`. Your PATH and POLRED_SPEC take precedence. If unset, PATH falls
back to the Sage directory above plus standard system directories, and
POLRED_SPEC falls back to the spec above. Set either to use another installation.
Missing gp causes a clear error before any tests run; center-field detection
needs PARI/gp through Polredabs. An unreadable spec also fails before testing.

A test passes only with exit status 0 **and zero output bytes**, counting stdout
and stderr together. Magma errors and dependency skip notices can produce output
without a failing exit status, so both checks are required. Tests are silent on
success; failures show their captured output. The runner reports elapsed seconds
for each test and a summary, and exits nonzero if any test fails.

Magma is licensed and cannot run on public GitHub runners; CI wiring is left to whoever has a licensed runner.
