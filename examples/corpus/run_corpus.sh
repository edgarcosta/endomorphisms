#!/usr/bin/env bash
#
# Shard a prepared corpus TSV across N Magma processes and collect the results.
#
#   ./run_corpus.sh -i g2.tsv -o results/g2.B50.tsv -B 50 -j 8
#
# Resumable: results are appended one flushed line at a time by run_corpus.m,
# and every invocation first harvests whatever the previous one left in the
# work directory, then runs only the ids that have no row yet at this B.
# Killing the run (or a chunk hitting its timeout) costs at most the curves in
# flight; re-running the same output file at a larger B runs the pass again.
#
# Requires a working PARI/gp on PATH; see README.md. The preflight below
# refuses to start without one, because a missing gp does not fail loudly, it
# silently degrades every answer.
#
set -euo pipefail

SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
HERE="$(dirname "$SELF")"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

# Both spec paths are overridable; the endomorphisms one defaults off this
# script's own location so the harness moves with the checkout.
: "${ENDO_SPEC:=$REPO_ROOT/endomorphisms/magma/spec}"
: "${POLRED_SPEC:=}"
: "${MAGMA:=magma}"
export ENDO_SPEC POLRED_SPEC

DRIVER="$HERE/run_corpus.m"

# ---------------------------------------------------------------------------
# Internal mode: run exactly one chunk. Re-invokes this script so that both
# GNU parallel and xargs can dispatch it without exported shell functions.
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "--run-chunk" ]]; then
    chunk="$2"
    base="$(basename "$chunk" .tsv)"
    dir="$(dirname "$chunk")"
    part="$dir/out.${base#chunk.}.tsv"
    log="$dir/log.${base#chunk.}.txt"
    rc=0
    CORPUS_INPUT="$chunk" CORPUS_OUTPUT="$part" CORPUS_B="${CORPUS_B:-50}" \
        timeout --foreground "${CHUNK_TIMEOUT:-3600}" \
        "$MAGMA" -b "$DRIVER" >"$log" 2>&1 || rc=$?
    if [[ $rc -eq 2 ]]; then
        echo "run_corpus.sh: $base aborted (gp guard or setup); see $log" >&2
    elif [[ $rc -eq 124 ]]; then
        echo "run_corpus.sh: $base timed out after ${CHUNK_TIMEOUT:-3600}s; partial results kept" >&2
    elif [[ $rc -ne 0 ]]; then
        echo "run_corpus.sh: $base exited $rc; see $log" >&2
    fi
    exit 0
fi

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
usage() {
    cat <<'EOF'
Usage: run_corpus.sh -i INPUT.tsv -o RESULTS.tsv [options]

  -i, --input FILE     prepared TSV from prepare_input.py            (required)
  -o, --output FILE    results TSV, created or resumed               (required)
  -B, --bound N        prime bound for RealRepresentationBound       (50)
  -j, --jobs N         concurrent Magma processes                    (nproc)
  -c, --chunk N        curves per Magma process                      (200)
  -t, --timeout SEC    wall-clock limit per chunk                    (3600)
  -w, --workdir DIR    scratch dir, deleted at the end                (OUTPUT.work)
      --ids FILE       only run ids listed in FILE, one per line
                       (second pass: report.py --list-nonsharp)
      --keep-work      do not delete the work directory when finished
  -h, --help           this message

Environment: ENDO_SPEC, POLRED_SPEC, MAGMA.
EOF
}

INPUT=""; OUTPUT=""; B=50; JOBS="$(nproc 2>/dev/null || echo 4)"
CHUNK=200; TIMEOUT_S=3600; WORKDIR=""; IDS=""; KEEP_WORK=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)   INPUT="$2"; shift 2 ;;
        -o|--output)  OUTPUT="$2"; shift 2 ;;
        -B|--bound)   B="$2"; shift 2 ;;
        -j|--jobs)    JOBS="$2"; shift 2 ;;
        -c|--chunk)   CHUNK="$2"; shift 2 ;;
        -t|--timeout) TIMEOUT_S="$2"; shift 2 ;;
        -w|--workdir) WORKDIR="$2"; shift 2 ;;
        --ids)        IDS="$2"; shift 2 ;;
        --keep-work)  KEEP_WORK=1; shift ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "run_corpus.sh: unknown option $1" >&2; usage >&2; exit 64 ;;
    esac
done

[[ -n "$INPUT"  ]] || { echo "run_corpus.sh: --input is required" >&2;  exit 64; }
[[ -n "$OUTPUT" ]] || { echo "run_corpus.sh: --output is required" >&2; exit 64; }
[[ -r "$INPUT"  ]] || { echo "run_corpus.sh: cannot read $INPUT" >&2;   exit 66; }
[[ -r "$DRIVER" ]] || { echo "run_corpus.sh: cannot read $DRIVER" >&2;  exit 66; }

mkdir -p "$(dirname "$OUTPUT")"
WORKDIR="${WORKDIR:-$OUTPUT.work}"
mkdir -p "$WORKDIR"

# The work directory is deleted wholesale when the pass completes, so a results
# file inside it is unsafe at every exit path; refusing before any work starts
# is the only check that cannot itself lose data.
OUTPUT_ABS="$(cd "$(dirname "$OUTPUT")" && pwd -P)/$(basename "$OUTPUT")"
WORKDIR_ABS="$(cd "$WORKDIR" && pwd -P)"
if [[ "$OUTPUT_ABS" == "${WORKDIR_ABS%/}/"* ]]; then
    echo "run_corpus.sh: --output $OUTPUT lies inside --workdir $WORKDIR," >&2
    echo "  which is deleted when the pass completes. Point -w somewhere else." >&2
    exit 64
fi

export CORPUS_B="$B" CHUNK_TIMEOUT="$TIMEOUT_S"

RESULT_HEADER=$'id\tkind\tB\texpected\tgot\tcputime\tstatus'
[[ -s "$OUTPUT" ]] || printf '%s\n' "$RESULT_HEADER" > "$OUTPUT"

# ---------------------------------------------------------------------------
# Fold any chunk output left on disk into $OUTPUT, keeping one row per
# (id, kind, B). Called before planning and after running, so a killed run is
# picked up by the next invocation.
# ---------------------------------------------------------------------------
harvest() {
    shopt -s nullglob
    local parts=("$WORKDIR"/out.*.tsv)
    shopt -u nullglob
    [[ ${#parts[@]} -gt 0 ]] || return 0

    local batch tmp ok=1
    batch="$(mktemp "$WORKDIR/.batch.XXXXXX")"
    tmp="$(mktemp "$OUTPUT.harvest.XXXXXX")"

    # awk, not cat: PrintFile writes a row with one write(2) whose short count
    # Magma ignores, so a part file can end mid-row with no newline, and cat
    # would fuse that partial row onto the next part's first row. Cut before
    # the row's first tab, the fusion is a well-formed 7-field row carrying a
    # different curve's id, which no field-level check can detect. awk ends a
    # record at every file boundary, so the two rows stay separate.
    # xargs, not one argv entry per part: 40k parts overflow ARG_MAX.
    printf '%s\0' "${parts[@]}" | xargs -0 awk 1 > "$batch" || ok=0

    # Rewrite rather than append, and read the staged rows as a separate awk
    # operand: $OUTPUT is then replaced atomically and can never itself be left
    # ending mid-row for the next harvest to fuse onto.
    if [[ $ok -eq 1 ]]; then
        {
            printf '%s\n' "$RESULT_HEADER"
            # The status check also drops a line truncated by a failed write,
            # so that id is simply re-run on the next pass.
            # For duplicate keys, status precedence is ok > unsound > error.
            # A later ok repairs either failure, while a later error cannot
            # erase unsound; equal statuses keep the first row.
            awk -F'\t' '
                function rank(status) {
                    if (status == "ok") return 3
                    if (status == "unsound") return 2
                    if (status == "error") return 1
                    return 0
                }
                NF == 7 {
                    priority = rank($7)
                    if (priority == 0) next
                    key = $1 FS $2 FS $3
                    if (!(key in position)) {
                        position[key] = ++count
                        keys[count] = key
                    }
                    if (!(key in priorities) || priority > priorities[key]) {
                        rows[key] = $0
                        priorities[key] = priority
                    }
                }
                END { for (i = 1; i <= count; i++) print rows[keys[i]] }
            ' "$OUTPUT" "$batch"
        } > "$tmp" || ok=0
    fi

    if [[ $ok -eq 0 ]]; then
        rm -f "$batch" "$tmp"
        echo "run_corpus.sh: harvest failed; chunk output left in $WORKDIR" >&2
        return 1
    fi

    # mktemp makes the replacement 0600; the results file keeps its own mode.
    chmod --reference="$OUTPUT" "$tmp"
    mv "$tmp" "$OUTPUT"
    rm -f "$batch"
    # Only now, with the new results file in place, are the parts redundant.
    printf '%s\0' "${parts[@]}" | xargs -0 rm -f --
}

# Never fatal: by the time this runs the pass has already reported its result,
# and a work directory that outlives it is a nuisance, not a failure.
drop_workdir() {
    rm -rf "$WORKDIR" || echo "run_corpus.sh: could not remove $WORKDIR" >&2
}

harvest

# ---------------------------------------------------------------------------
# Preflight: run the driver over an empty chunk. That exercises the AttachSpec
# paths and the gp round-trip guard exactly as the workers will, and costs one
# Magma startup.
# ---------------------------------------------------------------------------
printf '%s\n' $'id\tkind\tdata\texpected' > "$WORKDIR/empty.tsv"
if ! CORPUS_INPUT="$WORKDIR/empty.tsv" CORPUS_OUTPUT="$WORKDIR/preflight.out" \
        "$MAGMA" -b "$DRIVER"; then
    echo "run_corpus.sh: preflight failed; refusing to run the corpus." >&2
    exit 2
fi
rm -f "$WORKDIR/preflight.out"

# ---------------------------------------------------------------------------
# Plan: everything in $INPUT that does not already have a row in $OUTPUT AT
# THIS B. Keying the done set on the id alone would make a second pass at a
# larger B into a silent no-op on any output file that already holds the cheap
# rows, and would disagree with harvest(), which dedups on (id, kind, B).
# ---------------------------------------------------------------------------
PENDING="$WORKDIR/pending.tsv"
if [[ -n "$IDS" ]]; then
    [[ -r "$IDS" ]] || { echo "run_corpus.sh: cannot read $IDS" >&2; exit 66; }
    awk -F'\t' -v want="$IDS" -v bound="$B" '
        BEGIN { while ((getline line < want) > 0) { sub(/\r$/, "", line); if (line != "") keep[line] = 1 } }
        FNR == NR { if (FNR > 1 && NF >= 7 && $7 != "error") done[$1 SUBSEP $3] = 1; next }
        FNR == 1 && $1 == "id" { next }
        ($1 in keep) && !(($1 SUBSEP bound) in done)
    ' "$OUTPUT" "$INPUT" > "$PENDING"
else
    awk -F'\t' -v bound="$B" '
        FNR == NR { if (FNR > 1 && NF >= 7 && $7 != "error") done[$1 SUBSEP $3] = 1; next }
        FNR == 1 && $1 == "id" { next }
        !(($1 SUBSEP bound) in done)
    ' "$OUTPUT" "$INPUT" > "$PENDING"
fi

TOTAL=$(wc -l < "$PENDING")
DONE_ALREADY=$(awk -F'\t' -v bound="$B" \
    'FNR > 1 && NF >= 7 && $3 == bound && $7 != "error" { n++ } END { print n + 0 }' "$OUTPUT")
echo "run_corpus.sh: $DONE_ALREADY already done at B=$B, $TOTAL to run on $JOBS jobs"

if [[ "$TOTAL" -eq 0 ]]; then
    [[ "$KEEP_WORK" -eq 1 ]] || drop_workdir
    exit 0
fi

# find, not a glob: at small chunk sizes the previous pass can leave tens of
# thousands of chunk files, which overflow rm's argument list.
find "$WORKDIR" -maxdepth 1 -type f -name 'chunk.*.tsv' -delete
split -l "$CHUNK" -d -a 6 --additional-suffix=.tsv "$PENDING" "$WORKDIR/chunk."

shopt -s nullglob
CHUNKS=("$WORKDIR"/chunk.*.tsv)
shopt -u nullglob
echo "run_corpus.sh: ${#CHUNKS[@]} chunks of up to $CHUNK curves"

# ---------------------------------------------------------------------------
# Dispatch. Failures are per chunk and never abort the run; whatever a chunk
# managed to write is already on disk, so the status is remembered and reported
# after the harvest rather than aborting before it.
# ---------------------------------------------------------------------------
DISPATCH_RC=0
if command -v parallel >/dev/null 2>&1; then
    # --quote, because parallel builds a shell command line out of its
    # arguments: without it a $SELF containing a space word-splits and no
    # worker starts at all.
    printf '%s\n' "${CHUNKS[@]}" \
        | parallel --will-cite --quote -j "$JOBS" --halt never --line-buffer \
              "$SELF" --run-chunk {} || DISPATCH_RC=$?
else
    echo "run_corpus.sh: GNU parallel not found, falling back to xargs -P" >&2
    printf '%s\n' "${CHUNKS[@]}" \
        | xargs -P "$JOBS" -I{} "$SELF" --run-chunk {} || DISPATCH_RC=$?
fi

harvest

FINAL=$(( $(wc -l < "$OUTPUT") - 1 ))
MISSING=$(awk -F'\t' -v bound="$B" '
    FNR == NR { if (FNR > 1 && NF >= 7 && $7 != "error") done[$1 SUBSEP $3] = 1; next }
    FNR == 1 && $1 == "id" { next }
    !(($1 SUBSEP bound) in done) { n++ }
    END { print n + 0 }
' "$OUTPUT" "$PENDING")
echo "run_corpus.sh: $FINAL rows in $OUTPUT ($MISSING of this pass still missing)"

if [[ "$DISPATCH_RC" -ne 0 ]]; then
    echo "run_corpus.sh: dispatch reported failures (exit $DISPATCH_RC); see the logs in $WORKDIR" >&2
fi

if [[ "$KEEP_WORK" -eq 1 || "$MISSING" -gt 0 || "$DISPATCH_RC" -ne 0 ]]; then
    echo "run_corpus.sh: work directory kept at $WORKDIR"
else
    drop_workdir
fi

# A chunk that failed on its own is normal and already reported above, but a
# dispatcher that could not start the chunks at all must not look like success.
[[ "$DISPATCH_RC" -eq 0 ]] || exit 1
