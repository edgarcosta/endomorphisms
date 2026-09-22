#!/usr/bin/env bash

set -uo pipefail

fake_magma() {
    if [[ "$(head -n 1 "$CORPUS_INPUT")" == $'id\tkind\tdata\texpected' ]]; then
        : > "$CORPUS_OUTPUT"
        return 0
    fi

    local id kind data expected got
    while IFS=$'\t' read -r id kind data expected; do
        [[ -n "$id" ]] || continue
        got="none"
        [[ "$FAKE_STATUS" == "ok" ]] && got="computed"
        printf '%s\t%s\t%s\t%s\t%s\t0.01\t%s\n' \
            "$id" "$kind" "$CORPUS_B" "$expected" "$got" "$FAKE_STATUS" \
            >> "$CORPUS_OUTPUT"
        printf '%s\n' "$id" >> "$FAKE_RECORD"
    done < "$CORPUS_INPUT"
}

if [[ "${RUN_CORPUS_FAKE_MAGMA:-0}" == "1" ]]; then
    fake_magma
    exit 0
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$HERE/$(basename "${BASH_SOURCE[0]}")"
TARGET="${1:-$HERE/run_corpus.sh}"

if [[ ! -r "$TARGET" ]]; then
    printf 'test_run_corpus.sh: cannot read %s\n' "$TARGET" >&2
    exit 1
fi

ROOT="$(mktemp -d "${TMPDIR:-/tmp}/test-run-corpus.XXXXXX")"
trap 'rm -rf "$ROOT"' EXIT
RUNNER="$ROOT/run_corpus.sh"
cp "$TARGET" "$RUNNER"
cp "$HERE/run_corpus.m" "$ROOT/run_corpus.m"
chmod +x "$RUNNER"

failures=0

fail() {
    printf '%s\n' "$1" >&2
    failures=$((failures + 1))
}

assert_contents() {
    local label="$1" file="$2" expected="$3"
    if [[ ! -f "$file" ]]; then
        fail "$label: missing $file"
        return
    fi
    if ! cmp -s <(printf '%s' "$expected") "$file"; then
        fail "$label: output differed"
        diff -u <(printf '%s' "$expected") "$file" >&2 || true
    fi
}

assert_empty() {
    local label="$1" file="$2"
    if [[ -s "$file" ]]; then
        fail "$label: expected no planned ids, got:"
        sed 's/^/  /' "$file" >&2
    fi
}

run_target() {
    local label="$1" status="$2" record="$3" log="$4"
    shift 4
    local rc=0
    RUN_CORPUS_FAKE_MAGMA=1 FAKE_STATUS="$status" FAKE_RECORD="$record" \
        MAGMA="$SELF" "$RUNNER" "$@" > "$log" 2>&1 || rc=$?
    if [[ "$rc" -ne 0 ]]; then
        fail "$label: run exited $rc"
        sed 's/^/  /' "$log" >&2
    fi
}

header=$'id\tkind\tB\texpected\tgot\tcputime\tstatus\n'

unsound_dir="$ROOT/unsound"
mkdir -p "$unsound_dir/work"
printf 'id\tkind\tdata\texpected\nunsound-id\tg2\t[]\ttruth\n' > "$unsound_dir/input.tsv"
printf 'unsound-id\tg2\t50\ttruth\tbound\t0.01\tunsound\n' \
    > "$unsound_dir/work/out.000000.tsv"
: > "$unsound_dir/record"
run_target "unsound harvest" ok "$unsound_dir/record" "$unsound_dir/run.log" \
    -i "$unsound_dir/input.tsv" -o "$unsound_dir/output.tsv" \
    -w "$unsound_dir/work" -B 50 -j 1
assert_contents "unsound harvest" "$unsound_dir/output.tsv" \
    "${header}"$'unsound-id\tg2\t50\ttruth\tbound\t0.01\tunsound\n'
assert_empty "unsound harvest" "$unsound_dir/record"

retry_dir="$ROOT/retry"
mkdir -p "$retry_dir"
printf 'id\tkind\tdata\texpected\nretry-id\tg2\t[]\ttruth\n' > "$retry_dir/input.tsv"
: > "$retry_dir/record"
run_target "initial error" error "$retry_dir/record" "$retry_dir/first.log" \
    -i "$retry_dir/input.tsv" -o "$retry_dir/output.tsv" -B 50 -j 1
assert_contents "initial error" "$retry_dir/output.tsv" \
    "${header}"$'retry-id\tg2\t50\ttruth\tnone\t0.01\terror\n'

: > "$retry_dir/record"
run_target "error retry" ok "$retry_dir/record" "$retry_dir/second.log" \
    -i "$retry_dir/input.tsv" -o "$retry_dir/output.tsv" -B 50 -j 1
assert_contents "error retry plan" "$retry_dir/record" $'retry-id\n'
assert_contents "successful retry precedence" "$retry_dir/output.tsv" \
    "${header}"$'retry-id\tg2\t50\ttruth\tcomputed\t0.01\tok\n'

boundary_dir="$ROOT/boundary"
mkdir -p "$boundary_dir/work"
printf 'id\tkind\tdata\texpected\nvictim\tg2\t[]\ttruth\n' > "$boundary_dir/input.tsv"
printf 'forged-prefix' > "$boundary_dir/work/out.000001.tsv"
printf 'victim\tg2\t50\ttruth\tbound\t0.01\tok\n' \
    > "$boundary_dir/work/out.000002.tsv"
: > "$boundary_dir/record"
run_target "part boundary" ok "$boundary_dir/record" "$boundary_dir/run.log" \
    -i "$boundary_dir/input.tsv" -o "$boundary_dir/output.tsv" \
    -w "$boundary_dir/work" -B 50 -j 1
assert_contents "part boundary" "$boundary_dir/output.tsv" \
    "${header}"$'victim\tg2\t50\ttruth\tbound\t0.01\tok\n'
assert_empty "part boundary" "$boundary_dir/record"

guard_dir="$ROOT/guard"
mkdir -p "$guard_dir/work"
printf 'id\tkind\tdata\texpected\nguarded\tg2\t[]\ttruth\n' > "$guard_dir/input.tsv"
: > "$guard_dir/record"
guard_rc=0
RUN_CORPUS_FAKE_MAGMA=1 FAKE_STATUS=ok FAKE_RECORD="$guard_dir/record" \
    MAGMA="$SELF" "$RUNNER" -i "$guard_dir/input.tsv" \
    -o "$guard_dir/work/output.tsv" -w "$guard_dir/work" -B 50 -j 1 \
    > "$guard_dir/run.log" 2>&1 || guard_rc=$?
if [[ "$guard_rc" -ne 64 ]]; then
    fail "output guard: expected exit 64, got $guard_rc"
    sed 's/^/  /' "$guard_dir/run.log" >&2
fi
if [[ -e "$guard_dir/work/output.tsv" ]]; then
    fail "output guard: created an output file inside the work directory"
fi
assert_empty "output guard" "$guard_dir/record"

[[ "$failures" -eq 0 ]] || exit 1
