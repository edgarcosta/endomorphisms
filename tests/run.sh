#!/usr/bin/env bash
set -euo pipefail

# Bash supplies an unexported default PATH when none was inherited.
if [[ $(declare -p PATH) != "declare -x "* ]]; then
    PATH=/home/sage/sage-10.8/local/bin:/usr/local/bin:/usr/bin:/bin
fi
export PATH
export POLRED_SPEC="${POLRED_SPEC-/home/edgarcosta/projects/CHIMP/CHIMP/MagmaPolred/spec}"

if ! command -v gp >/dev/null 2>&1; then
    printf 'PARI/gp is required: add the directory containing gp to PATH.\n' >&2
    exit 1
fi
if ! command -v magma >/dev/null 2>&1; then
    printf 'Magma is required: add the directory containing magma to PATH.\n' >&2
    exit 1
fi
if [[ ! -f "$POLRED_SPEC" || ! -r "$POLRED_SPEC" ]]; then
    printf 'POLRED_SPEC must name a readable MagmaPolred spec file: %s\n' "$POLRED_SPEC" >&2
    exit 1
fi
if [[ "$POLRED_SPEC" != /* ]]; then
    export POLRED_SPEC="$PWD/$POLRED_SPEC"
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
utils="$script_dir/../examples/Utils"
test_dir="$script_dir/upperbounds"
tests=()
if (( $# )); then
    for name in "$@"; do
        name=${name%.m}
        if [[ "$name" == */* || ! -f "$test_dir/$name.m" ]]; then
            printf 'Unknown test: %s (use names such as unit or corpus-cells).\n' "$name" >&2
            exit 2
        fi
        tests+=("$test_dir/$name.m")
    done
else
    shopt -s nullglob
    tests=("$test_dir"/*.m)
fi
if (( ${#tests[@]} == 0 )); then
    printf 'No upperbounds tests found in %s.\n' "$test_dir" >&2
    exit 1
fi

output_dir=$(mktemp -d "${TMPDIR:-/tmp}/endomorphisms-tests.XXXXXX")
trap 'rm -rf -- "$output_dir"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
cd -- "$test_dir"

elapsed() {
    awk -v start="$1" -v end="$(date +%s.%N)" 'BEGIN { printf "%.3f", end - start }'
}

passed=0
failed=0
total_start=$(date +%s.%N)
for test in "${tests[@]}"; do
    name=${test##*/}
    name=${name%.m}
    output="$output_dir/$name.log"
    start=$(date +%s.%N)
    code=0
    magma -b "$utils/SetQuitOnErrortrue.m" "$test" "$utils/Exit.m" >"$output" 2>&1 || code=$?
    seconds=$(elapsed "$start")
    bytes=$(wc -c < "$output")
    if (( code == 0 )) && [[ ! -s "$output" ]]; then
        printf 'PASS %s %ss (exit=%d, bytes=%s)\n' "$name" "$seconds" "$code" "$bytes"
        passed=$((passed + 1))
    else
        printf 'FAIL %s %ss (exit=%d, bytes=%s)\n' "$name" "$seconds" "$code" "$bytes"
        cat -- "$output"
        failed=$((failed + 1))
    fi
done
printf 'Summary: %d passed, %d failed, %d total in %ss\n' \
    "$passed" "$failed" "${#tests[@]}" "$(elapsed "$total_start")"
(( failed == 0 ))
