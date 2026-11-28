#!/bin/sh
set -eu

bench=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
EDLIB_SOURCE=${EDLIB_SOURCE:-$bench/../../target/edlib-upstream}
case "$EDLIB_SOURCE" in /*) ;; *) EDLIB_SOURCE="$PWD/$EDLIB_SOURCE" ;; esac
RUSTUP_TOOLCHAIN=${RUSTUP_TOOLCHAIN:-stable}
export EDLIB_SOURCE RUSTUP_TOOLCHAIN

smoke=false
against=
while [ "$#" -gt 0 ]; do
    case "$1" in
        --smoke) smoke=true; shift ;;
        --against)
            if [ "$#" -lt 2 ] || [ -z "$2" ]; then
                echo 'usage: run.sh [--smoke] [--against <git-ref>] [Criterion arguments...]' >&2
                exit 2
            fi
            against=$2
            shift 2
            ;;
        *) break ;;
    esac
done
if "$smoke"; then
    filter='^(bounded-global/substitutions/global/score/64x64/k=[78]/'
    filter="$filter|anchored-score/free-flank/(prefix|suffix)/score/65x81/k=8/"
    filter="$filter|infix-score/late-exact/|alignment/mixed-edits-1%/global/(score|path)/[0-9]+x500/)"
    set -- "$filter" --sample-size 10 --warm-up-time 0.005 --measurement-time 0.02 "$@"
fi

# Commands that don't save new estimates (or explicitly request native output)
# keep Criterion's normal stdout instead of producing a comparison table.
table=true
for arg do
    case "$arg" in
        --) break ;;
        -h|--help|--list|--test|--profile-time|--profile-time=*|\
        --discard-baseline|--load-baseline|--load-baseline=*|--output-format|--output-format=*)
            table=false
            break
            ;;
    esac
done

if ! "$table"; then
    if [ -n "$against" ]; then
        echo '--against requires a Criterion run that saves new estimates' >&2
        exit 2
    fi
    exec cargo bench --locked --manifest-path "$bench/Cargo.toml" --bench compare -- "$@"
fi

target_dir=${CARGO_TARGET_DIR:-$bench/target}
case "$target_dir" in /*) ;; *) target_dir="$PWD/$target_dir" ;; esac
CRITERION_HOME=${CRITERION_HOME:-$target_dir/criterion}
case "$CRITERION_HOME" in /*) ;; *) CRITERION_HOME="$PWD/$CRITERION_HOME" ;; esac
export CRITERION_HOME

scratch=$(mktemp -d "${TMPDIR:-/tmp}/edlib-bench.XXXXXXXX")
trap 'rm -rf "$scratch"' 0

if [ -n "$against" ]; then
    repo=$bench/../..
    revision=$(git -C "$repo" rev-parse --verify --end-of-options "$against^{commit}")
    baseline=$scratch/baseline
    baseline_bench=$baseline/tools/edlib-bench
    mkdir -p "$baseline_bench"
    # Benchmark the committed library with this checkout's harness and settings.
    git -C "$repo" archive --output="$scratch/baseline.tar" "$revision" -- Cargo.toml src
    tar -xf "$scratch/baseline.tar" -C "$baseline"
    cp -R "$bench/Cargo.toml" "$bench/Cargo.lock" "$bench/build.rs" \
        "$bench/upstream.cpp" "$bench/benches" "$baseline_bench/"
    printf 'Benchmarking rust-edlib baseline %s\n' "$revision" >&2
    (
        unset EDLIB_BENCH_SUMMARY
        # Cargo may need to update the Rust package version in the scratch lockfile.
        CARGO_TARGET_DIR="$target_dir" CRITERION_HOME="$scratch/criterion" EDLIB_BENCH_RUST_ONLY=1 \
            cargo bench --manifest-path "$baseline_bench/Cargo.toml" --bench compare -- "$@"
    ) >&2
    export EDLIB_BENCH_BASELINE_HOME="$scratch/criterion" EDLIB_BENCH_BASELINE_REV="$revision"
    printf '\nBenchmarking the current checkout and upstream Edlib\n' >&2
fi

summary=$scratch/summary.md
EDLIB_BENCH_SUMMARY="$summary" \
    cargo bench --locked --manifest-path "$bench/Cargo.toml" --bench compare -- "$@" >&2
cat "$summary"
