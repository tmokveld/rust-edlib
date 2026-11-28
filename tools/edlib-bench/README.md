# Edlib benchmarks

The [Criterion](https://criterion-rs.github.io/book/) suite compares the working
checkout of `rust-edlib` with upstream C++ Edlib 1.3.9 and, optionally, with a
previous `rust-edlib` revision. Run these commands from the repository root.

## Setup and quick runs

Use stable Rust and a C++ compiler. Clone upstream Edlib once:

```sh
git clone --depth 1 --branch python-1.3.9 \
  https://github.com/Martinsos/edlib.git target/edlib-upstream
sh tools/edlib-bench/run.sh --smoke
```

Set `EDLIB_SOURCE` if your upstream checkout is elsewhere. Omit `--smoke` for
the full panel, which includes million-byte inputs. Pass a Criterion filter to
select a group (`bounded-global`, `anchored-score`, `infix-score`, or `alignment`):

```sh
sh tools/edlib-bench/run.sh anchored-score --sample-size 10
```

## Compare with a Git revision

```sh
sh tools/edlib-bench/run.sh --smoke --against 09c70ae
sh tools/edlib-bench/run.sh --against HEAD alignment --sample-size 10
```

Put `--against <ref>` before any Criterion arguments. The ref must resolve in
the local Git repository and expose the alignment and score-only APIs used by
the current suite. `09c70ae` is the pinned CI baseline; `HEAD` is useful for
comparing uncommitted changes with the current commit. The runner builds the
selected revision in a temporary directory without changing the checkout and
uses the same harness, cases, and release settings for both Rust builds.

The runner prints a Markdown table of median nanoseconds per call to stdout;
Criterion progress and confidence intervals go to stderr. With `--against`,
the table records the resolved commit and adds `Baseline ns/call` and
`Baseline / Rust`. A ratio above 1× favors the current checkout. The
[CI benchmark job](../../.github/workflows/ci.yml) appends this table to the job
summary.

## Native Criterion output

For native Criterion output or HTML reports without a Git comparison, invoke
Cargo directly:

```sh
EDLIB_SOURCE="$PWD/target/edlib-upstream" \
  cargo +stable bench --locked --manifest-path tools/edlib-bench/Cargo.toml \
  --bench compare --features html_reports -- anchored-score
```
