# rust-edlib

Exact, unit-cost Levenshtein alignment over byte slices, implemented in safe Rust
with Edlib semantics.

## Installation

Use Rust 1.85 or later. Add the checked-out library to your project's
`Cargo.toml`, adjusting the path relative to that manifest:

```toml
[dependencies]
rust-edlib = { path = "../rust-edlib" }
```

## Example

For an exact global score, `rust_edlib::distance(query, target)?` returns a
`usize` without computing locations or an edit path.

For a bounded search, other alignment modes, or an edit path, use `align`:

```rust
use rust_edlib::{AlignConfig, AlignError, AlignmentMode, AlignmentTask, CigarFormat, align};

fn main() -> Result<(), AlignError> {
    let config = AlignConfig {
        mode: AlignmentMode::Infix,
        task: AlignmentTask::Path,
        max_distance: Some(0),
        ..Default::default()
    };
    let result = align(b"AC", b"TACGAC", &config)?;

    assert_eq!(result.edit_distance(), Some(0));
    assert!(result.locations().unwrap().eq([1..3, 4..6]));
    println!("{}", result.cigar(CigarFormat::Extended).unwrap()); // 2=
    Ok(())
}
```

## Modes, tasks, and options

Every mode consumes the entire query:

| Mode | Target span |
|---|---|
| `Global` (NW) | The entire target |
| `Prefix` (SHW) | A prefix; the trailing flank is free |
| `Infix` (HW) | A substring; both flanks are free |

| Task | Result on a match |
|---|---|
| `Distance` | Distance and all reported optimal ends |
| `Locations` | Also the corresponding starts |
| `Path` | Also one optimal edit path for the first reported span |

`AlignConfig::default()` selects `Global`, `Distance`, no bound, and no extra
equalities. `max_distance: Some(k)` accepts a distance of `k` or less. If no
alignment satisfies the bound, `align` returns `Ok`, but `edit_distance()`
returns `None`, the ends are empty, and starts and path are absent. This
differs from a zero distance and from `AlignError::InputTooLong`, which
means required sizes or band arithmetic exceed platform limits.

All 256 byte values are valid, including NUL and bytes that are not valid UTF-8.
Pass a slice of `EqualityPair` values as `additional_equalities` to treat
different bytes as matches. For example, `EqualityPair::new(b'A', b'a')`
adds equality in both directions. Every byte equals itself; pairs are not
transitive, so adding `A = B` and `B = C` does not add `A = C`.

## Spans, paths, and CIGAR

`locations()` yields half-open target byte ranges in increasing end order;
`start_locations()` gives the starts alone. `end_offsets()` yields exclusive
end offsets even with the `Distance` task. For Edlib compatibility,
`end_locations()` yields inclusive ends, where `-1` means end offset zero.
Global and prefix spans start at zero. For each reported infix end,
`locations()` uses the earliest optimal start.

`path()` covers only the first reported span. `Insert` consumes a query byte,
`Delete` consumes a target byte, and `Match` and `Mismatch` consume one byte
from each. Extended CIGAR uses `=`, `X`, `I`, and `D`; standard CIGAR merges
adjacent matches and mismatches into `M` runs. Neither includes free target
flanks. The result owns its buffers and is `Send + Sync`; accessors borrow
without copying, and `cigar()` returns a formatter that writes directly to
its destination.

## Compatibility and edge cases

- Like `edlib-rs`, the API takes a borrowed configuration. Unlike the bindings,
  it uses `usize` for distances and offsets, `isize` for inclusive ends, and
  exposes no native types or C ABI. It does not interpret Unicode or support
  affine gap penalties.
- In `Prefix` and `Infix` modes, when the target is nonempty and the query
  length is a positive multiple of 64, the optimized fork omits the
  boundary-zero end even when it is optimal. This crate preserves that
  convention for every task without changing the optimal distance.
- Empty inputs still honor the task and bound. An empty global query deletes the
  entire target, an empty prefix query matches only `0..0`, and an empty infix
  query matches every boundary. An empty target consumes the query with
  insertions. A computed empty path is `Some(&[])` and formats as an empty
  CIGAR, rather than returning `None`.
- Paths are deterministic. Backward traceback tries `Insert`, then `Delete`,
  then a diagonal step. Hirschberg tries interior split rows in increasing
  order before the boundary rows. Different trace memory thresholds can
  lead to different, equally optimal paths; identical CIGAR across
  implementations is not guaranteed.

## Algorithm and memory

The distance kernel uses Myers bit vectors with Ukkonen band pruning. Queries
of up to 64 bytes use a single 64-bit word. Without an explicit bound, the
search starts at 64 or the largest feasible distance if smaller, and doubles
the bound as needed up to that limit. The kernel reads inputs without copying
them, including when traversing in reverse. Infix start recovery runs a
reverse prefix search per reported end, so many ties can make `Locations` and
`Path` more expensive than `Distance`.

Path reconstruction starts from the query profile built for the distance search.
Short path problems keep packed traceback columns. Once the optimal distance is
known, a single-word band cannot lose its only word, so queries of up to 64
bytes record their columns without band updates. Traceback tests insertion
steps directly from the stored vertical deltas and shares score decoding between
neighboring cells. Larger problems use banded Hirschberg reconstruction and
release split workspace before recursing.
Working memory is linear in the query length apart from a bounded trace buffer
and the recursion stack. Returned locations and edit operations take additional
space proportional to their respective output sizes.

## Build, test, and docs

From this crate's root, build the library, run the tests, and generate API
docs with `CC` and `CXX` set to fail if invoked:

```sh
CC=false CXX=false cargo build --locked
CC=false CXX=false cargo test --locked
CC=false CXX=false cargo test --locked --release
CC=false CXX=false cargo doc --locked --no-deps
```

Open `target/doc/rust_edlib/index.html` to browse the generated API docs.

## CI and upstream benchmark

[GitHub Actions](.github/workflows/ci.yml) tests debug and release builds on Rust
1.85 with the C and C++ compilers disabled. A separate release job checks scores
against [Martin Šošić's Edlib 1.3.9](https://github.com/Martinsos/edlib/tree/f8afceb49ab0095c852e0b8b488ae2c88e566afd)
and reports median time per call on the same runner. It covers short and long
global searches and bounded prefix hits and misses. The global Rust score-only
call does not collect ends; Edlib's distance task always does. Score mismatches
fail the job; timing differences do not.

To run the comparison locally with a C++ compiler:

```sh
git clone --depth 1 --branch python-1.3.9 \
  https://github.com/Martinsos/edlib.git target/edlib-upstream
EDLIB_SOURCE="$PWD/target/edlib-upstream" \
  cargo run --locked --release --manifest-path tools/edlib-bench/Cargo.toml
```

The benchmark is a separate Cargo workspace; it adds no native dependency to
the library.

## License and references

MIT; see [LICENSE](LICENSE). The banding and path reconstruction algorithms
are translated from [Martin Šošić's Edlib](https://github.com/Martinsos/edlib),
with additional optimizations introduced.