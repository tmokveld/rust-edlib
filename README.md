# rust-edlib

Exact, unit-cost Levenshtein alignment over byte slices in safe Rust.

## Installation

Use Rust 1.85 or later. Add the crate from crates.io:

```sh
cargo add rust-edlib
```

## Examples

Each distance function scores the whole query. `rust_edlib::distance(query, target)`
compares it with the whole target. `prefix_distance` and `suffix_distance`
minimize over target prefixes and suffixes, respectively. `infix_distance`
minimizes over all contiguous target substrings, including the empty one.
All four return `Result<usize, AlignError>` without collecting locations or
building an edit path. The suffix search reads the slices backwards without
copying them.

```rust
use rust_edlib::{AlignError, distance, infix_distance, prefix_distance, suffix_distance};

fn main() -> Result<(), AlignError> {
    let target = b"CAGCAA";
    assert_eq!(distance(b"CAG", target)?, 3);
    assert_eq!(prefix_distance(b"CAG", target)?, 0);
    assert_eq!(suffix_distance(b"CAA", target)?, 0);
    assert_eq!(suffix_distance(b"CAG", target)?, 1); // Closest suffix: CAA.
    assert_eq!(infix_distance(b"AGC", target)?, 0);
    Ok(())
}
```

For score-only searches with an inclusive cutoff `k`, use `distance_at_most`,
`prefix_distance_at_most`, `suffix_distance_at_most`, or `infix_distance_at_most`.
Each returns `Ok(Some(d))` when its minimum score `d` is at most `k`, or `Ok(None)`
otherwise. A sizing error returns `Err(AlignError::InputTooLong)`, not a cutoff miss.

```rust
use rust_edlib::{AlignError, distance_at_most};

fn main() -> Result<(), AlignError> {
    assert_eq!(distance_at_most(b"ACGT", b"AGGGA", 3)?, Some(3));
    assert_eq!(distance_at_most(b"ACGT", b"AGGGA", 2)?, None);
    Ok(())
}
```

Use `align` for optimal ends, custom byte equalities, or an edit path. Even with
`AlignmentTask::Distance`, infix `align` collects all reported optimal ends.
`infix_distance` needs only the score and stops at the first exact match.

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
transitive, so adding `A = G` and `G = C` does not add `A = C`.

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

## Algorithm and memory

The distance functions use banded Myers bit vectors. Score-only searches need
at most linear space in the query length and store no locations. For infix
search, the first exact match settles the score; `align` keeps scanning for
tied ends.

`align` also needs space for the locations and edit operations it returns.
For long alignments, it reconstructs the path with Hirschberg's algorithm
instead of storing the full alignment matrix.

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
1.85 with the C and C++ compilers disabled. A separate release job compares
against [Martin Šošić's Edlib 1.3.9](https://github.com/Martinsos/edlib/tree/f8afceb49ab0095c852e0b8b488ae2c88e566afd).

To run the comparison locally with a C++ compiler:

```sh
git clone --depth 1 --branch python-1.3.9 \
  https://github.com/Martinsos/edlib.git target/edlib-upstream
EDLIB_SOURCE="$PWD/target/edlib-upstream" \
  cargo run --locked --release --manifest-path tools/edlib-bench/Cargo.toml
```

## License and references

MIT; see [LICENSE](LICENSE). The banding and path reconstruction algorithms
are translated from [Martin Šošić's Edlib](https://github.com/Martinsos/edlib),
with additional optimizations introduced.