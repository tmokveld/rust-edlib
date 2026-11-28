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
building an edit path.

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
If you have an estimate rather than a cutoff, `distance_with_hint(query, target, hint)`
returns the exact global score even when `hint` is too small.

```rust
use rust_edlib::{
    AlignError, distance_at_most, infix_distance_at_most, prefix_distance_at_most,
    suffix_distance_at_most,
};

fn main() -> Result<(), AlignError> {
    assert_eq!(distance_at_most(b"ACGT", b"AGGGA", 3)?, Some(3));
    assert_eq!(distance_at_most(b"ACGT", b"AGGGA", 2)?, None);

    // ATG is one edit from ACG; each mode ignores a different target flank.
    assert_eq!(prefix_distance_at_most(b"ACG", b"ATGTT", 1)?, Some(1));
    assert_eq!(prefix_distance_at_most(b"ACG", b"ATGTT", 0)?, None);

    assert_eq!(suffix_distance_at_most(b"ACG", b"TTATG", 1)?, Some(1));
    assert_eq!(suffix_distance_at_most(b"ACG", b"TTATG", 0)?, None);

    assert_eq!(infix_distance_at_most(b"ACG", b"TTATGTT", 1)?, Some(1));
    assert_eq!(infix_distance_at_most(b"ACG", b"TTATGTT", 0)?, None);
    Ok(())
}
```

Use `align` for optimal ends, custom byte equalities, or an edit path. Even with
`AlignmentTask::Distance`, infix `align` collects all reported optimal ends.
`infix_distance` needs only the score and stops at the first exact match.
When you also need the range but not a path, `prefix_locate`, `suffix_locate`, and
`infix_locate` return one optimal half-open target span and its score using
literal byte equality. The prefix starts at zero; the suffix ends at the target
length and is the longest optimal suffix. Spans can be empty. The corresponding
`*_locate_at_most` functions return `Ok(None)` on an inclusive cutoff miss.

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

Use `align_runs` when you need one optimal target range and a coalesced path,
not every tied location or a byte per operation:

```rust
use rust_edlib::{AlignError, AlignmentMode, CigarFormat, EditOperation, Run, align_runs};

fn main() -> Result<(), AlignError> {
    let hit = align_runs(b"TAG", b"xxTAG", AlignmentMode::Suffix, Some(0), &[])?.unwrap();
    assert_eq!(hit.target, 2..5);
    assert_eq!(hit.runs, [Run { op: EditOperation::Match, len: 3 }]);
    assert_eq!(hit.cigar(CigarFormat::Extended).to_string(), "3=");
    Ok(())
}
```

## Modes, tasks, and options

Every mode consumes the entire query:

| Mode | Target span |
|---|---|
| `Global` (NW) | The entire target |
| `Prefix` (SHW) | A prefix; the trailing flank is free |
| `Suffix` (Reverse SHW) | A suffix; the leading flank is free |
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

## Spans, paths, and CIGAR

`locations()` yields half-open target byte ranges in increasing end order,
`start_locations()` gives the starts alone. `end_offsets()` yields exclusive
end offsets even with the `Distance` task. For Edlib compatibility,
`end_locations()` yields inclusive ends, where `-1` means end offset zero.
Global and prefix spans start at zero. For each reported infix end,
`locations()` uses the earliest optimal start.
`path()` covers only the first reported span. `Insert` consumes a query byte,
`Delete` consumes a target byte, and `Match` and `Mismatch` consume one byte
from each. Extended CIGAR uses `=`, `X`, `I`, and `D`, standard CIGAR merges
adjacent matches and mismatches into `M` runs.

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

## CI and benchmarks

[GitHub Actions](.github/workflows/ci.yml) tests debug and release builds on Rust
1.85 with the C and C++ compilers disabled.

See the [benchmark guide](tools/edlib-bench/README.md) to compare performance.

## License and references

MIT; see [LICENSE](LICENSE). The banding and path reconstruction algorithms
are translated from [Martin Šošić's Edlib](https://github.com/Martinsos/edlib),
with additional optimizations introduced.