#[path = "support/report.rs"]
mod report;

use criterion::{BenchmarkGroup, BenchmarkId, Criterion, measurement::WallTime};
use report::Table;
use rust_edlib::{
    AlignConfig, AlignError, AlignmentMode, AlignmentTask, EditOperation, align, distance,
    distance_at_most, infix_distance, infix_distance_at_most, prefix_distance, prefix_distance_at_most,
    suffix_distance, suffix_distance_at_most,
};
use std::{
    cell::LazyCell,
    ffi::{c_char, c_int},
    fmt,
    hint::black_box,
    io,
    path::Path,
    time::Duration,
};

unsafe extern "C" {
    fn upstream_align(
        query: *const c_char,
        query_len: c_int,
        target: *const c_char,
        target_len: c_int,
        k: c_int,
        mode: c_int,
        path: bool,
        check_path: bool,
    ) -> c_int;
}

fn upstream_mode(mode: AlignmentMode) -> (c_int, &'static str) {
    // The baseline crate may expose only the modes shared with upstream.
    #[allow(unreachable_patterns)]
    match mode {
        AlignmentMode::Global => (0, "global"), // EDLIB_MODE_NW
        AlignmentMode::Prefix => (1, "prefix"), // EDLIB_MODE_SHW
        AlignmentMode::Infix => (2, "infix"),   // EDLIB_MODE_HW
        _ => panic!("alignment mode is not supported by upstream Edlib"),
    }
}

fn upstream(
    query: &[u8],
    target: &[u8],
    k: c_int,
    mode: AlignmentMode,
    path: bool,
    check_path: bool,
) -> Option<usize> {
    let (mode, _) = upstream_mode(mode);
    // SAFETY: Both slices are readable for the supplied lengths, and the shim
    // consumes them synchronously without writing to or retaining them.
    let score = unsafe {
        upstream_align(
            query.as_ptr().cast(),
            query.len().try_into().unwrap(),
            target.as_ptr().cast(),
            target.len().try_into().unwrap(),
            k,
            mode,
            path,
            check_path,
        )
    };
    match score {
        -1 => None,
        0.. => Some(score as usize),
        _ => panic!("upstream Edlib failed with code {score}"),
    }
}

struct Suite<'a> {
    name: &'static str,
    group: BenchmarkGroup<'a, WallTime>,
    table: Option<&'a mut Table>,
    rust_only: bool,
}

impl<'a> Suite<'a> {
    fn new(criterion: &'a mut Criterion, name: &'static str, table: Option<&'a mut Table>) -> Self {
        Self {
            name,
            group: criterion.benchmark_group(name),
            table,
            rust_only: std::env::var_os("EDLIB_BENCH_RUST_ONLY").is_some(),
        }
    }
}

struct Outcome {
    cutoff: Option<usize>,
    score: Option<usize>,
}

fn bench_pair<R, E>(
    suite: &mut Suite<'_>,
    case: &Case<'_>,
    prepare: impl FnOnce() -> (Outcome, R, E),
) where
    R: Fn() -> Option<usize>,
    E: Fn() -> Option<usize>,
{
    // Validate once before the first selected implementation is timed. A filter
    // or --list should not run alignments for unselected benchmark IDs.
    let calls = LazyCell::new(prepare);
    let id = case.to_string();
    let dir = suite.table.as_ref().map(|table| table.case_dir(suite.name, &id));
    let rust = measure(&mut suite.group, &id, "rust", dir.as_deref(), || &calls.1);
    if suite.rust_only {
        return;
    }
    let edlib = measure(&mut suite.group, &id, "edlib", dir.as_deref(), || &calls.2);
    if let Some(table) = &mut suite.table {
        if rust.is_some() || edlib.is_some() {
            table
                .write_row(suite.name, case, &calls.0, rust, edlib)
                .expect("writing benchmark comparison table");
        }
    }
}

fn measure<'a, R: Fn() -> Option<usize> + 'a>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    case: &str,
    engine: &str,
    dir: Option<&Path>,
    run: impl Fn() -> &'a R,
) -> Option<f64> {
    let mut pending = dir.map(|dir| report::estimate_path(dir, engine));
    let mut started = None;
    group.bench_function(BenchmarkId::new(case, engine), |b| {
        // A filtered-out benchmark never invokes this callback. Clear only a
        // selected engine's old estimate, once, outside the measured loop.
        if let Some(path) = pending.take() {
            report::clear_estimate(&path);
            started = Some(path);
        }
        b.iter(run());
    });
    started.as_deref().map(report::median)
}

struct Case<'a> {
    input: &'a str,
    mode: &'static str,
    task: &'static str,
    query_len: usize,
    target_len: usize,
    bound_id: String,
}

impl fmt::Display for Case<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}/{}/{}/{}x{}/k={}",
            self.input, self.mode, self.task, self.query_len, self.target_len, self.bound_id,
        )
    }
}

fn benchmark_case<'a>(
    input: &'a str,
    (query, target): (&[u8], &[u8]),
    mode: &'static str,
    task: &'static str,
    bound_id: String,
) -> Case<'a> {
    Case {
        input,
        mode,
        task,
        query_len: query.len(),
        target_len: target.len(),
        bound_id,
    }
}

fn bound_label(limit: Option<usize>) -> String {
    limit.map_or_else(|| "unbounded".into(), |k| k.to_string())
}

fn upstream_bound(limit: Option<usize>) -> c_int {
    limit.map_or(-1, |k| k.try_into().unwrap())
}

fn next_random(state: &mut u32) -> u32 {
    *state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
    *state
}

fn substitution_pair(len: usize, edits: usize) -> (Vec<u8>, Vec<u8>) {
    let mut state = 1_u32;
    let query: Vec<u8> = (0..len)
        .map(|_| b"ACG"[((next_random(&mut state) >> 16) % 3) as usize])
        .collect();
    let mut target = query.clone();
    // T is absent from the query, so every planted substitution costs one edit.
    for i in 0..edits {
        target[i * len / edits] = b'T';
    }
    (query, target)
}

fn compare_short(
    group: &mut Suite<'_>,
    len: usize,
    edits: usize,
    mode: AlignmentMode,
    limit: Option<usize>,
) {
    let (query, mut target) = substitution_pair(len, edits);
    if mode == AlignmentMode::Prefix {
        target.resize(len + len / 4, b'T');
    }

    let config = AlignConfig {
        mode,
        max_distance: limit,
        ..Default::default()
    };
    let expected = limit.is_none_or(|k| edits <= k).then_some(edits);
    compare(group, "substitutions", (&query, &target), config, || expected);
}

fn anchored_scores(c: &mut Criterion, table: Option<&mut Table>) {
    let mut group = Suite::new(c, "anchored-score", table);
    for (len, edits) in [(63, 8), (64, 8), (65, 8), (512, 24)] {
        let (query, mut prefix_target) = substitution_pair(len, edits);
        let flank = len / 4;
        let mut suffix_target = Vec::with_capacity(len + flank);
        suffix_target.resize(flank, b'T');
        suffix_target.extend_from_slice(&prefix_target);
        prefix_target.resize(len + flank, b'T');
        // Edlib has SHW but no suffix mode. Reverse both slices outside timing.
        let reversed_query: Vec<_> = query.iter().rev().copied().collect();
        let reversed_target: Vec<_> = suffix_target.iter().rev().copied().collect();

        // The target has only len - edits bytes from the query alphabet.
        // At least `edits` query bytes must be changed or deleted; the core
        // attains that score regardless of which target end is free.
        for limit in [None, Some(edits), Some(edits - 1)] {
            let prefix = (query.as_slice(), prefix_target.as_slice());
            compare_score(
                &mut group,
                "free-flank",
                prefix,
                prefix,
                ("prefix", AlignmentMode::Prefix, limit),
                edits,
                scorer(prefix_distance, prefix_distance_at_most, limit),
            );
            compare_score(
                &mut group,
                "free-flank",
                (&query, &suffix_target),
                (&reversed_query, &reversed_target),
                ("suffix", AlignmentMode::Prefix, limit),
                edits,
                scorer(suffix_distance, suffix_distance_at_most, limit),
            );
        }
    }
}

fn scorer(
    exact: impl Fn(&[u8], &[u8]) -> Result<usize, AlignError>,
    bounded: impl Fn(&[u8], &[u8], usize) -> Result<Option<usize>, AlignError>,
    limit: Option<usize>,
) -> impl Fn(&[u8], &[u8]) -> Option<usize> {
    move |q, t| match limit {
        None => Some(exact(q, t).unwrap()),
        Some(k) => bounded(q, t, k).unwrap(),
    }
}

fn compare_score(
    group: &mut Suite<'_>,
    input: &str,
    sequences: (&[u8], &[u8]),
    (edlib_query, edlib_target): (&[u8], &[u8]),
    (mode, edlib_mode, limit): (&'static str, AlignmentMode, Option<usize>),
    score: usize,
    score_rust: impl Fn(&[u8], &[u8]) -> Option<usize>,
) {
    let case = benchmark_case(input, sequences, mode, "score", bound_label(limit));
    let (query, target) = sequences;
    bench_pair(group, &case, || {
        let k = upstream_bound(limit);
        let expected = limit.is_none_or(|k| score <= k).then_some(score);
        let rust = move || score_rust(black_box(query), black_box(target));
        let edlib = move || {
            upstream(
                black_box(edlib_query),
                black_box(edlib_target),
                k,
                edlib_mode,
                false,
                false,
            )
        };
        assert_eq!(rust(), expected, "Rust score: {case}");
        assert_eq!(edlib(), expected, "Upstream Edlib: {case}");
        (Outcome { cutoff: limit, score: expected }, rust, edlib)
    });
}

fn compare_infix_score(
    group: &mut Suite<'_>,
    input: &str,
    sequences: (&[u8], &[u8]),
    limit: Option<usize>,
    score: usize,
) {
    compare_score(
        group,
        input,
        sequences,
        sequences,
        ("infix", AlignmentMode::Infix, limit),
        score,
        scorer(infix_distance, infix_distance_at_most, limit),
    );
}

fn infix_scores(c: &mut Criterion, table: Option<&mut Table>) {
    let mut group = Suite::new(c, "infix-score", table);
    for (len, edits) in [(63, 8), (64, 8), (65, 8), (512, 24)] {
        let (query, core) = substitution_pair(len, edits);
        let mut target = vec![b'T'; 3 * len];
        target[len..2 * len].copy_from_slice(&core);
        // No flank byte matches. Only len - edits target bytes can match
        // anywhere, and substituting the T's in the core attains the bound.
        for limit in [None, Some(edits), Some(edits - 1)] {
            let sequences = (query.as_slice(), target.as_slice());
            compare_infix_score(&mut group, "free-flanks", sequences, limit, edits);
        }
    }

    let config = AlignConfig {
        mode: AlignmentMode::Infix,
        ..Default::default()
    };
    let (query, mut target) = substitution_pair(65, 1);
    target.resize(query.len() + 2048, b'T');
    target.extend_from_slice(&query);
    let late = "late-exact";
    let sequences = (query.as_slice(), target.as_slice());
    compare_infix_score(&mut group, late, sequences, Some(1), 0);
    compare(
        &mut group,
        late,
        sequences,
        AlignConfig {
            task: AlignmentTask::Path,
            max_distance: Some(1),
            ..config
        },
        || Some(0),
    );

    let target = vec![b'a'; 1_000_000];
    let miss = (b"b".as_slice(), target.as_slice());
    compare_infix_score(&mut group, "no-exact", miss, Some(0), 1);
    let ties = "tied-exact";
    let sequences = (b"a".as_slice(), target.as_slice());
    for limit in [None, Some(0)] {
        compare_infix_score(&mut group, ties, sequences, limit, 0);
    }
    compare(&mut group, ties, sequences, config, || Some(0));
}

const DNA: &[u8; 4] = b"ACGT";

fn dna_reference(len: usize) -> Vec<u8> {
    let mut state = 1_u32;
    (0..len)
        .map(|_| DNA[(next_random(&mut state) >> 30) as usize])
        .collect()
}

fn mutated_query(source: &[u8], mutations_per_million: u32) -> Vec<u8> {
    let mut mutation_state = 42_u32;
    let mut query = Vec::with_capacity(source.len());
    let cutoff = (u64::from(mutations_per_million) << 32) / 1_000_000;
    for &base in source {
        if u64::from(next_random(&mut mutation_state)) >= cutoff {
            query.push(base);
            continue;
        }
        match (next_random(&mut mutation_state) >> 16) % 3 {
            0 => {
                let original = DNA.iter().position(|&b| b == base).unwrap();
                let offset = ((next_random(&mut mutation_state) >> 16) % 3) as usize + 1;
                query.push(DNA[(original + offset) % DNA.len()]);
            }
            1 => {
                query.push(DNA[(next_random(&mut mutation_state) >> 30) as usize]);
                query.push(base);
            }
            _ => {} // Delete the base from the query.
        }
    }
    query
}

fn check_rust_path(query: &[u8], target: &[u8], path: &[EditOperation], score: usize) {
    let (mut q, mut t, mut cost) = (0, 0, 0);
    for &op in path {
        match op {
            EditOperation::Match | EditOperation::Mismatch => {
                assert_eq!(query[q] == target[t], op == EditOperation::Match);
                q += 1;
                t += 1;
            }
            EditOperation::Insert => q += 1,
            EditOperation::Delete => t += 1,
        }
        cost += usize::from(op != EditOperation::Match);
    }
    assert_eq!((q, t, cost), (query.len(), target.len(), score));
}

fn bench_alignment<'a>(
    group: &mut Suite<'_>,
    case: &Case<'_>,
    (query, target): (&[u8], &[u8]),
    prepare: impl FnOnce() -> (AlignConfig<'a>, Option<usize>),
) {
    bench_pair(group, case, || {
        let (config, expected) = prepare();
        let path = config.task == AlignmentTask::Path;
        let k = upstream_bound(config.max_distance);
        let score_only = !path && config.mode == AlignmentMode::Global;
        let rust = move |q: &[u8], t: &[u8]| {
            if score_only {
                match config.max_distance {
                    Some(k) => distance_at_most(q, t, k).unwrap(),
                    None => Some(distance(q, t).unwrap()),
                }
            } else {
                let result = align(q, t, &config).unwrap();
                black_box(&result).edit_distance()
            }
        };
        if path {
            let result = align(query, target, &config).unwrap();
            assert_eq!(result.edit_distance(), expected, "Rust path: {case}");
            let span = result.locations().unwrap().next().unwrap();
            if matches!(config.mode, AlignmentMode::Global | AlignmentMode::Prefix) {
                assert_eq!(span.start, 0, "Rust path start: {case}");
            }
            if config.mode == AlignmentMode::Global {
                assert_eq!(span.end, target.len(), "Rust anchored path end: {case}");
            }
            check_rust_path(
                query,
                &target[span],
                result.path().unwrap(),
                expected.unwrap(),
            );
        } else {
            assert_eq!(rust(query, target), expected, "Rust score: {case}");
        }
        assert_eq!(
            upstream(query, target, k, config.mode, path, path),
            expected,
            "Upstream Edlib: {case}"
        );

        (
            Outcome { cutoff: config.max_distance, score: expected },
            move || rust(black_box(query), black_box(target)),
            move || upstream(black_box(query), black_box(target), k, config.mode, path, false),
        )
    });
}

fn compare(
    group: &mut Suite<'_>,
    input: &str,
    sequences: (&[u8], &[u8]),
    config: AlignConfig<'_>,
    expected: impl FnOnce() -> Option<usize>,
) {
    let task = if config.task == AlignmentTask::Path {
        "path"
    } else if config.mode == AlignmentMode::Global {
        "score"
    } else {
        "distance-ends"
    };
    let (_, mode) = upstream_mode(config.mode);
    let case = benchmark_case(input, sequences, mode, task, bound_label(config.max_distance));
    bench_alignment(group, &case, sequences, || (config, expected()));
}

fn compare_score_and_path(
    group: &mut Suite<'_>,
    input: &str,
    (query, target): (&[u8], &[u8]),
    config: AlignConfig<'_>,
) {
    let expected = LazyCell::new(|| {
        let k = upstream_bound(config.max_distance);
        let expected = upstream(query, target, k, config.mode, false, false);
        assert!(expected.is_some(), "Expected a match: {input}, k={k}");
        expected
    });
    for task in [AlignmentTask::Distance, AlignmentTask::Path] {
        compare(
            group,
            input,
            (query, target),
            AlignConfig { task, ..config },
            || *expected,
        );
    }
}

fn bounded_global(c: &mut Criterion, table: Option<&mut Table>) {
    let mut group = Suite::new(c, "bounded-global", table);
    let mode = AlignmentMode::Global;
    for (len, edits) in [(64, 8), (512, 24)] {
        for k in [edits, edits - 1] {
            compare_short(&mut group, len, edits, mode, Some(k));
        }
    }
    for len in [10_000, 100_000, 1_000_000] {
        let target = dna_reference(len);
        let query = mutated_query(&target, 1_000);
        let sequences = (query.as_slice(), target.as_slice());
        // Mixed edits need not add up to the optimal score. Set both cutoffs
        // from an independent exact distance before timing either scorer.
        let score = LazyCell::new(|| {
            let d = upstream(&query, &target, -1, mode, false, false).unwrap();
            assert!(d > 0, "the one-below-cutoff case requires a nonzero distance");
            d
        });
        for (bound, miss) in [("score", false), ("score-1", true)] {
            let case = benchmark_case("mixed-edits-0.1%", sequences, "global", "score", bound.into());
            bench_alignment(&mut group, &case, sequences, || {
                let config = AlignConfig {
                    max_distance: Some(*score - usize::from(miss)),
                    ..Default::default()
                };
                (config, (!miss).then_some(*score))
            });
        }
    }
}

fn alignments(c: &mut Criterion, table: Option<&mut Table>) {
    use AlignmentMode::{Global, Prefix};

    let mut group = Suite::new(c, "alignment", table);
    compare_short(&mut group, 64, 8, Global, None);
    compare_short(&mut group, 4096, 160, Global, None);
    compare_short(&mut group, 512, 24, Prefix, Some(24));
    compare_short(&mut group, 512, 24, Prefix, Some(23));

    for (len, mutations_per_million) in [
        (500, 10_000),
        (1_000, 10_000),
        (10_000, 10_000),
        (100_000, 10_000),
        (100_000, 200_000),
        (1_000_000, 100),
        (1_000_000, 1_000),
    ] {
        let target = dna_reference(len);
        let query = mutated_query(&target, mutations_per_million);
        let rate = f64::from(mutations_per_million) / 10_000.0;
        compare_score_and_path(
            &mut group,
            &format!("mixed-edits-{rate}%"),
            (&query, &target),
            AlignConfig::default(),
        );
    }

    // Use the same uneven pair to measure charging for the extra target
    // sequence (global) versus skipping it (prefix).
    let target = dna_reference(100_000);
    let query = mutated_query(&target[..25_000], 10_000);
    let input = "quarter-length-1%";
    let sequences = (query.as_slice(), target.as_slice());
    compare_score_and_path(&mut group, input, sequences, AlignConfig::default());
    let prefix = AlignConfig {
        mode: AlignmentMode::Prefix,
        max_distance: Some(512),
        ..Default::default()
    };
    compare_score_and_path(&mut group, input, sequences, prefix);
}

fn main() -> io::Result<()> {
    let mut criterion = Criterion::default()
        .warm_up_time(Duration::from_millis(100))
        .measurement_time(Duration::from_millis(250))
        .sample_size(30)
        .nresamples(10_000)
        .configure_from_args();
    let mut table = Table::from_env()?;
    bounded_global(&mut criterion, table.as_mut());
    anchored_scores(&mut criterion, table.as_mut());
    infix_scores(&mut criterion, table.as_mut());
    alignments(&mut criterion, table.as_mut());
    criterion.final_summary();
    if let Some(table) = table {
        table.finish()?;
    }
    Ok(())
}
