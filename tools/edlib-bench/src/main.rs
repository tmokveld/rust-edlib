use rust_edlib::{
    AlignConfig, AlignmentMode, AlignmentTask, EditOperation, align, distance, distance_at_most,
    infix_distance, infix_distance_at_most, prefix_distance, prefix_distance_at_most,
    suffix_distance, suffix_distance_at_most,
};
use std::{
    ffi::{c_char, c_int},
    fmt::Display,
    hint::black_box,
    time::Instant,
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

fn upstream(
    query: &[u8], target: &[u8], k: c_int, mode: AlignmentMode, path: bool, check_path: bool,
) -> Option<usize> {
    // SAFETY: Both slices are readable for the supplied lengths, and the shim
    // consumes them synchronously without writing to or retaining them.
    let score = unsafe {
        upstream_align(
            query.as_ptr().cast(),
            query.len().try_into().unwrap(),
            target.as_ptr().cast(),
            target.len().try_into().unwrap(),
            k,
            match mode {
                AlignmentMode::Global => 0, // EDLIB_MODE_NW
                AlignmentMode::Prefix => 1, // EDLIB_MODE_SHW
                AlignmentMode::Infix => 2,  // EDLIB_MODE_HW
            },
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

fn sample(iterations: u32, run: &impl Fn() -> Option<usize>) -> f64 {
    let start = Instant::now();
    for _ in 0..iterations {
        black_box(run());
    }
    start.elapsed().as_secs_f64() * 1e9 / f64::from(iterations)
}

fn compare_calls(
    (rust_iterations, cpp_iterations): (u32, u32),
    rust: impl Fn() -> Option<usize>,
    cpp: impl Fn() -> Option<usize>,
) -> (f64, f64) {
    let mut rust_ns = [0.; 5];
    let mut cpp_ns = [0.; 5];
    for (i, (r, c)) in rust_ns.iter_mut().zip(&mut cpp_ns).enumerate() {
        if i % 2 == 0 {
            *r = sample(rust_iterations, &rust);
            *c = sample(cpp_iterations, &cpp);
        } else {
            *c = sample(cpp_iterations, &cpp);
            *r = sample(rust_iterations, &rust);
        }
    }
    rust_ns.sort_unstable_by(f64::total_cmp);
    cpp_ns.sort_unstable_by(f64::total_cmp);
    (rust_ns[2], cpp_ns[2])
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

fn compare_short(len: usize, edits: usize, mode: AlignmentMode, limit: Option<usize>, iterations: u32) {
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
    compare("Substitutions", &query, &target, config, expected, iterations);
}

fn compare_anchored_scores() {
    for (len, edits, iterations) in [(63, 8, 1000), (64, 8, 1000), (65, 8, 1000), (512, 24, 200)] {
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
            let k = limit.map_or(-1, |k| k.try_into().unwrap());
            let expected = limit.is_none_or(|k| edits <= k).then_some(edits);
            let prefix = || {
                let (q, t) = (black_box(query.as_slice()), black_box(prefix_target.as_slice()));
                match limit {
                    None => Some(prefix_distance(q, t).unwrap()),
                    Some(k) => prefix_distance_at_most(q, t, k).unwrap(),
                }
            };
            let suffix = || {
                let (q, t) = (black_box(query.as_slice()), black_box(suffix_target.as_slice()));
                match limit {
                    None => Some(suffix_distance(q, t).unwrap()),
                    Some(k) => suffix_distance_at_most(q, t, k).unwrap(),
                }
            };
            let edlib = |q: &[u8], t: &[u8]| {
                upstream(black_box(q), black_box(t), k, AlignmentMode::Prefix, false, false)
            };
            let cpp_prefix = || edlib(&query, &prefix_target);
            let cpp_suffix = || edlib(&reversed_query, &reversed_target);
            assert_eq!(prefix(), expected, "Rust prefix, len={len}, k={k}");
            assert_eq!(cpp_prefix(), expected, "Edlib prefix, len={len}, k={k}");
            assert_eq!(suffix(), expected, "Rust suffix, len={len}, k={k}");
            assert_eq!(cpp_suffix(), expected, "Edlib reversed prefix, len={len}, k={k}");
            let report = |name, timings| print_row(
                "Substitutions + free flank", name, "Score only", len, len + flank,
                limit, expected, timings,
            );
            report("Prefix", compare_calls((iterations, iterations), prefix, cpp_prefix));
            report("Suffix", compare_calls((iterations, iterations), suffix, cpp_suffix));
        }
    }
}

fn compare_infix_score(
    input: &str,
    query: &[u8],
    target: &[u8],
    distance: usize,
    limit: Option<usize>,
    iterations: (u32, u32),
) {
    let k = limit.map_or(-1, |k| k.try_into().unwrap());
    let expected = limit.is_none_or(|k| distance <= k).then_some(distance);
    let rust = || {
        let (q, t) = (black_box(query), black_box(target));
        match limit {
            None => Some(infix_distance(q, t).unwrap()),
            Some(k) => infix_distance_at_most(q, t, k).unwrap(),
        }
    };
    let cpp = || upstream(black_box(query), black_box(target), k, AlignmentMode::Infix, false, false);
    assert_eq!(rust(), expected, "Rust infix: {input}, k={k}");
    assert_eq!(cpp(), expected, "Edlib HW: {input}, k={k}");
    print_row(
        input, "Infix", "Score only", query.len(), target.len(), limit, expected,
        compare_calls(iterations, rust, cpp),
    );
}

fn compare_infix_scores() {
    for (len, edits, iterations) in [(63, 8, 500), (64, 8, 500), (65, 8, 500), (512, 24, 100)] {
        let (query, core) = substitution_pair(len, edits);
        let mut target = vec![b'T'; 3 * len];
        target[len..2 * len].copy_from_slice(&core);
        // No flank byte matches. Only len - edits target bytes can match
        // anywhere, and substituting the T's in the core attains the bound.
        for limit in [None, Some(edits), Some(edits - 1)] {
            compare_infix_score(
                "Substitutions + two free flanks", &query, &target, edits, limit,
                (iterations, iterations),
            );
        }
    }

    let config = AlignConfig { mode: AlignmentMode::Infix, ..Default::default() };
    let (query, mut target) = substitution_pair(65, 1);
    target.resize(query.len() + 2048, b'T');
    target.extend_from_slice(&query);
    let late = "One-edit hit, then exact";
    compare_infix_score(late, &query, &target, 0, Some(1), (100, 100));
    compare(
        late, &query, &target,
        AlignConfig { task: AlignmentTask::Path, max_distance: Some(1), ..config },
        Some(0), 100,
    );

    let target = vec![b'a'; 1_000_000];
    compare_infix_score("No exact match", b"b", &target, 1, Some(0), (1, 1));
    let ties = "1M tied exact ends";
    // Batch the one-column Rust calls; Edlib materializes a million ends per call.
    for limit in [None, Some(0)] {
        compare_infix_score(ties, b"a", &target, 0, limit, (10_000, 1));
    }
    compare(ties, b"a", &target, config, Some(0), 1);
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

fn compare(
    input: impl Display,
    query: &[u8],
    target: &[u8],
    config: AlignConfig<'_>,
    expected: Option<usize>,
    iterations: u32,
) {
    let path = config.task == AlignmentTask::Path;
    let k = config
        .max_distance
        .map_or(-1, |limit| limit.try_into().unwrap());
    let score_only = !path && config.mode == AlignmentMode::Global;
    let rust = |q: &[u8], t: &[u8]| {
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
    let cpp = |q: &[u8], t: &[u8]| upstream(q, t, k, config.mode, path, false);
    if path {
        let result = align(query, target, &config).unwrap();
        assert_eq!(result.edit_distance(), expected, "Rust path: {input}");
        let span = result.locations().unwrap().next().unwrap();
        if config.mode != AlignmentMode::Infix {
            assert_eq!(span.start, 0, "Rust path start: {input}");
        }
        if config.mode == AlignmentMode::Global {
            assert_eq!(span.end, target.len(), "Rust global path end: {input}");
        }
        check_rust_path(
            query,
            &target[span],
            result.path().unwrap(),
            expected.unwrap(),
        );
    } else {
        assert_eq!(rust(query, target), expected, "Rust score: {input}, k={k}");
    }
    assert_eq!(
        upstream(query, target, k, config.mode, path, path),
        expected,
        "Upstream Edlib: {input}, k={k}, path={path}"
    );

    let timings = compare_calls(
        (iterations, iterations),
        || rust(black_box(query), black_box(target)),
        || cpp(black_box(query), black_box(target)),
    );
    let name = match config.mode {
        AlignmentMode::Global => "Global",
        AlignmentMode::Prefix => "Prefix",
        AlignmentMode::Infix => "Infix",
    };
    let task = if path { "Path" } else if score_only { "Score only" } else { "Distance + ends" };
    print_row(input, name, task, query.len(), target.len(), config.max_distance, expected, timings);
}

fn print_row(
    input: impl Display,
    name: &str,
    task: &str,
    query_len: usize,
    target_len: usize,
    limit: Option<usize>,
    expected: Option<usize>,
    (rust_ns, cpp_ns): (f64, f64),
) {
    let bound = limit.map_or_else(|| "unbounded".into(), |limit| limit.to_string());
    let score = expected.map_or_else(|| "no match".into(), |score| score.to_string());
    println!(
        "| {input} | {name} | {task} | {query_len} | {target_len} | {bound} | {score} | {rust_ns:.0} | {cpp_ns:.0} | {:.2}× |",
        cpp_ns / rust_ns
    );
}

fn compare_score_and_path(
    input: impl Display,
    query: &[u8],
    target: &[u8],
    config: AlignConfig<'_>,
    iterations: u32,
) {
    let k = config
        .max_distance
        .map_or(-1, |limit| limit.try_into().unwrap());
    let expected = upstream(query, target, k, config.mode, false, false);
    assert!(expected.is_some(), "Expected a match: {input}, k={k}");
    for task in [AlignmentTask::Distance, AlignmentTask::Path] {
        compare(
            &input,
            query,
            target,
            AlignConfig { task, ..config },
            expected,
            iterations,
        );
    }
}

fn main() {
    use AlignmentMode::{Global, Prefix};

    println!(
        "\n| Input | Search | Rust task | Query bytes | Target bytes | Bound | Score | Rust ns/call | Edlib ns/call | Edlib / Rust |"
    );
    println!("|---|---|---|---:|---:|---:|---:|---:|---:|---:|");
    let suite = std::env::args().nth(1);
    match suite.as_deref() {
        Some("--anchored-score") => return compare_anchored_scores(),
        Some("--infix-score") => return compare_infix_scores(),
        _ => {}
    }
    for (len, edits, iterations) in [(64, 8, 10_000), (512, 24, 1_000)] {
        for k in [edits, edits - 1] {
            compare_short(len, edits, Global, Some(k), iterations);
        }
    }
    for (len, iterations) in [(10_000, 20), (100_000, 2), (1_000_000, 1)] {
        let target = dna_reference(len);
        let query = mutated_query(&target, 1_000);
        // Mixed edits need not add up to the optimal score. Set both cutoffs
        // from an independent exact distance before timing either scorer.
        let d = upstream(&query, &target, -1, Global, false, false).unwrap();
        assert!(d > 0, "the one-below-cutoff case requires a nonzero distance");
        for (k, expected) in [(d, Some(d)), (d - 1, None)] {
            compare(
                "Mixed edits (0.1%)",
                &query,
                &target,
                AlignConfig { max_distance: Some(k), ..Default::default() },
                expected,
                iterations,
            );
        }
    }
    if suite.as_deref() == Some("--bounded-global") {
        return;
    }
    compare_anchored_scores();
    compare_infix_scores();
    compare_short(64, 8, Global, None, 20_000);
    compare_short(4096, 160, Global, None, 200);
    compare_short(512, 24, Prefix, Some(24), 2_000);
    compare_short(512, 24, Prefix, Some(23), 2_000);

    // Batch the short calls; one call per sample is enough for long alignments.
    for (len, mutations_per_million, iterations) in [
        (500, 10_000, 1_000),
        (1_000, 10_000, 400),
        (10_000, 10_000, 20),
        (100_000, 10_000, 1),
        (100_000, 200_000, 1),
        (1_000_000, 100, 1),
        (1_000_000, 1_000, 1),
    ] {
        let target = dna_reference(len);
        let query = mutated_query(&target, mutations_per_million);
        let rate = f64::from(mutations_per_million) / 10_000.0;
        compare_score_and_path(
            format_args!("Mixed edits ({rate}%)"),
            &query,
            &target,
            AlignConfig::default(),
            iterations,
        );
    }

    // Use the same uneven pair to measure charging for the extra target
    // sequence (global) versus skipping it (prefix).
    let target = dna_reference(100_000);
    let query = mutated_query(&target[..25_000], 10_000);
    let input = "Quarter-length query (1% mixed edits)";
    compare_score_and_path(input, &query, &target, AlignConfig::default(), 1);
    let prefix = AlignConfig {
        mode: AlignmentMode::Prefix,
        max_distance: Some(512),
        ..Default::default()
    };
    compare_score_and_path(input, &query, &target, prefix, 10);
}
