use rust_edlib::{AlignConfig, AlignmentMode, AlignmentTask, EditOperation, align, distance};
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
        path: bool,
        check_path: bool,
    ) -> c_int;
}

fn upstream(query: &[u8], target: &[u8], k: c_int, path: bool, check_path: bool) -> Option<usize> {
    // SAFETY: Both slices are readable for the supplied lengths, and the shim
    // consumes them synchronously without writing to or retaining them.
    let score = unsafe {
        upstream_align(
            query.as_ptr().cast(),
            query.len().try_into().unwrap(),
            target.as_ptr().cast(),
            target.len().try_into().unwrap(),
            k,
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

fn sample(
    query: &[u8],
    target: &[u8],
    iterations: u32,
    run: impl Fn(&[u8], &[u8]) -> Option<usize>,
) -> f64 {
    let start = Instant::now();
    for _ in 0..iterations {
        black_box(run(black_box(query), black_box(target)));
    }
    start.elapsed().as_secs_f64() * 1e9 / f64::from(iterations)
}

fn next_random(state: &mut u32) -> u32 {
    *state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
    *state
}

fn compare_short(len: usize, edits: usize, prefix_limit: Option<usize>, iterations: u32) {
    let mut state = 1_u32;
    let query: Vec<u8> = (0..len)
        .map(|_| b"abcd"[(next_random(&mut state) >> 30) as usize])
        .collect();
    let mut target = query.clone();
    // These bytes cannot match any query byte, so the planted substitution
    // count is also a lower bound. The prefix tail is free to skip.
    for i in 0..edits {
        target[i * len / edits] = b'x';
    }
    if prefix_limit.is_some() {
        target.resize(len + len / 4, b'z');
    }

    let config = AlignConfig {
        mode: if prefix_limit.is_some() {
            AlignmentMode::Prefix
        } else {
            AlignmentMode::Global
        },
        max_distance: prefix_limit,
        ..Default::default()
    };
    let expected = prefix_limit
        .is_none_or(|limit| edits <= limit)
        .then_some(edits);
    compare("Substitutions", &query, &target, config, expected, iterations);
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
            Some(distance(q, t).unwrap())
        } else {
            let result = align(q, t, &config).unwrap();
            black_box(&result).edit_distance()
        }
    };
    let cpp = |q: &[u8], t: &[u8]| upstream(q, t, k, path, false);
    if path {
        let result = align(query, target, &config).unwrap();
        assert_eq!(result.edit_distance(), expected, "Rust path: {input}");
        let span = result.locations().unwrap().next().unwrap();
        assert_eq!(span.start, 0, "Rust path start: {input}");
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
        upstream(query, target, k, path, path),
        expected,
        "Upstream Edlib: {input}, k={k}, path={path}"
    );

    let mut rust_ns = [0.; 5];
    let mut cpp_ns = [0.; 5];
    for (i, (r, c)) in rust_ns.iter_mut().zip(&mut cpp_ns).enumerate() {
        if i % 2 == 0 {
            *r = sample(query, target, iterations, rust);
            *c = sample(query, target, iterations, cpp);
        } else {
            *c = sample(query, target, iterations, cpp);
            *r = sample(query, target, iterations, rust);
        }
    }
    rust_ns.sort_unstable_by(f64::total_cmp);
    cpp_ns.sort_unstable_by(f64::total_cmp);
    let (rust_ns, cpp_ns) = (rust_ns[2], cpp_ns[2]);
    let name = if config.mode == AlignmentMode::Prefix {
        "Prefix"
    } else {
        "Global"
    };
    let task = if path { "Path" } else { "Score" };
    let bound = config
        .max_distance
        .map_or_else(|| "unbounded".into(), |limit| limit.to_string());
    let score = expected.map_or_else(|| "no match".into(), |score| score.to_string());
    println!(
        "| {input} | {name} | {task} | {} | {} | {bound} | {score} | {rust_ns:.0} | {cpp_ns:.0} | {:.2}× |",
        query.len(),
        target.len(),
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
    let expected = upstream(query, target, k, false, false);
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
    println!(
        "\n| Input | Search | Task | Query bytes | Target bytes | Bound | Score | Rust ns/call | Edlib ns/call | Edlib / Rust |"
    );
    println!("|---|---|---|---:|---:|---:|---:|---:|---:|---:|");
    compare_short(64, 8, None, 20_000);
    compare_short(4096, 160, None, 200);
    compare_short(512, 24, Some(24), 2_000);
    compare_short(512, 24, Some(23), 2_000);

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
