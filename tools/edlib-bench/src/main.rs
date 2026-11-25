use rust_edlib::{AlignConfig, AlignmentMode, align, distance};
use std::{
    ffi::{c_char, c_int},
    hint::black_box,
    time::Instant,
};

unsafe extern "C" {
    fn upstream_distance(
        query: *const c_char,
        query_len: c_int,
        target: *const c_char,
        target_len: c_int,
        k: c_int,
    ) -> c_int;
}

fn upstream(query: &[u8], target: &[u8], k: c_int) -> Option<usize> {
    // SAFETY: Both slices are readable for the supplied lengths, and the shim
    // consumes them synchronously without writing to or retaining them.
    let score = unsafe {
        upstream_distance(
            query.as_ptr().cast(),
            query.len().try_into().unwrap(),
            target.as_ptr().cast(),
            target.len().try_into().unwrap(),
            k,
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

fn compare(len: usize, edits: usize, prefix_limit: Option<usize>, iterations: u32) {
    let mut state = 1_u32;
    let query: Vec<u8> = (0..len)
        .map(|_| {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            b"abcd"[(state >> 30) as usize]
        })
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
        mode: AlignmentMode::Prefix,
        max_distance: prefix_limit,
        ..Default::default()
    };
    let k = prefix_limit.map_or(-1, |limit| limit.try_into().unwrap());
    let rust = |q: &[u8], t: &[u8]| match prefix_limit {
        None => Some(distance(q, t).unwrap()),
        Some(_) => align(q, t, &config).unwrap().edit_distance(),
    };
    let cpp = |q: &[u8], t: &[u8]| upstream(q, t, k);
    let expected = prefix_limit
        .is_none_or(|limit| edits <= limit)
        .then_some(edits);
    assert_eq!(
        rust(&query, &target),
        expected,
        "Rust: {len}, {prefix_limit:?}"
    );
    assert_eq!(
        cpp(&query, &target),
        expected,
        "Edlib: {len}, {prefix_limit:?}"
    );

    let mut rust_ns = [0.; 5];
    let mut cpp_ns = [0.; 5];
    for (i, (r, c)) in rust_ns.iter_mut().zip(&mut cpp_ns).enumerate() {
        if i % 2 == 0 {
            *r = sample(&query, &target, iterations, rust);
            *c = sample(&query, &target, iterations, cpp);
        } else {
            *c = sample(&query, &target, iterations, cpp);
            *r = sample(&query, &target, iterations, rust);
        }
    }
    rust_ns.sort_unstable_by(f64::total_cmp);
    cpp_ns.sort_unstable_by(f64::total_cmp);
    let (rust_ns, cpp_ns) = (rust_ns[2], cpp_ns[2]);
    let name = if prefix_limit.is_some() {
        "Prefix"
    } else {
        "Global"
    };
    let bound = prefix_limit.map_or_else(|| "unbounded".into(), |limit| limit.to_string());
    let score = expected.map_or_else(|| "no match".into(), |score| score.to_string());
    println!(
        "| {name} | {len} | {} | {bound} | {score} | {rust_ns:.0} | {cpp_ns:.0} | {:.2}× |",
        target.len(),
        cpp_ns / rust_ns
    );
}

fn main() {
    println!(
        "\n| Search | Query bytes | Target bytes | Bound | Score | Rust ns/call | Edlib ns/call | Edlib / Rust |"
    );
    println!("|---|---:|---:|---:|---:|---:|---:|---:|");
    compare(64, 8, None, 20_000);
    compare(4096, 160, None, 200);
    compare(512, 24, Some(24), 2_000);
    compare(512, 24, Some(23), 2_000);
}
