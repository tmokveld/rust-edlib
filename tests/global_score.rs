mod common;

use common::quadratic;
use rust_edlib::{AlignmentMode, GlobalProbe, distance, distance_at_most, distance_with_hint};

fn scalar(a: &[u8], b: &[u8]) -> usize {
    quadratic(a, b, AlignmentMode::Global, &[]).distance
}

fn check(a: &[u8], b: &[u8], d: usize) {
    for (a, b) in [(a, b), (b, a)] {
        assert_eq!(distance(a, b), Ok(d));
        for hint in [0, d.saturating_sub(1), d, d + 1, usize::MAX] {
            assert_eq!(distance_with_hint(a, b, hint), Ok(d), "hint={hint}");
        }
        for k in [0, d.saturating_sub(1), d, d + 1, usize::MAX] {
            assert_eq!(distance_at_most(a, b, k), Ok((d <= k).then_some(d)), "k={k}");
        }
    }
}

#[test]
fn short_optimum_can_cross_the_middle_at_either_query_edge() {
    for n in [2, 31, 32, 63, 64] {
        let query: Vec<_> = (0..n).map(|i| (199 * i) as u8).collect();
        for extra in [n + 2, n + 3] {
            let pad = vec![255; extra];
            for mut target in [[&pad[..], &query].concat(), [&query[..], &pad].concat()] {
                let pos = target.iter().position(|&b| b == query[n / 2]).unwrap();
                target[pos] = 254;
                check(&query, &target, scalar(&query, &target));
            }
        }
    }
}

#[test]
fn tight_cutoffs_reject_short_misses_but_accept_late_matches() {
    let query = [0xff; 64];
    let miss = [0; 64];
    for k in [0, 1, 3, 63] {
        assert_eq!(distance_at_most(&query, &miss, k), Ok(None));
    }
    let mut target = vec![0; 257];
    target.extend(query);
    check(&query, &target, 257);
    check(&query[..1], &target, target.len() - 1);
    check(&query[..1], &miss, miss.len());
    check(&[], &target, target.len());
}

#[test]
fn small_hints_handle_indel_excursions_at_both_edges() {
    let a: Vec<u8> = (0..265).map(|i| ((37 * i + i / 7) % 251) as u8).collect();
    for n in [65, 72, 127, 193, 265] {
        let a = &a[..n];
        for gap in [0, 1, 7, 28, 55, 56] {
            // These excursions use the extreme diagonals of the sliding window,
            // including a partially shifted last tile.
            let deleted = (56 - gap) / 2;
            let inserted = gap + deleted;
            for near_start in [false, true] {
                let mut b = a.to_vec();
                if near_start {
                    b.drain(1..1 + deleted);
                    b.splice(b.len() - 1..b.len() - 1, std::iter::repeat_n(255, inserted));
                } else {
                    b.drain(n - 1 - deleted..n - 1);
                    b.splice(1..1, std::iter::repeat_n(255, inserted));
                }
                check(a, &b, scalar(a, &b));
            }
        }
    }
}

#[test]
fn underestimated_hint_still_finds_the_exact_score() {
    let a: Vec<u8> = (0..265).map(|i| ((37 * i + i / 7) % 251) as u8).collect();
    // The optimum leaves the sliding window, which would overestimate the score.
    let mut b = [&a[..36], &[255; 48], &a[36..198], &a[246..]].concat();
    b[..32].fill(254);
    b[240..].fill(253);
    check(&a, &b, scalar(&a, &b));

    // Below-gap hints start at the mandatory gap, not at gap plus the hint.
    // Finding this optimum also requires growing zero extra slack after a miss.
    let target = [&[255; 129], &a[..261], &[254; 4]].concat();
    check(&a, &target, 133);
    check(&a, &a, 0);
}

#[test]
fn trimmed_global_scores_keep_the_ordered_gap_lower_bound() {
    let head = [0; 65];
    let tail = [255; 97];
    let query = [&head[..], b"ABC", &tail].concat();
    let target = [&head[..], b"_ABC_", &tail].concat();
    check(&query, &target, 2);
    // The multiset contains all query bytes but they are in the wrong order.
    let target = [&head[..], b"_CBA_", &tail].concat();
    assert_eq!(distance_at_most(&query, &target, 2), Ok(None));
    assert_eq!(distance_at_most(&target, &query, 2), Ok(None));
    check(&query, &target, scalar(&query, &target));
    check(&query, &query, 0);
    check(b"", b"ABC", 3);
}

#[test]
fn disjoint_mismatch_proof_gives_exact_scores_but_overlap_is_inconclusive() {
    let mut query = [b'x'; 384];
    query[128..192].fill(b'y');
    let mut target = query;
    target[40..64].fill(b'y');
    target[320..344].fill(b'y');
    let probe = GlobalProbe::new(&query, &target);
    assert_eq!(probe.diagonal_score(), Some(48));
    assert_eq!(probe.distance_at_most(47), Ok(None));
    assert_eq!(probe.distance_at_most(48), Ok(Some(48)));
    check(&query, &target, 48);

    let query: Vec<_> = b"ab".iter().copied().cycle().take(384).collect();
    let target: Vec<_> = b"ba".iter().copied().cycle().take(384).collect();
    let probe = GlobalProbe::new(&query, &target);
    assert_eq!(probe.diagonal_score(), None);
    assert_eq!(probe.distance_at_most(1), Ok(None));
    assert_eq!(probe.distance_at_most(2), Ok(Some(2)));
    assert_eq!(probe.distance_at_most(64), Ok(Some(2)));
    check(&query, &target, 2);
    assert_eq!(GlobalProbe::new(b"x", b"xx").diagonal_score(), None);
    assert_eq!(GlobalProbe::new(b"xx", b"xx").diagonal_score(), Some(0));
}
