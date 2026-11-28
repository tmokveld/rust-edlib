mod common;

use common::{MODES, binary_strings, equal, limits, quadratic};
use rust_edlib::{
    AlignConfig, AlignmentMode, AlignmentTask, EditOperation, EqualityPair, GlobalProbe, Run,
    RunAlignment, align, align_runs,
};

fn replay(query: &[u8], target: &[u8], pairs: &[EqualityPair], alignment: &RunAlignment) {
    let mut query = query;
    let mut target = &target[alignment.target.clone()];
    let mut score = 0;
    let mut previous = None;
    for &Run { op, len } in &alignment.runs {
        assert!(len > 0);
        assert_ne!(previous.replace(op), Some(op), "runs must be coalesced");
        match op {
            EditOperation::Match | EditOperation::Mismatch => {
                for (&a, &b) in query[..len].iter().zip(&target[..len]) {
                    assert_eq!(equal(a, b, pairs), op == EditOperation::Match);
                }
                query = &query[len..];
                target = &target[len..];
            }
            EditOperation::Insert => query = &query[len..],
            EditOperation::Delete => target = &target[len..],
        }
        if op != EditOperation::Match {
            score += len;
        }
    }
    assert!(query.is_empty() && target.is_empty());
    assert_eq!(alignment.score, score);
}

fn check(query: &[u8], target: &[u8], pairs: &[EqualityPair]) {
    for mode in MODES {
        let expected = quadratic(query, target, mode, pairs);
        for bound in limits(expected.distance) {
            let result = align_runs(query, target, mode, bound, pairs).unwrap();
            assert_eq!(result.as_ref().map(|hit| hit.score),
                bound.is_none_or(|k| expected.distance <= k).then_some(expected.distance),
                "mode={mode:?}, bound={bound:?}, query={query:?}, target={target:?}");
            let Some(alignment) = result else { continue };
            match mode {
                AlignmentMode::Global => assert_eq!(alignment.target, 0..target.len()),
                AlignmentMode::Prefix => assert_eq!(alignment.target.start, 0),
                AlignmentMode::Infix => {}
                AlignmentMode::Suffix => {
                    assert_eq!(alignment.target, expected.starts[0]..target.len());
                }
            }
            replay(query, target, pairs, &alignment);
        }
    }
}

#[test]
fn one_optimal_path_for_all_modes_and_bounds() {
    let inputs = binary_strings(3);
    for query in &inputs {
        for target in &inputs {
            check(query, target, &[]);
        }
    }
    // Byte equalities are symmetric, but they do not chain across pairs.
    let pairs = [EqualityPair::new(0, 255), EqualityPair::new(255, 128)];
    for len in [1, 63, 64, 65, 127, 128, 129] {
        let query: Vec<_> = [0, 255, 0, 1].into_iter().cycle().take(len).collect();
        let mut target = vec![2, 2, 2];
        target.extend(query.iter().map(|&b| if b == 0 { 255 } else { b }));
        target[len / 2 + 3] = 128;
        target.extend([2, 2]);
        check(&query, &target, &pairs);
    }
}

#[test]
fn handles_empty_spans_and_tied_ends() {
    let query = [b'a'; 64];
    let target = b"bbbb";
    for mode in [AlignmentMode::Prefix, AlignmentMode::Infix] {
        assert_eq!(align_runs(&query, target, mode, Some(usize::MAX), &[]).unwrap(),
            Some(RunAlignment {
                score: query.len(), target: 0..0,
                runs: vec![Run { op: EditOperation::Insert, len: query.len() }],
            }));
    }
    let suffix = align_runs(&query, target, AlignmentMode::Suffix, None, &[]).unwrap().unwrap();
    assert_eq!(suffix.target, 0..target.len()); // Longest of five tied suffixes.
    replay(&query, target, &[], &suffix);
    assert_eq!(align_runs(b"", target, AlignmentMode::Suffix, Some(0), &[]).unwrap(),
        Some(RunAlignment { score: 0, target: target.len()..target.len(), runs: vec![] }));

    let target = [b'a'; 4096];
    let hit = align_runs(&query, &target, AlignmentMode::Infix, None, &[]).unwrap().unwrap();
    assert_eq!(hit, RunAlignment {
        score: 0, target: 0..query.len(),
        runs: vec![Run { op: EditOperation::Match, len: query.len() }],
    });
}

#[test]
fn zero_cutoff_infix_finds_overlapping_binary_matches() {
    use rust_edlib::infix_distance_at_most;

    for len in [7, 8, 9, 25] {
        let query: Vec<_> = [0, 255].into_iter().cycle().take(len - 1).chain([128]).collect();
        for repeats in [1, 64] {
            let mut target = [0, 255].repeat(repeats);
            let start = target.len();
            // The unique hit overlaps failed candidates and ends at the last byte.
            target.extend_from_slice(&query);
            assert_eq!(infix_distance_at_most(&query, &target, 0), Ok(Some(0)));
            let hit = align_runs(&query, &target, AlignmentMode::Infix, Some(0), &[]).unwrap().unwrap();
            assert_eq!((hit.score, &hit.target), (0, &(start..target.len())));
            assert_eq!(hit.runs, [Run { op: EditOperation::Match, len }]);

            *target.last_mut().unwrap() = 129;
            assert_eq!(infix_distance_at_most(&query, &target, 0), Ok(None));
            assert_eq!(align_runs(&query, &target, AlignmentMode::Infix, Some(0), &[]), Ok(None));
            let pairs = [EqualityPair::new(128, 129)];
            assert_eq!(align_runs(&query, &target, AlignmentMode::Infix, Some(0), &pairs), Ok(Some(hit)));
        }
    }
}

#[test]
fn coalesces_a_run_across_checkpoint_and_traceback_leaves() {
    let query = vec![b'a'; 8192];
    let mut target = query.clone();
    target[0] = b'b';
    target[query.len() - 1] = b'b';
    let hit = align_runs(&query, &target, AlignmentMode::Global, Some(2), &[]).unwrap().unwrap();
    assert_eq!(hit, RunAlignment {
        score: 2, target: 0..target.len(),
        runs: vec![
            Run { op: EditOperation::Mismatch, len: 1 },
            Run { op: EditOperation::Match, len: query.len() - 2 },
            Run { op: EditOperation::Mismatch, len: 1 },
        ],
    });
}

#[test]
fn skewed_trimmed_runs_hit_the_gap_and_retry_a_failed_embedding() {
    let head = [0; 65];
    let tail = [255; 65];
    let query = [&head[..], b"ABC", &tail].concat();
    let middle = [b'_'; 33];
    let target = [&head[..], b"_A", &middle, b"BC_", &tail].concat();
    let miss = [&head[..], b"_A", &middle, b"BX_", &tail].concat();
    let gap = target.len() - query.len();
    for (query, target, miss_query, miss_target, indel) in [
        (&query, &target, &query, &miss, EditOperation::Delete),
        (&target, &query, &miss, &query, EditOperation::Insert),
    ] {
        for bound in [Some(gap), None] {
            let hit = align_runs(query, target, AlignmentMode::Global, bound, &[]).unwrap().unwrap();
            assert_eq!((hit.score, hit.target.clone()), (gap, 0..target.len()));
            replay(query, target, &[], &hit);
            assert_eq!(hit.runs.iter().filter(|r| r.op == indel).map(|r| r.len).sum::<usize>(), gap);
        }
        // Replacing C by X defeats the embedding in both directions.
        assert_eq!(align_runs(miss_query, miss_target, AlignmentMode::Global, Some(gap), &[]), Ok(None));
        for bound in [Some(gap + 1), None] {
            let hit = align_runs(miss_query, miss_target, AlignmentMode::Global, bound, &[]).unwrap().unwrap();
            assert_eq!(hit.score, gap + 1);
            replay(miss_query, miss_target, &[], &hit);
        }
    }
}

#[test]
fn literal_shortcuts_preserve_legacy_ties_and_do_not_use_custom_equalities() {
    use EditOperation::{Delete, Match};

    let query = b"ba";
    let target = b"bbaa";
    let config = AlignConfig { task: AlignmentTask::Path, ..Default::default() };
    let legacy = align(query, target, &config).unwrap();
    assert_eq!(legacy.path(), Some(&[Match, Delete, Match, Delete][..]));
    let hit = align_runs(query, target, AlignmentMode::Global, Some(2), &[]).unwrap().unwrap();
    assert_eq!(hit.score, 2);
    replay(query, target, &[], &hit);

    let query = [b'A'; 65];
    let target = [&[b'_'; 33][..], &[b'a'; 65], b"_"].concat();
    let gap = target.len() - query.len();
    let pairs = [EqualityPair::new(b'A', b'a')];
    for (query, target) in [(&query[..], &target[..]), (&target[..], &query[..])] {
        assert_eq!(align_runs(query, target, AlignmentMode::Global, Some(gap), &[]), Ok(None));
        let hit = align_runs(query, target, AlignmentMode::Global, Some(gap), &pairs).unwrap().unwrap();
        assert_eq!(hit.score, gap);
        replay(query, target, &pairs, &hit);
    }
    let hit = align_runs(&[b'A'; 384], &[b'a'; 384], AlignmentMode::Global, Some(32), &pairs)
        .unwrap().unwrap();
    assert_eq!(hit.runs, [Run { op: Match, len: 384 }]);
}

#[test]
fn cached_diagonal_or_myers_builds_one_full_span() {
    let query = [b'a'; 384];
    let mut target = query;
    target[20..40].fill(b'b');
    target[340..360].fill(b'b');
    let probe = GlobalProbe::new(&query, &target);
    assert_eq!(probe.diagonal_score(), Some(40));
    assert_eq!(probe.align_runs(Some(39)), Ok(None));
    let hit = probe.align_runs(None).unwrap().unwrap();
    assert_eq!(hit.target, 0..target.len());
    assert_eq!(hit.runs.len(), 5);
    replay(&query, &target, &[], &hit);
    assert_eq!(align_runs(&query, &target, AlignmentMode::Global, Some(40), &[]), Ok(Some(hit)));

    // Pin matching ends around a shifted, rather than substituted, core.
    let query = [&b"start"[..], &b"ab".repeat(192), b"end"].concat();
    let target = [&b"start"[..], &b"ba".repeat(192), b"end"].concat();
    let probe = GlobalProbe::new(&query, &target);
    assert_eq!(probe.diagonal_score(), None);
    assert_eq!(probe.align_runs(Some(1)), Ok(None));
    let hit = probe.align_runs(Some(2)).unwrap().unwrap();
    assert_eq!((hit.score, hit.target.clone()), (2, 0..target.len()));
    replay(&query, &target, &[], &hit);
}
