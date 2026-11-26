mod common;

use common::*;
use rust_edlib::aligner::{
    AlignConfig, Alignment, AlignmentMode, AlignmentTask, CigarFormat, EditOperation, EqualityPair,
    align,
};

fn exercise(query: &[u8], target: &[u8], pairs: &[EqualityPair]) {
    for mode in MODES {
        let expected = quadratic(query, target, mode, pairs);
        for task in TASKS {
            for limit in limits(expected.distance) {
                let config = AlignConfig {
                    mode,
                    task,
                    max_distance: limit,
                    additional_equalities: pairs,
                };
                let result = align(query, target, &config).unwrap_or_else(|error| {
                    panic!("seed={SEED:#x} {config:?} query={query:?} target={target:?}: {error}")
                });
                assert_snapshot(
                    query,
                    target,
                    &config,
                    &expected,
                    &Snapshot::capture(&result),
                );
            }
        }
    }
}

#[test]
fn exhaustive_binary_inputs_all_modes_tasks_and_distance_bounds() {
    let sequences = binary_strings(4);
    for query in &sequences {
        for target in &sequences {
            exercise(query, target, &[]);
        }
    }
}

#[test]
fn score_only_distances_match_scalar_dp() {
    use rust_edlib::{
        distance, prefix_distance, prefix_distance_at_most, suffix_distance, suffix_distance_at_most,
    };
    let check = |q: &[u8], t: &[u8]| {
        let reverse = |s: &[u8]| s.iter().rev().copied().collect::<Vec<_>>();
        let oracle = |a: &[u8], b: &[u8]| quadratic(a, b, AlignmentMode::Prefix, &[]).distance;
        let global = quadratic(q, t, AlignmentMode::Global, &[]).distance;
        assert_eq!(distance(q, t), Ok(global));
        assert_eq!(distance(t, q), Ok(global));
        let prefix = oracle(q, t);
        let suffix = oracle(&reverse(q), &reverse(t));
        assert_eq!(prefix_distance(q, t), Ok(prefix));
        assert_eq!(suffix_distance(q, t), Ok(suffix));
        let scorers = [prefix_distance_at_most, suffix_distance_at_most];
        for (score, d) in scorers.into_iter().zip([prefix, suffix]) {
            for k in limits(d).into_iter().flatten() {
                assert_eq!(score(q, t, k), Ok((d <= k).then_some(d)));
            }
        }
    };
    let words = binary_strings(5);
    for a in &words {
        for b in &words {
            check(a, b);
        }
    }
    let mut rng = Rng::new(SEED);
    for n in BOUNDARY_LENGTHS {
        let query = rng.bytes(n);
        for m in [63, 64, 65, 129, n, 3 * n] {
            check(&query, &rng.bytes(m));
        }
        check(&query, &mutated(&mut rng, &query));
    }
    // An expensive pinned end and a cheap remainder give different retry
    // predictions in opposite directions. Neither can fix an endpoint.
    let query = rng.bytes(513);
    let mut target = query.clone();
    target[1..193].fill(255);
    target.splice(257..257, [0; 193]);
    target[0] ^= 1;
    *target.last_mut().unwrap() ^= 1;
    check(&query, &target);
    check(&target, &query);
}

#[test]
fn score_only_anchored_block_pruning_tracks_the_last_real_row() {
    // Cross the 2048-column strong-pruning boundary and force multiple
    // adaptive bands. Every query A absent from the target needs an edit,
    // even if the target end is free. Substituting all B's attains the bound.
    for len in [4095, 4096, 4097] {
        let query = vec![b'A'; len];
        let mut target = query.clone();
        for b in target.iter_mut().step_by(23) {
            *b = b'B';
        }
        target[len - 1] = b'B';
        let differences = target.iter().filter(|&&b| b == b'B').count();
        assert_eq!(rust_edlib::prefix_distance(&query, &target), Ok(differences));
        assert_eq!(rust_edlib::suffix_distance(&query, &target), Ok(differences));
        for score in [rust_edlib::prefix_distance_at_most, rust_edlib::suffix_distance_at_most] {
            assert_eq!(score(&query, &target, differences - 1), Ok(None));
            assert_eq!(score(&query, &target, differences), Ok(Some(differences)));
        }
    }
}

#[test]
fn exhaustive_nontransitive_and_symmetric_equalities() {
    // 0~1 and 1~2 never imply 0~2, including after reversing a path.
    let pairs = [
        EqualityPair::new(0, 1),
        EqualityPair::new(1, 2),
        EqualityPair::new(1, 0),
        EqualityPair::new(255, 255),
    ];
    let sequences: Vec<Vec<u8>> = (0..=2)
        .flat_map(|len| {
            (0..3usize.pow(len)).map(move |mut digits| {
                (0..len)
                    .map(|_| {
                        let byte = (digits % 3) as u8;
                        digits /= 3;
                        byte
                    })
                    .collect()
            })
        })
        .collect();
    for query in &sequences {
        for target in &sequences {
            exercise(query, target, &pairs);
        }
    }
    assert_eq!(
        align(
            &[0],
            &[2],
            &AlignConfig {
                additional_equalities: &pairs,
                ..Default::default()
            }
        )
        .unwrap()
        .edit_distance(),
        Some(1)
    );
}

#[test]
fn every_byte_and_equality_direction_is_representable() {
    let pairs: Vec<_> = (0..=255u8)
        .map(|byte| EqualityPair::new(byte, byte.wrapping_add(1)))
        .collect();
    for query in 0..=255u8 {
        for target in 0..=255u8 {
            let config = AlignConfig {
                additional_equalities: &pairs,
                ..Default::default()
            };
            let expected = usize::from(
                query != target
                    && query.wrapping_add(1) != target
                    && target.wrapping_add(1) != query,
            );
            assert_eq!(
                align(&[query], &[target], &config).unwrap().edit_distance(),
                Some(expected),
                "query={query} target={target}"
            );
        }
    }
    let ascending: Vec<u8> = (0..=255).collect();
    let descending: Vec<u8> = (0..=255).rev().collect();
    exercise(&ascending, &descending, &pairs);
    exercise(&ascending, &ascending, &[]);
}

#[test]
fn all_word_boundaries_and_target_boundaries_match_scalar_dp() {
    let mut rng = Rng::new(SEED);
    for query_len in BOUNDARY_LENGTHS {
        let query = rng.bytes(query_len);
        for target_len in BOUNDARY_LENGTHS {
            let target = rng.bytes(target_len);
            exercise(&query, &target, &[]);
        }
        exercise(&query, &mutated(&mut rng, &query), &[]);
        let mut embedded = vec![255, 0, 128];
        embedded.extend_from_slice(&query);
        embedded.extend([128, 0, 255]);
        exercise(&query, &embedded, &[]);
    }
}

#[test]
fn seeded_all_byte_inputs_with_changing_equalities_and_configs() {
    let mut rng = Rng::new(SEED ^ 0xa11b_17e5);
    let pairs = [
        EqualityPair::new(0, 255),
        EqualityPair::new(255, 128),
        EqualityPair::new(128, 17),
        EqualityPair::new(b'N', b'A'),
    ];
    for case in 0..192 {
        let query_len = (rng.next() % 321) as usize;
        let query = rng.bytes(query_len);
        let target = if case % 3 == 0 {
            mutated(&mut rng, &query)
        } else {
            let target_len = (rng.next() % 321) as usize;
            rng.bytes(target_len)
        };
        exercise(&query, &target, if case % 2 == 0 { &pairs } else { &[] });
    }
}

#[test]
fn exact_empty_prefix_locations_preserve_full_word_reporting_convention() {
    for len in [1, 63, 64, 65, 127, 128, 129] {
        let query = vec![b'A'; len];
        let target = vec![b'B'; 3];
        for mode in [AlignmentMode::Prefix, AlignmentMode::Infix] {
            for task in TASKS {
                let config = AlignConfig {
                    mode,
                    task,
                    ..Default::default()
                };
                let result = align(&query, &target, &config).unwrap();
                let first: isize = if len % 64 == 0 { 0 } else { -1 };
                let last = if mode == AlignmentMode::Prefix {
                    len.min(target.len())
                } else {
                    target.len()
                } as isize
                    - 1;
                assert_eq!(result.edit_distance(), Some(len));
                assert_eq!(
                    result.end_locations(),
                    &(first..=last).collect::<Vec<_>>(),
                    "len={len} {config:?}"
                );
                if task != AlignmentTask::Distance {
                    let starts: Vec<_> = (first..=last)
                        .map(|end| {
                            if mode == AlignmentMode::Infix {
                                ((end + 1) as usize).saturating_sub(len)
                            } else {
                                0
                            }
                        })
                        .collect();
                    assert_eq!(result.start_locations(), Some(starts.as_slice()));
                }
                assert_snapshot(
                    &query,
                    &target,
                    &config,
                    &quadratic(&query, &target, mode, &[]),
                    &Snapshot::capture(&result),
                );
            }
        }
    }
}

#[test]
fn empty_inputs_paths_and_every_empty_infix_boundary_are_distinct_from_no_match() {
    for query in [b"".as_slice(), b"AC"] {
        for target in [b"".as_slice(), b"AC"] {
            exercise(query, target, &[]);
        }
    }
    for mode in MODES {
        let config = AlignConfig {
            mode,
            task: AlignmentTask::Path,
            max_distance: Some(0),
            ..Default::default()
        };
        let result = align(b"", b"", &config).unwrap();
        assert_eq!(result.path(), Some([].as_slice()));
        assert_eq!(result.cigar(CigarFormat::Extended).unwrap().to_string(), "");
    }
    let result = align(
        b"",
        b"AC",
        &AlignConfig {
            mode: AlignmentMode::Infix,
            task: AlignmentTask::Path,
            max_distance: Some(0),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.end_locations(), &[-1, 0, 1]);
    assert_eq!(result.start_locations(), Some([0, 1, 2].as_slice()));
    assert_eq!(
        result.locations().unwrap().collect::<Vec<_>>(),
        vec![0..0, 1..1, 2..2]
    );
}

#[test]
fn starts_choose_earliest_optimum_and_path_covers_only_first_span() {
    let config = AlignConfig {
        mode: AlignmentMode::Infix,
        task: AlignmentTask::Path,
        ..Default::default()
    };
    let result = align(b"AC", b"TACGAC", &config).unwrap();
    assert_eq!(result.edit_distance(), Some(0));
    assert_eq!(result.start_locations(), Some([1, 4].as_slice()));
    assert_eq!(result.end_locations(), &[2, 5]);
    assert_eq!(
        result.locations().unwrap().collect::<Vec<_>>(),
        vec![1..3, 4..6]
    );
    assert_eq!(result.path(), Some([EditOperation::Match; 2].as_slice()));
    // At the first endpoint both inserting and mismatching can be optimal.
    // Earlier starts choose the longer span, not an arbitrary traceback start.
    exercise(b"AA", b"BA", &[]);
    exercise(b"ABC", b"XXABXXABCABC", &[]);
}

#[test]
fn cigar_literals_preserve_orientation_and_merge_standard_diagonals() {
    for (query, target, extended, standard) in [
        (b"ACGT".as_slice(), b"AGGT".as_slice(), "1=1X2=", "4M"),
        (b"ACGT", b"AGT", "1=1I2=", "1M1I2M"),
        (b"AGT", b"ACGT", "1=1D2=", "1M1D2M"),
        (b"", b"AC", "2D", "2D"),
        (b"AC", b"", "2I", "2I"),
        (b"", b"", "", ""),
    ] {
        let result = align(
            query,
            target,
            &AlignConfig {
                task: AlignmentTask::Path,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            result.cigar(CigarFormat::Extended).unwrap().to_string(),
            extended
        );
        assert_eq!(
            result.cigar(CigarFormat::Standard).unwrap().to_string(),
            standard
        );
    }
    let pairs = [EqualityPair::new(255, 128)];
    let result = align(
        &[0, 255],
        &[0, 128],
        &AlignConfig {
            task: AlignmentTask::Path,
            additional_equalities: &pairs,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        result.cigar(CigarFormat::Extended).unwrap().to_string(),
        "2="
    );
}

#[test]
fn adaptive_distance_and_strong_band_transitions_have_closed_form_answers() {
    // Cross adaptive doublings and the strong reduction at target columns
    // 0/2048/4096, without running quadratic DP on every bound and task.
    for distance in [1, 31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 256, 257] {
        let (query, target, expected) = known_substitutions(4_097, distance);
        for mode in MODES {
            for task in TASKS {
                for limit in limits(distance) {
                    let config = AlignConfig {
                        mode,
                        task,
                        max_distance: limit,
                        ..Default::default()
                    };
                    let result = align(&query, &target, &config).unwrap();
                    assert_snapshot(
                        &query,
                        &target,
                        &config,
                        &expected,
                        &Snapshot::capture(&result),
                    );
                }
            }
        }
    }
    // Length imbalance forces first-block movement; free-end modes still find
    // the interior exact match across a strong-reduction boundary.
    let mut rng = Rng::new(SEED ^ 2048);
    let query = rng.bytes(129);
    let mut target = vec![254; 2_049];
    target.extend_from_slice(&query);
    target.extend([253; 67]);
    exercise(&query, &target, &[]);
}

#[test]
fn twenty_thousand_byte_paths_are_optimal_replayable_and_deterministic() {
    // A 1,177-wide edit band over 20k columns exceeds the trace-memory cutoff
    // and exercises Hirschberg without a 400-million-cell scalar oracle.
    let (query, target, expected) = known_substitutions(20_000, 1_177);
    for mode in MODES {
        let config = AlignConfig {
            mode,
            task: AlignmentTask::Path,
            ..Default::default()
        };
        let first = Snapshot::capture(&align(&query, &target, &config).unwrap());
        assert_snapshot(&query, &target, &config, &expected, &first);
        let second = Snapshot::capture(&align(&query, &target, &config).unwrap());
        assert_eq!(first, second, "mode={mode:?} seed={SEED:#x}");
        let miss = align(
            &query,
            &target,
            &AlignConfig {
                max_distance: Some(1_176),
                ..config
            },
        )
        .unwrap();
        assert_eq!(miss.edit_distance(), None);
        assert!(miss.path().is_none());
    }
}

#[test]
fn owned_results_are_send_sync_and_outlive_inputs_equalities_and_config() {
    fn send_sync<T: Send + Sync>() {}
    send_sync::<Alignment>();
    send_sync::<EqualityPair>();
    send_sync::<AlignConfig<'static>>();
    let result = {
        let query = vec![0, 255, b'N'];
        let target = vec![0, 128, b'A'];
        let pairs = [EqualityPair::new(255, 128), EqualityPair::new(b'N', b'A')];
        let config = AlignConfig {
            task: AlignmentTask::Path,
            additional_equalities: &pairs,
            ..Default::default()
        };
        align(&query, &target, &config).unwrap()
    };
    let result = std::thread::spawn(move || {
        assert_eq!(result.edit_distance(), Some(0));
        result
    })
    .join()
    .unwrap();
    std::thread::scope(|scope| {
        scope.spawn(|| assert_eq!(result.path(), Some([EditOperation::Match; 3].as_slice())));
        scope.spawn(|| {
            assert_eq!(
                result.cigar(CigarFormat::Extended).unwrap().to_string(),
                "3="
            )
        });
    });
}

#[test]
fn repeated_calls_with_reused_mutating_configs_do_not_retain_state() {
    let pairs = [EqualityPair::new(0, 255), EqualityPair::new(255, 128)];
    let mut config = AlignConfig::default();
    for iteration in 0..96 {
        config.mode = MODES[iteration % MODES.len()];
        config.task = TASKS[(iteration / 3) % TASKS.len()];
        config.max_distance = [Some(0), None, Some(usize::MAX), Some(1)][iteration % 4];
        config.additional_equalities = if iteration % 2 == 0 { &pairs } else { &[] };
        for (query, target) in [
            (b"".as_slice(), b"".as_slice()),
            (&[0, 128][..], &[255, 128][..]),
            (b"ABC", b"ZABCDABC"),
        ] {
            let expected = quadratic(query, target, config.mode, config.additional_equalities);
            let first = Snapshot::capture(&align(query, target, &config).unwrap());
            assert_snapshot(query, target, &config, &expected, &first);
            assert_eq!(
                first,
                Snapshot::capture(&align(query, target, &config).unwrap())
            );
        }
    }
}

#[test]
fn twenty_thousand_byte_hirschberg_paths_handle_both_indel_orientations() {
    let (query, target, expected) = known_insertions(20_000, 23);
    for (query, target, expected) in [
        (query.as_slice(), target.as_slice(), expected.clone()),
        (
            target.as_slice(),
            query.as_slice(),
            Expected {
                ends: vec![query.len() as isize - 1],
                ..expected
            },
        ),
    ] {
        for limit in [None, Some(expected.distance), Some(expected.distance - 1)] {
            let config = AlignConfig {
                task: AlignmentTask::Path,
                max_distance: limit,
                ..Default::default()
            };
            let result = Snapshot::capture(&align(query, target, &config).unwrap());
            assert_snapshot(query, target, &config, &expected, &result);
        }
    }
}
