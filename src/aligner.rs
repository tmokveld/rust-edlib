//! Owned alignment results and borrowed configuration.

use crate::kernel::{self, Profile, Sequence};
use crate::myers;
use std::{fmt, ops::Range};

mod global;
mod infix;
pub use global::GlobalProbe;

/// Which target gaps are free.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AlignmentMode {
    /// Align the entire query to the entire target (NW).
    #[default]
    Global,
    /// Align the query to a target prefix; trailing target gaps are free (SHW).
    Prefix,
    /// Align the query anywhere in the target; both target flanks are free (HW).
    Infix,
    /// Align the query to a target suffix; leading target gaps are free.
    /// The only end is the target length, and ties choose the earliest start.
    Suffix,
}

/// More detailed tasks require more computation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AlignmentTask {
    /// Compute distance and every reported optimal end location.
    #[default]
    Distance,
    /// Also compute the corresponding start locations.
    Locations,
    /// Also compute one optimal path for the first reported span.
    Path,
}

/// A symmetric additional byte equality. Equalities are not transitively closed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EqualityPair {
    pub(crate) first: u8,
    pub(crate) second: u8,
}

impl EqualityPair {
    /// Declare two bytes equal in either order.
    pub const fn new(first: u8, second: u8) -> Self {
        Self { first, second }
    }
    /// First byte in the pair.
    pub const fn first(self) -> u8 {
        self.first
    }
    /// Second byte in the pair.
    pub const fn second(self) -> u8 {
        self.second
    }
}

/// Borrowed configuration reusable across independent alignments.
#[derive(Clone, Copy, Debug, Default)]
pub struct AlignConfig<'a> {
    /// Alignment mode; defaults to global.
    pub mode: AlignmentMode,
    /// Requested result detail; defaults to distance and end locations.
    pub task: AlignmentTask,
    /// Inclusive distance bound; `None` uses adaptive band expansion.
    pub max_distance: Option<usize>,
    /// Symmetric, reflexive, non-transitive extensions of byte equality.
    pub additional_equalities: &'a [EqualityPair],
}

/// Invalid request, distinct from a successful bounded search without a match.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignError {
    /// Required address, buffer-size, or signed band arithmetic cannot fit.
    InputTooLong,
}

impl fmt::Display for AlignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("alignment dimensions exceed addressable size or band arithmetic")
    }
}

impl std::error::Error for AlignError {}

fn check_lengths(m: usize, n: usize) -> Result<(), AlignError> {
    // Band expressions combine the two dimensions, the bound, and word padding.
    m.checked_add(n)
        .and_then(|x| x.checked_add(m.max(n)))
        .and_then(|x| x.checked_add(128))
        .filter(|&x| x <= isize::MAX as usize)
        .ok_or(AlignError::InputTooLong)?;
    Ok(())
}

pub(crate) fn search_distance(
    upper: usize, limit: Option<usize>, mut search: impl FnMut(usize) -> Option<usize>,
) -> Option<usize> {
    if let Some(k) = limit {
        return search(k);
    }
    let mut k = 64.min(upper);
    loop {
        if let Some(distance) = search(k) {
            return Some(distance);
        }
        // `upper` is a feasible edit distance, so the search succeeds by then.
        k = k.saturating_mul(2).min(upper);
    }
}

/// Exact global Levenshtein distance between two byte slices, without a bound.
/// Computes only the score, not locations or an edit path.
///
/// ```
/// assert_eq!(rust_edlib::distance(b"ACGT", b"AGGGA")?, 3);
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn distance(query: &[u8], target: &[u8]) -> Result<usize, AlignError> {
    global::distance(query, target, None)
}

/// Exact global [`distance`] using an initial estimate of the total edit
/// distance, including the mandatory length difference. The hint is not a
/// cutoff: if a narrow search misses, the search widens until it finds the
/// exact score. Hints outside the feasible score range are clamped.
/// Computes only the score, not locations or an edit path.
///
/// # Errors
/// Returns [`AlignError::InputTooLong`] if the required dimensions or buffers
/// cannot be accommodated.
///
/// ```
/// assert_eq!(rust_edlib::distance_with_hint(b"ACGT", b"AGGGA", 1)?, 3);
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn distance_with_hint(query: &[u8], target: &[u8], hint: usize) -> Result<usize, AlignError> {
    global::distance(query, target, Some(hint))
}

/// Exact global [`distance`] up to the inclusive cutoff `k`, without collecting
/// locations or an edit path. Returns `Ok(None)` if the distance exceeds `k`.
///
/// # Errors
/// Returns [`AlignError::InputTooLong`] if the required dimensions or buffers
/// cannot be accommodated.
///
/// ```
/// use rust_edlib::distance_at_most;
///
/// assert_eq!(distance_at_most(b"ACGT", b"AGGGA", 3)?, Some(3));
/// assert_eq!(distance_at_most(b"ACGT", b"AGGGA", 2)?, None);
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn distance_at_most(
    query: &[u8],
    target: &[u8],
    k: usize,
) -> Result<Option<usize>, AlignError> {
    global::distance_at_most(query, target, k)
}

/// Exact Levenshtein distance from the entire `query` to the best prefix of
/// `target`, including its empty prefix. The trailing target bytes are free.
/// Computes only the score, not locations or an edit path.
///
/// # Errors
/// Returns [`AlignError::InputTooLong`] if the search dimensions or required
/// buffers cannot be accommodated.
///
/// ```
/// assert_eq!(rust_edlib::prefix_distance(b"AC", b"ACTG")?, 0);
/// assert_eq!(rust_edlib::prefix_distance(b"ACTG", b"AC")?, 2);
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn prefix_distance(query: &[u8], target: &[u8]) -> Result<usize, AlignError> {
    Ok(prefix_distance_at_most(query, target, query.len())?
        .expect("the empty prefix is within the query-length bound"))
}

/// Exact [`prefix_distance`] up to the inclusive cutoff `k`, without collecting
/// endpoints or a path. Returns `Ok(None)` if even the best prefix costs more.
///
/// # Errors
/// Returns [`AlignError::InputTooLong`] if the required dimensions or buffers
/// cannot be accommodated.
///
/// ```
/// assert_eq!(rust_edlib::prefix_distance_at_most(b"ACG", b"ATGCA", 1)?, Some(1));
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn prefix_distance_at_most(
    query: &[u8], target: &[u8], k: usize,
) -> Result<Option<usize>, AlignError> {
    anchored_distance::<false>(query, target, k)
}

/// Exact Levenshtein distance from the entire `query` to the best suffix of
/// `target`, including its empty suffix. The leading target bytes are free.
/// Computes only the score, not locations or an edit path.
///
/// # Errors
/// Returns [`AlignError::InputTooLong`] if the search dimensions or required
/// buffers cannot be accommodated.
///
/// ```
/// assert_eq!(rust_edlib::suffix_distance(b"TG", b"ACTG")?, 0);
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn suffix_distance(query: &[u8], target: &[u8]) -> Result<usize, AlignError> {
    Ok(suffix_distance_at_most(query, target, query.len())?
        .expect("the empty suffix is within the query-length bound"))
}

/// Exact [`suffix_distance`] up to the inclusive cutoff `k`, without collecting
/// endpoints or a path. Returns `Ok(None)` if even the best suffix costs more.
///
/// # Errors
/// Returns [`AlignError::InputTooLong`] if the required dimensions or buffers
/// cannot be accommodated.
///
/// ```
/// assert_eq!(rust_edlib::suffix_distance_at_most(b"ACG", b"TTATG", 0)?, None);
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn suffix_distance_at_most(
    query: &[u8], target: &[u8], k: usize,
) -> Result<Option<usize>, AlignError> {
    anchored_distance::<true>(query, target, k)
}

fn anchored_distance<const REVERSE: bool>(
    query: &[u8],
    target: &[u8],
    k: usize,
) -> Result<Option<usize>, AlignError> {
    if query.len().saturating_sub(target.len()) > k {
        return Ok(None);
    }
    let shared = if REVERSE {
        query.iter().rev().zip(target.iter().rev()).take_while(|(a, b)| a == b).count()
    } else {
        query.iter().zip(target).take_while(|(a, b)| a == b).count()
    };
    let m = query.len() - shared;
    let n = target.len() - shared;
    if m == 0 || n == 0 {
        return Ok(Some(m));
    }
    if k == 0 {
        return Ok(None); // The first remaining bytes at the pinned end differ.
    }
    // Deleting the remaining query costs m. Any longer target flank than
    // m + min(k, m) loses on length difference alone.
    let limit = k.min(m);
    let extent = n.min(m.saturating_add(limit));
    let (query, target) = if REVERSE {
        (&query[..m], &target[n - extent..n])
    } else {
        (&query[shared..], &target[shared..shared + extent])
    };
    // Both slices are nonempty and their first bytes in traversal order differ.
    myers::distance::<REVERSE>(query, target, m.saturating_sub(extent).max(1), limit)
}

/// Minimum Levenshtein distance from the entire `query` to any contiguous
/// substring of `target`, including the empty substring. Both target flanks
/// are free; bytes match only when equal.
///
/// Unlike [`align`] with [`AlignmentTask::Distance`], this search does not
/// collect endpoints. It stops at the first exact match.
///
/// # Errors
/// Returns [`AlignError::InputTooLong`] if the required dimensions or buffers
/// cannot be accommodated.
///
/// ```
/// assert_eq!(rust_edlib::infix_distance(b"AC", b"TACGAC")?, 0);
/// assert_eq!(rust_edlib::infix_distance(b"AC", b"")?, 2);
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn infix_distance(query: &[u8], target: &[u8]) -> Result<usize, AlignError> {
    Ok(infix_distance_at_most(query, target, query.len())?
        .expect("the empty substring is within the query-length bound"))
}

/// Exact [`infix_distance`] up to the inclusive cutoff `k`. Returns `Ok(Some(d))`
/// for the minimum score `d` if it is within the cutoff, or `Ok(None)` otherwise.
/// A zero cutoff searches for an exact byte substring without building a Myers profile.
///
/// # Errors
/// Returns [`AlignError::InputTooLong`] if the required dimensions or buffers
/// cannot be accommodated.
///
/// ```
/// use rust_edlib::infix_distance_at_most;
///
/// assert_eq!(infix_distance_at_most(b"ab", b"ax---ab", 1)?, Some(0));
/// assert_eq!(infix_distance_at_most(b"ab", b"ax", 0)?, None);
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn infix_distance_at_most(
    query: &[u8],
    target: &[u8],
    k: usize,
) -> Result<Option<usize>, AlignError> {
    let (m, n) = (query.len(), target.len());
    check_lengths(m, n)?;
    let gap = m.saturating_sub(n);
    if gap > k {
        return Ok(None);
    }
    if m == 0 || n == 0 {
        return Ok(Some(m));
    }
    if k == 0 {
        return Ok(infix::find_exact(query, target)?.map(|_| 0));
    }
    let profile = Profile::new(Sequence::new(query), &[])?;
    let empty = (k >= m).then_some(m);
    let limit = k.min(m - 1); // Only search for scores better than the empty substring.
    let mut bound = if empty.is_some() {
        (gap + 64).min(limit)
    } else {
        limit
    };
    let mut blocks = Vec::new();
    loop {
        let score = kernel::semiglobal(
            &profile, target.iter().copied(), bound, true, None, &mut blocks,
        );
        if score.is_some() || bound == limit {
            return Ok(score.or(empty));
        }
        // A narrow miss rules out this budget, not a match. Replay from column zero.
        bound = bound.saturating_mul(2).min(limit);
    }
}

/// Align arbitrary bytes without copying the inputs. The result borrows neither
/// inputs nor configuration. A bound exceeded is `Ok` with distance `None`.
///
/// Even [`AlignmentTask::Distance`] collects optimal ends. For an infix score
/// alone with ordinary byte equality, use [`infix_distance`] or [`infix_distance_at_most`].
///
/// [`AlignmentMode::Suffix`] reports only the end of `target`. Locations and
/// paths use the longest optimal suffix, in the original target coordinates.
/// An empty query uses the empty span at the end of `target`.
///
/// End locations are inclusive for binding compatibility; prefer `locations()`
/// and `end_offsets()`. In particular, inclusive end `-1` maps to offset zero.
///
/// ```
/// use rust_edlib::{align, AlignConfig, AlignmentTask, CigarFormat};
/// let config = AlignConfig { task: AlignmentTask::Path, ..Default::default() };
/// let result = align(b"ACGT", b"AGGT", &config)?;
/// assert_eq!(result.edit_distance(), Some(1));
/// assert_eq!(result.locations().unwrap().next(), Some(0..4));
/// assert_eq!(result.cigar(CigarFormat::Extended).unwrap().to_string(), "1=1X2=");
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn align(
    query: &[u8],
    target: &[u8],
    config: &AlignConfig<'_>,
) -> Result<Alignment, AlignError> {
    let (m, n) = (query.len(), target.len());
    let suffix = config.mode == AlignmentMode::Suffix;
    let upper = if config.mode == AlignmentMode::Global {
        m.max(n)
    } else {
        m
    };
    let limit = config.max_distance.map(|k| k.min(upper));
    if suffix {
        // Longer suffixes cannot meet the cutoff on length difference alone.
        let extent = n.min(m.saturating_add(limit.unwrap_or(upper)));
        check_lengths(m, extent)?;
    } else {
        check_lengths(m, n)?;
        kernel::check_buffer::<isize>(n.checked_add(1).ok_or(AlignError::InputTooLong)?)?;
    }
    let mut result = Alignment {
        distance: None,
        ends: Vec::new(),
        starts: None,
        path: None,
    };
    if m == 0 || n == 0 {
        if limit.is_some_and(|k| upper > k) {
            return Ok(result);
        }
        result.distance = Some(upper);
        let empty_infix = m == 0 && config.mode == AlignmentMode::Infix;
        if empty_infix {
            result.ends = (0..=n).map(|i| i as isize - 1).collect();
        } else {
            result.ends.push(if suffix || config.mode == AlignmentMode::Global {
                n as isize - 1
            } else {
                -1
            });
        }
        if config.task != AlignmentTask::Distance {
            result.starts = Some(if empty_infix {
                (0..=n).collect()
            } else if suffix {
                vec![n]
            } else {
                vec![0]
            });
        }
        if config.task == AlignmentTask::Path {
            result.path = Some(vec![
                if m == 0 {
                    EditOperation::Delete
                } else {
                    EditOperation::Insert
                };
                upper
            ]);
        }
        return Ok(result);
    }
    let lower = if config.mode == AlignmentMode::Global {
        m.abs_diff(n)
    } else {
        m.saturating_sub(n)
    };
    if limit.is_some_and(|k| k < lower) {
        return Ok(result);
    }

    if suffix {
        let mut consumed = 0;
        result.distance = {
            let reverse = Profile::new(Sequence::reversed(query), config.additional_equalities)?;
            let mut blocks = Vec::new();
            let search = |k| {
                if k < lower {
                    return None;
                }
                let window = target.iter().rev().take(m + k).copied();
                if config.task == AlignmentTask::Distance {
                    kernel::semiglobal(&reverse, window, k, false, None, &mut blocks)
                } else {
                    let (distance, len) = kernel::anchored_last(&reverse, window, k, &mut blocks)?;
                    consumed = len;
                    Some(distance)
                }
            };
            search_distance(upper, limit, search)
        };
        if let Some(distance) = result.distance {
            result.ends.push(n as isize - 1);
            if config.task != AlignmentTask::Distance {
                let start = n - consumed;
                result.starts = Some(vec![start]);
                if config.task == AlignmentTask::Path {
                    result.path = Some(crate::path::reconstruct(
                        query,
                        &target[start..],
                        config.additional_equalities,
                        distance,
                        None,
                    )?);
                }
            }
        }
        return Ok(result);
    }

    if config.mode == AlignmentMode::Global && config.task == AlignmentTask::Path {
        if let Some((distance, path)) =
            crate::path::find_global(query, target, config.additional_equalities, limit, 0)?
        {
            result.distance = Some(distance);
            result.ends.push(n as isize - 1);
            result.starts = Some(vec![0]);
            result.path = Some(path);
        }
        return Ok(result);
    }
    let profile = Profile::new(Sequence::new(query), config.additional_equalities)?;
    let mut blocks = Vec::new();
    let search = |k| {
        if config.mode == AlignmentMode::Global {
            kernel::distance_global(&profile, Sequence::new(target), k, &mut blocks)
        } else {
            kernel::semiglobal(
                &profile,
                target.iter().copied(),
                k,
                config.mode == AlignmentMode::Infix,
                Some(&mut result.ends),
                &mut blocks,
            )
        }
    };
    result.distance = search_distance(upper, limit, search);
    let Some(distance) = result.distance else {
        return Ok(result);
    };
    if config.mode == AlignmentMode::Global {
        result.ends.push(n as isize - 1);
    }
    if config.task != AlignmentTask::Distance {
        let mut starts = Vec::with_capacity(result.ends.len());
        if config.mode == AlignmentMode::Infix {
            let reverse = Profile::new(Sequence::reversed(query), config.additional_equalities)?;
            let mut reverse_ends = Vec::new();
            for &end in &result.ends {
                let offset = (end + 1) as usize;
                if offset == 0 {
                    starts.push(0);
                    continue;
                }
                let score = kernel::semiglobal(
                    &reverse,
                    target[..offset].iter().rev().copied(),
                    distance,
                    false,
                    Some(&mut reverse_ends),
                    &mut blocks,
                );
                debug_assert_eq!(score, Some(distance));
                // Last reverse endpoint selects the earliest optimal start.
                let consumed = (reverse_ends.last().expect("optimal reverse prefix") + 1) as usize;
                starts.push(offset - consumed);
            }
        } else {
            starts.resize(result.ends.len(), 0);
        }
        result.starts = Some(starts);
    }
    if config.task == AlignmentTask::Path {
        drop(blocks); // The search band is not needed for path reconstruction.
        let start = result.starts.as_ref().unwrap()[0];
        let end = (result.ends[0] + 1) as usize;
        // The forward profile depends only on the query and equalities. The
        // reconstruction's top level reuses it and releases it before recursing.
        result.path = Some(crate::path::reconstruct(
            query,
            &target[start..end],
            config.additional_equalities,
            distance,
            Some(profile),
        )?);
    }
    Ok(result)
}

/// Aligns the entire query to one optimal target span and returns forward,
/// coalesced runs. Does not collect all tied spans or build a unit-operation path.
///
/// The cutoff is inclusive; `Ok(None)` means no alignment meets it. For
/// [`AlignmentMode::Suffix`], ties choose the longest optimal suffix (an empty
/// query uses the empty span at the end). Other tied paths and spans are not
/// canonical.
///
/// # Errors
/// Returns [`AlignError::InputTooLong`] if the required dimensions or buffers
/// cannot be accommodated.
pub fn align_runs(
    query: &[u8],
    target: &[u8],
    mode: AlignmentMode,
    max_distance: Option<usize>,
    additional_equalities: &[EqualityPair],
) -> Result<Option<RunAlignment>, AlignError> {
    if mode == AlignmentMode::Global && additional_equalities.is_empty() {
        return global::align_runs(query, target, max_distance);
    }

    let (m, n) = (query.len(), target.len());
    let upper = if mode == AlignmentMode::Global { m.max(n) } else { m };
    let limit = max_distance.map(|k| k.min(upper));
    let extent = match mode {
        // A longer anchored span would cost more than the cutoff on length alone.
        AlignmentMode::Prefix | AlignmentMode::Suffix => n.min(m.saturating_add(limit.unwrap_or(upper))),
        _ => n,
    };
    check_lengths(m, extent)?;

    if m == 0 || n == 0 {
        if limit.is_some_and(|k| k < upper) {
            return Ok(None);
        }
        let span = match mode {
            AlignmentMode::Global => 0..n,
            AlignmentMode::Suffix => n..n,
            _ => 0..0,
        };
        let runs = crate::path::reconstruct(query, &target[span.clone()], additional_equalities, upper, None)?;
        return Ok(Some(RunAlignment { score: upper, target: span, runs }));
    }

    let lower = if mode == AlignmentMode::Global { m.abs_diff(n) } else { m.saturating_sub(n) };
    if limit.is_some_and(|k| k < lower) {
        return Ok(None);
    }
    if limit == Some(0) && mode == AlignmentMode::Infix && additional_equalities.is_empty() {
        let Some(start) = infix::find_exact(query, target)? else { return Ok(None) };
        let mut runs = Vec::new();
        runs.try_reserve_exact(1).map_err(|_| AlignError::InputTooLong)?;
        runs.push(Run { op: EditOperation::Match, len: m });
        return Ok(Some(RunAlignment { score: 0, target: start..start + m, runs }));
    }
    if mode == AlignmentMode::Global {
        let hit = crate::path::find_global(query, target, additional_equalities, limit, 0)?;
        return Ok(hit.map(|(score, runs)| RunAlignment { score, target: 0..n, runs }));
    }

    let suffix = mode == AlignmentMode::Suffix;
    let infix = mode == AlignmentMode::Infix;
    let sequence = if suffix { Sequence::reversed(query) } else { Sequence::new(query) };
    let profile = Profile::new(sequence, additional_equalities)?;
    let mut blocks = Vec::new();
    let mut consumed = 0;
    let score = search_distance(upper, limit, |k| {
        if k < lower {
            return None;
        }
        let (score, end) = if suffix {
            let window = target.iter().rev().take(m + k).copied();
            kernel::anchored_last(&profile, window, k, &mut blocks)
        } else {
            let end = if infix { n } else { n.min(m + k) };
            kernel::semiglobal_first(&profile, target[..end].iter().copied(), k, infix, &mut blocks)
        }?;
        consumed = end;
        Some(score)
    });
    let Some(score) = score else { return Ok(None) };
    let (span, forward) = if suffix {
        drop(profile);
        (n - consumed..n, None)
    } else {
        let start = if !infix || consumed == 0 {
            0
        } else if score == 0 {
            consumed - m // An exact hit has only one possible start.
        } else {
            let reverse = Profile::new(Sequence::reversed(query), additional_equalities)?;
            let window = target[..consumed].iter().rev().take(m + score).copied();
            let (reverse_score, len) = kernel::anchored_last(&reverse, window, score, &mut blocks)
                .expect("the selected optimal end has an optimal start");
            debug_assert_eq!(reverse_score, score);
            consumed - len
        };
        (start..consumed, Some(profile))
    };
    drop(blocks);
    let runs = crate::path::reconstruct(query, &target[span.clone()], additional_equalities, score, forward)?;
    Ok(Some(RunAlignment { score, target: span, runs }))
}

/// A positive number of consecutive operations of the same kind.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Run {
    /// Operation repeated by this run.
    pub op: EditOperation,
    /// Number of unit operations.
    pub len: usize,
}

/// One optimal alignment of the entire query, with no free flanks in its path.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RunAlignment {
    /// Exact unit-cost Levenshtein distance: the sum of the non-match run lengths.
    pub score: usize,
    /// Half-open byte range in the original target.
    pub target: Range<usize>,
    /// Positive-length, coalesced operations in forward order, including matches.
    pub runs: Vec<Run>,
}

/// Owns all result buffers. Accessors borrow without copying; `Send + Sync`.
#[derive(Clone, Debug)]
pub struct Alignment {
    distance: Option<usize>,
    ends: Vec<isize>,
    starts: Option<Vec<usize>>,
    path: Option<Vec<EditOperation>>,
}
impl Alignment {
    /// `None` means no alignment satisfied the configured bound.
    pub fn edit_distance(&self) -> Option<usize> {
        self.distance
    }
    /// Inclusive target ends, ordered increasingly; `-1` is boundary zero.
    pub fn end_locations(&self) -> &[isize] {
        &self.ends
    }
    /// Starts corresponding to the ends, or `None` when not computed.
    pub fn start_locations(&self) -> Option<&[usize]> {
        self.starts.as_deref()
    }
    /// Exclusive target end offsets, including zero for inclusive end `-1`.
    pub fn end_offsets(&self) -> impl ExactSizeIterator<Item = usize> + '_ {
        self.ends.iter().map(|&end| (end + 1) as usize)
    }
    /// Half-open spans. The path covers only the first span.
    pub fn locations(&self) -> Option<impl ExactSizeIterator<Item = Range<usize>> + '_> {
        Some(
            self.starts
                .as_ref()?
                .iter()
                .zip(self.end_offsets())
                .map(|(&start, end)| start..end),
        )
    }
    /// First-span path. A computed empty path is `Some(&[])`.
    pub fn path(&self) -> Option<&[EditOperation]> {
        self.path.as_deref()
    }
    /// Borrowed, allocation-free CIGAR formatter, absent if no path was computed.
    pub fn cigar(&self, format: CigarFormat) -> Option<Cigar<'_>> {
        self.path().map(|path| Cigar { path, format })
    }
}

/// Unit-cost edit operations. With query=read and target=reference, I/D follow SAM.
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EditOperation {
    /// Consume equal query and target bytes, cost zero.
    Match = 0,
    /// Consume a query byte only (insertion into target), cost one.
    Insert = 1,
    /// Consume a target byte only (deletion from target), cost one.
    Delete = 2,
    /// Consume unequal query and target bytes, cost one.
    Mismatch = 3,
}
/// How to format matches and mismatches.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarFormat {
    /// Merge adjacent matches and mismatches into `M` runs.
    Standard,
    /// Use `=` for matches and `X` for mismatches.
    Extended,
}
/// Formats directly into the caller's formatter, without an intermediate buffer.
#[derive(Clone, Copy, Debug)]
pub struct Cigar<'a> {
    path: &'a [EditOperation],
    format: CigarFormat,
}
impl fmt::Display for Cigar<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut codes = self.path.iter().map(|op| match (op, self.format) {
            (EditOperation::Insert, _) => 'I',
            (EditOperation::Delete, _) => 'D',
            (EditOperation::Match, CigarFormat::Extended) => '=',
            (EditOperation::Mismatch, CigarFormat::Extended) => 'X',
            _ => 'M',
        });
        let Some(mut previous) = codes.next() else {
            return Ok(());
        };
        let mut count = 1usize;
        for code in codes {
            if code == previous {
                count += 1;
            } else {
                write!(f, "{count}{previous}")?;
                previous = code;
                count = 1;
            }
        }
        write!(f, "{count}{previous}")
    }
}
