//! Owned alignment results and borrowed configuration.

use crate::kernel::{self, Profile, Sequence};
use crate::myers;
use std::{fmt, ops::Range};

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
}

/// More detailed tasks require more computation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AlignmentTask {
    /// Compute distance and every reference-reported optimal end location.
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

fn unbounded_distance(upper: usize, mut search: impl FnMut(usize) -> Option<usize>) -> usize {
    let mut k = 64.min(upper);
    loop {
        if let Some(distance) = search(k) {
            return distance;
        }
        // `upper` is a feasible edit distance, so the search succeeds by then.
        k = k.saturating_mul(2).min(upper);
    }
}

/// Exact global Levenshtein distance between two byte slices, without a bound.
/// Computes only the score, not locations or an edit path.
///
/// ```
/// assert_eq!(rust_edlib::distance(b"kitten", b"sitting")?, 3);
/// # Ok::<(), rust_edlib::AlignError>(())
/// ```
pub fn distance(query: &[u8], target: &[u8]) -> Result<usize, AlignError> {
    let (short, long) = if query.len() <= target.len() {
        (query, target)
    } else {
        (target, query)
    };
    check_lengths(short.len(), long.len())?;
    if short.is_empty() {
        return Ok(long.len());
    }
    let profile = Profile::new(Sequence::new(short), &[])?;
    let target = Sequence::new(long);
    Ok(unbounded_distance(long.len(), |k| {
        kernel::distance_global(&profile, target, k)
    }))
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
    anchored_distance::<false>(query, target)
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
    anchored_distance::<true>(query, target)
}

fn anchored_distance<const REVERSE: bool>(
    query: &[u8],
    target: &[u8],
) -> Result<usize, AlignError> {
    let shared = if REVERSE {
        query.iter().rev().zip(target.iter().rev()).take_while(|(a, b)| a == b).count()
    } else {
        query.iter().zip(target).take_while(|(a, b)| a == b).count()
    };
    let m = query.len() - shared;
    let n = target.len() - shared;
    if m == 0 || n == 0 {
        return Ok(m);
    }
    // The empty prefix/suffix costs m; consuming over 2m target bytes costs
    // more than m from the length gap alone.
    let extent = n.min(m.saturating_mul(2));
    let (query, target) = if REVERSE {
        (&query[..m], &target[n - extent..n])
    } else {
        (&query[shared..], &target[shared..shared + extent])
    };
    // Both slices are nonempty and their first bytes in traversal order differ.
    myers::distance::<REVERSE>(query, target, m.saturating_sub(extent).max(1))
}

/// Align arbitrary bytes without copying the inputs. The result borrows neither
/// inputs nor configuration. A bound exceeded is `Ok` with distance `None`.
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
    check_lengths(m, n)?;
    kernel::check_buffer::<isize>(n.checked_add(1).ok_or(AlignError::InputTooLong)?)?;
    let upper = if config.mode == AlignmentMode::Global {
        m.max(n)
    } else {
        m
    };
    let limit = config.max_distance.map(|k| k.min(upper));
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
            result.ends.push(if config.mode == AlignmentMode::Global {
                n as isize - 1
            } else {
                -1
            });
        }
        if config.task != AlignmentTask::Distance {
            result.starts = Some(if empty_infix {
                (0..=n).collect()
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

    let profile = Profile::new(Sequence::new(query), config.additional_equalities)?;
    let mut search = |k| {
        if config.mode == AlignmentMode::Global {
            kernel::distance_global(&profile, Sequence::new(target), k)
        } else {
            kernel::semiglobal(
                &profile,
                Sequence::new(target),
                k,
                config.mode == AlignmentMode::Infix,
                &mut result.ends,
            )
        }
    };
    result.distance = match limit {
        Some(k) => search(k),
        None => Some(unbounded_distance(upper, search)),
    };
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
                    Sequence::reversed(&target[..offset]),
                    distance,
                    false,
                    &mut reverse_ends,
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
        let start = result.starts.as_ref().unwrap()[0];
        let end = (result.ends[0] + 1) as usize;
        // The forward profile depends only on the query and equalities. The
        // reconstruction's top level reuses it and releases it before recursing.
        result.path = Some(crate::path::reconstruct(
            query,
            &target[start..end],
            config.additional_equalities,
            distance,
            profile,
        )?);
    }
    Ok(result)
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
