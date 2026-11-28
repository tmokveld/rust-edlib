//! Literal global shortcuts and exact continuation after a diagonal probe.

use super::{AlignError, EditOperation, Run, RunAlignment, check_lengths};
use crate::kernel::{self, Profile, Sequence};
use crate::path::{self, Output};
use crate::word;

/// Reusable probe for global Levenshtein alignment with literal byte equality.
///
/// Construction trims matching ends and tries to prove that the diagonal is
/// optimal. If inconclusive, callers can try another algorithm before using
/// [`Self::distance_at_most`] or [`Self::align_runs`] to finish the search.
pub struct GlobalProbe<'a> {
    core: Core<'a>,
    diagonal: Option<Diagonal>,
}

impl<'a> GlobalProbe<'a> {
    /// Borrows both inputs and eagerly tries the diagonal proof, without allocating.
    pub fn new(query: &'a [u8], target: &'a [u8]) -> Self {
        let core = Core::new(query, target);
        Self { core, diagonal: certified_diagonal::<true>(core.query, core.target) }
    }

    /// The proven optimal score, or `None` if the diagonal proof was inconclusive.
    /// Unequal lengths are inconclusive; `None` never means a cutoff miss.
    pub fn diagonal_score(&self) -> Option<usize> {
        self.diagonal.map(|diagonal| diagonal.score)
    }

    /// Exact global score up to the inclusive cutoff `k`, without repeating the proof.
    ///
    /// # Errors
    /// Returns [`AlignError::InputTooLong`] if the required dimensions or buffers
    /// cannot be accommodated.
    pub fn distance_at_most(&self, k: usize) -> Result<Option<usize>, AlignError> {
        if let Some(score) = self.diagonal_score() {
            return Ok((score <= k).then_some(score));
        }
        self.core.bounded::<false>(k)
    }

    /// One optimal global alignment as forward, coalesced runs, with an optional
    /// inclusive cutoff. Does not repeat the diagonal proof.
    ///
    /// # Errors
    /// Returns [`AlignError::InputTooLong`] if the required dimensions or buffers
    /// cannot be accommodated.
    pub fn align_runs(&self, k: Option<usize>) -> Result<Option<RunAlignment>, AlignError> {
        if let Some(diagonal) = self.diagonal {
            return self.core.finish_diagonal(diagonal, k);
        }
        self.core.align_runs::<false>(k)
    }
}

pub(super) fn distance(query: &[u8], target: &[u8], hint: Option<usize>) -> Result<usize, AlignError> {
    Core::new(query, target).distance(hint)
}

pub(super) fn distance_at_most(query: &[u8], target: &[u8], k: usize) -> Result<Option<usize>, AlignError> {
    if k == 0 {
        return Ok((query == target).then_some(0));
    }
    if query.len().abs_diff(target.len()) > k {
        return Ok(None);
    }
    Core::new(query, target).bounded::<true>(k)
}

pub(super) fn align_runs(
    query: &[u8], target: &[u8], k: Option<usize>,
) -> Result<Option<RunAlignment>, AlignError> {
    if k.is_some_and(|k| query.len().abs_diff(target.len()) > k) {
        return Ok(None);
    }
    Core::new(query, target).align_runs::<true>(k)
}

#[derive(Clone, Copy)]
struct Core<'a> {
    query: &'a [u8],
    target: &'a [u8],
    prefix: usize,
    suffix: usize,
}

impl<'a> Core<'a> {
    fn new(query: &'a [u8], target: &'a [u8]) -> Self {
        // An optimal literal global path can match common end bytes for free.
        let prefix = common_end::<false>(query, target);
        let query = &query[prefix..];
        let target = &target[prefix..];
        let suffix = common_end::<true>(query, target);
        Self {
            query: &query[..query.len() - suffix],
            target: &target[..target.len() - suffix],
            prefix,
            suffix,
        }
    }

    fn shorter_first(self) -> (&'a [u8], &'a [u8]) {
        if self.query.len() <= self.target.len() { (self.query, self.target) } else { (self.target, self.query) }
    }

    fn distance(self, hint: Option<usize>) -> Result<usize, AlignError> {
        let (short, long) = self.shorter_first();
        match short.len() {
            0 => return Ok(long.len()),
            1 => return Ok(word::short(short, long)),
            _ => {}
        }

        let gap = long.len() - short.len();
        let default_slack = kernel::initial_slack(gap, short.len()).min(short.len());
        let mut slack = hint.map_or(default_slack, |d| d.saturating_sub(gap).min(short.len()));
        if let Some(diagonal) = self.maybe_diagonal::<false>(gap + slack) {
            return Ok(diagonal.score);
        }
        check_lengths(short.len(), long.len())?;
        if short.len() <= word::WIDTH {
            return Ok(word::short(short, long));
        }
        if gap + slack <= word::BAND {
            if let Some(d) = word::sliding(short, long) {
                return Ok(d);
            }
            // The sliding pass ruled out the entire small band, not just the hint.
            slack = default_slack;
        }
        let profile = Profile::new(Sequence::new(short), &[])?;
        let target = Sequence::new(long);
        let mut blocks = Vec::new();
        loop {
            if let Some(d) = kernel::distance_global(&profile, target, gap + slack, &mut blocks) {
                return Ok(d);
            }
            // Grow only the allowance above the mandatory length gap.
            slack = slack.saturating_mul(2).max(default_slack).min(short.len());
        }
    }

    fn bounded<const AUTO_PROBE: bool>(self, k: usize) -> Result<Option<usize>, AlignError> {
        let (short, long) = self.shorter_first();
        let gap = long.len() - short.len();
        if k < gap {
            return Ok(None);
        }
        match short.len() {
            0 => return Ok(Some(gap)),
            1 => {
                let score = word::short(short, long);
                return Ok((score <= k).then_some(score));
            }
            _ => {}
        }
        if k == gap {
            // Attaining the gap allows only indels in one direction. The
            // nonempty equal-length core starts at a mismatch, so gap=0 misses.
            return Ok((gap != 0 && subsequence_blocks(short, long, |_, _| {}).is_some()).then_some(gap));
        }
        if AUTO_PROBE {
            if let Some(diagonal) = self.maybe_diagonal::<false>(k) {
                return Ok((diagonal.score <= k).then_some(diagonal.score));
            }
        }
        check_lengths(short.len(), long.len())?;
        if short.len() <= word::WIDTH {
            return Ok(word::short_at_most(short, long, k));
        }
        // No global score exceeds long.len(), but k itself need not fit in isize.
        let profile = Profile::new(Sequence::new(short), &[])?;
        Ok(kernel::distance_global(&profile, Sequence::new(long), k.min(long.len()), &mut Vec::new()))
    }

    fn align_runs<const AUTO_PROBE: bool>(self, limit: Option<usize>) -> Result<Option<RunAlignment>, AlignError> {
        let (m, n) = (self.query.len(), self.target.len());
        let gap = m.abs_diff(n);
        if limit.is_some_and(|k| k < gap) {
            return Ok(None);
        }
        if m == 0 || n == 0 {
            let op = if m == 0 { EditOperation::Delete } else { EditOperation::Insert };
            return self.framed_runs(gap, usize::from(gap != 0), |runs| {
                if gap != 0 {
                    runs.push(Run { op, len: gap });
                }
            }).map(Some);
        }
        if limit == Some(0) {
            return Ok(None); // The trimmed core starts at a mismatch.
        }
        if m.min(n) == 1 {
            let single = path::OneByte::literal(self.query, self.target);
            let score = single.score();
            if limit.is_some_and(|k| score > k) {
                return Ok(None);
            }
            return self.framed_runs(score, single.run_count(), |runs| runs.extend(single.runs())).map(Some);
        }
        if limit == Some(gap) || (gap > 32 && gap >= m.min(n).div_ceil(2)) {
            if let Some(hit) = self.subsequence_runs(gap)? {
                return Ok(Some(hit));
            }
            if limit == Some(gap) {
                return Ok(None);
            }
            // A failed embedding only rules out the gap, not a higher cutoff.
        }
        if AUTO_PROBE {
            if let Some(diagonal) = self.maybe_diagonal::<true>(limit.unwrap_or(m.max(n))) {
                return self.finish_diagonal(diagonal, limit);
            }
        }
        check_lengths(m, n)?;
        let Some((score, mut runs)) = path::find_global::<Vec<Run>>(
            self.query, self.target, &[], limit.map(|k| k.min(m.max(n))), self.prefix,
        )? else { return Ok(None) };
        runs.append_run(0, EditOperation::Match, self.suffix)?;
        Ok(Some(self.result(score, runs)))
    }

    fn maybe_diagonal<const COUNT_RUNS: bool>(self, budget: usize) -> Option<Diagonal> {
        // Avoid a whole-input scan on a short core or for a tight search.
        if self.query.len() >= 256 && budget >= 32 {
            certified_diagonal::<COUNT_RUNS>(self.query, self.target)
        } else {
            None
        }
    }

    fn finish_diagonal(self, diagonal: Diagonal, limit: Option<usize>) -> Result<Option<RunAlignment>, AlignError> {
        if limit.is_some_and(|k| diagonal.score > k) {
            return Ok(None);
        }
        // A nonempty trimmed core begins and ends with mismatches.
        let core_runs = (2 * diagonal.groups).saturating_sub(1);
        self.framed_runs(diagonal.score, core_runs, |runs| {
            let mut i = 0;
            while i < self.query.len() {
                // Trimming (or the preceding match block) already found this mismatch.
                let changed = 1 + self.query[i + 1..].iter().zip(&self.target[i + 1..]).take_while(|(x, y)| x != y).count();
                runs.push(Run { op: EditOperation::Mismatch, len: changed });
                i += changed;
                let matched = common_end::<false>(&self.query[i..], &self.target[i..]);
                if matched != 0 {
                    runs.push(Run { op: EditOperation::Match, len: matched });
                    i += matched;
                }
            }
        }).map(Some)
    }

    fn subsequence_runs(self, gap: usize) -> Result<Option<RunAlignment>, AlignError> {
        let (short, long) = self.shorter_first();
        let op = if self.query.len() < self.target.len() { EditOperation::Delete } else { EditOperation::Insert };
        // Confirm the embedding and count exact runs before allocating anything.
        let mut blocks = 0;
        let Some(trailing) = subsequence_blocks(short, long, |_, _| blocks += 1) else { return Ok(None) };
        debug_assert!(trailing > 0);
        self.framed_runs(gap, 2 * blocks + 1, |runs| {
            let end = subsequence_blocks(short, long, |skip, matched| {
                debug_assert!(skip > 0);
                runs.push(Run { op, len: skip });
                runs.push(Run { op: EditOperation::Match, len: matched });
            });
            debug_assert_eq!(end, Some(trailing));
            runs.push(Run { op, len: trailing });
        }).map(Some)
    }

    fn framed_runs(self, score: usize, core_runs: usize, emit: impl FnOnce(&mut Vec<Run>)) -> Result<RunAlignment, AlignError> {
        // Every proven run consumes at least one byte of the original longer input.
        let count = usize::from(self.prefix != 0) + core_runs + usize::from(self.suffix != 0);
        let mut runs = Vec::new();
        runs.try_reserve_exact(count).map_err(|_| AlignError::InputTooLong)?;
        if self.prefix != 0 {
            runs.push(Run { op: EditOperation::Match, len: self.prefix });
        }
        emit(&mut runs);
        if self.suffix != 0 {
            runs.push(Run { op: EditOperation::Match, len: self.suffix });
        }
        debug_assert_eq!(runs.len(), count);
        Ok(self.result(score, runs))
    }

    fn result(self, score: usize, runs: Vec<Run>) -> RunAlignment {
        RunAlignment { score, target: 0..self.prefix + self.target.len() + self.suffix, runs }
    }
}

#[derive(Clone, Copy)]
struct Diagonal {
    score: usize,
    groups: usize,
}

fn certified_diagonal<const COUNT_RUNS: bool>(query: &[u8], target: &[u8]) -> Option<Diagonal> {
    if query.len() != target.len() {
        return None;
    }
    // Matched positions cancel from the byte counts. With disjoint mismatch
    // alphabets, each changed source byte is surplus; one edit removes at most
    // one surplus byte, and the diagonal changes exactly that many bytes.
    let (mut from, mut to) = ([0u64; 4], [0u64; 4]);
    let (mut score, mut groups, mut previous_end) = (0, 0, 0);
    for (i, (&x, &y)) in query.iter().zip(target).enumerate() {
        if x == y {
            continue;
        }
        let (xi, yi) = (usize::from(x / 64), usize::from(y / 64));
        let (xb, yb) = (1u64 << (x % 64), 1u64 << (y % 64));
        if to[xi] & xb != 0 || from[yi] & yb != 0 {
            return None;
        }
        from[xi] |= xb;
        to[yi] |= yb;
        if COUNT_RUNS {
            groups += usize::from(score == 0 || i != previous_end);
            previous_end = i + 1;
        }
        score += 1;
    }
    Some(Diagonal { score, groups })
}

// Earliest matches leave the most room for all remaining bytes. In a nonempty
// trimmed core, a successful pass alternates skips and match blocks, then has
// a trailing skip because the last bytes disagree.
fn subsequence_blocks(mut short: &[u8], mut long: &[u8], mut block: impl FnMut(usize, usize)) -> Option<usize> {
    while let Some((&first, rest)) = short.split_first() {
        let spare = long.len() - short.len();
        let skip = long[..=spare].iter().position(|&b| b == first)?;
        let tail = &long[skip + 1..];
        let extra = common_end::<false>(rest, tail);
        block(skip, 1 + extra);
        short = &rest[extra..];
        long = &tail[extra..];
    }
    Some(long.len())
}

#[inline]
fn chunk<const REVERSE: bool>(bytes: &[u8]) -> Option<(&[u8; 32], &[u8])> {
    if REVERSE {
        bytes.split_last_chunk().map(|(rest, last)| (last, rest))
    } else {
        bytes.split_first_chunk()
    }
}

fn common_end<const REVERSE: bool>(mut a: &[u8], mut b: &[u8]) -> usize {
    // Skip equal blocks, then inspect only the first unequal or partial block.
    let mut matched = 0;
    while let (Some((x, tail_a)), Some((y, tail_b))) = (chunk::<REVERSE>(a), chunk::<REVERSE>(b)) {
        if x != y {
            break;
        }
        matched += x.len();
        (a, b) = (tail_a, tail_b);
    }
    matched + if REVERSE {
        a.iter().rev().zip(b.iter().rev()).take_while(|(x, y)| x == y).count()
    } else {
        a.iter().zip(b).take_while(|(x, y)| x == y).count()
    }
}
