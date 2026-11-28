use std::mem::size_of;

use crate::aligner::{AlignError, EditOperation, EqualityPair, Run, search_distance};
use crate::kernel::{Band, Profile, Sequence, Trace, distance_global, global, trace_single};

mod cuts;

const WORD_BITS: usize = u64::BITS as usize;
const TRACE_BYTES: usize = 1024 * 1024;

/// A monomorphized output of the same forward traversal and backward leaves.
pub(crate) trait Output: Sized {
    // Unit paths preserve the historical placement across Hirschberg leaves.
    const LEGACY_SPLIT_TIES: bool = false;

    type Tail;

    fn single(operation: EditOperation, count: usize) -> Self;
    fn buffer(m: usize, n: usize, distance: usize) -> Result<Self, AlignError>;
    fn begin_leaf(&mut self) -> (usize, Self::Tail);
    fn finish_leaf(&mut self, first: usize, tail: Self::Tail) -> Result<(), AlignError>;

    // Forward appends pass 0. A backward leaf may only coalesce operations
    // written since `first`, until its first forward run is known.
    fn append_run(&mut self, first: usize, op: EditOperation, len: usize) -> Result<(), AlignError>;

    #[inline]
    fn append_step(&mut self, first: usize, op: EditOperation) -> Result<(), AlignError> {
        self.append_run(first, op, 1)
    }
}

impl Output for Vec<EditOperation> {
    const LEGACY_SPLIT_TIES: bool = true;

    type Tail = ();

    fn single(operation: EditOperation, count: usize) -> Self {
        vec![operation; count]
    }

    fn buffer(m: usize, n: usize, distance: usize) -> Result<Self, AlignError> {
        operation_buffer(m, n, distance)
    }

    #[inline]
    fn begin_leaf(&mut self) -> (usize, ()) {
        (self.len(), ())
    }

    #[inline]
    fn finish_leaf(&mut self, first: usize, (): ()) -> Result<(), AlignError> {
        self[first..].reverse();
        Ok(())
    }

    fn append_run(&mut self, _: usize, op: EditOperation, len: usize) -> Result<(), AlignError> {
        let end = self.len().checked_add(len).ok_or(AlignError::InputTooLong)?;
        debug_assert!(end <= self.capacity());
        self.resize(end, op);
        Ok(())
    }

    #[inline]
    fn append_step(&mut self, _: usize, op: EditOperation) -> Result<(), AlignError> {
        self.push(op);
        Ok(())
    }
}

impl Output for Vec<Run> {
    type Tail = Option<Run>;

    fn single(op: EditOperation, len: usize) -> Self {
        if len == 0 { Vec::new() } else { vec![Run { op, len }] }
    }

    fn buffer(_: usize, _: usize, _: usize) -> Result<Self, AlignError> {
        // Allocate for the actual runs, not for the path length or score.
        Ok(Vec::new())
    }

    fn begin_leaf(&mut self) -> (usize, Option<Run>) {
        // Hold the preceding run aside so a coalescing boundary never
        // requires allocating a temporary extra slot just to merge into it.
        let tail = self.pop();
        (self.len(), tail)
    }

    fn finish_leaf(&mut self, first: usize, tail: Option<Run>) -> Result<(), AlignError> {
        if let Some(Run { op, len }) = tail {
            self.append_run(first, op, len)?;
        }
        self[first..].reverse();
        Ok(())
    }

    #[inline]
    fn append_run(&mut self, first: usize, op: EditOperation, len: usize) -> Result<(), AlignError> {
        if len == 0 {
            return Ok(());
        }
        if self.len() > first {
            let last = self.last_mut().unwrap();
            if last.op == op {
                // Every step consumes at least one byte of the checked dimensions.
                last.len += len;
                return Ok(());
            }
        }
        if self.len() == self.capacity() {
            self.try_reserve(1).map_err(|_| AlignError::InputTooLong)?;
        }
        self.push(Run { op, len });
        Ok(())
    }
}

/// A global traceback with one byte on the shorter side. It pairs that byte
/// once and uses only gaps of the same type on either side of the pairing.
pub(crate) struct OneByte {
    gap: EditOperation,
    diagonal: EditOperation,
    before: usize,
    after: usize,
}

impl OneByte {
    pub(crate) fn literal(query: &[u8], target: &[u8]) -> Self {
        Self::find(query, target, &[], false)
    }

    fn find(query: &[u8], target: &[u8], pairs: &[EqualityPair], split: bool) -> Self {
        let (byte, long, gap) = if query.len() == 1 {
            (query[0], target, EditOperation::Delete)
        } else {
            debug_assert_eq!(target.len(), 1);
            (target[0], query, EditOperation::Insert)
        };
        let (before, matched) = if pairs.is_empty() {
            one_byte_position(long, |&b| b == byte, split)
        } else {
            let mut equal = [false; 256];
            equal[usize::from(byte)] = true;
            for &EqualityPair { first, second } in pairs {
                if first == byte {
                    equal[usize::from(second)] = true;
                } else if second == byte {
                    equal[usize::from(first)] = true;
                }
            }
            one_byte_position(long, |&b| equal[usize::from(b)], split)
        };
        Self {
            gap,
            diagonal: if matched { EditOperation::Match } else { EditOperation::Mismatch },
            before,
            after: long.len() - 1 - before,
        }
    }

    pub(crate) fn score(&self) -> usize {
        self.before + self.after + usize::from(self.diagonal == EditOperation::Mismatch)
    }

    pub(crate) fn run_count(&self) -> usize {
        1 + usize::from(self.before != 0) + usize::from(self.after != 0)
    }

    pub(crate) fn runs(self) -> impl Iterator<Item = Run> {
        [
            Run { op: self.gap, len: self.before },
            Run { op: self.diagonal, len: 1 },
            Run { op: self.gap, len: self.after },
        ].into_iter().filter(|run| run.len != 0)
    }
}

fn one_byte_position(long: &[u8], equal: impl Fn(&u8) -> bool, split: bool) -> (usize, bool) {
    // Backward traceback prefers gaps over a diagonal, leaving the first
    // match (or the first byte on a miss) inside a leaf.
    if !split {
        return long.iter().position(equal).map_or((0, false), |i| (i, true));
    }

    // Hirschberg tries row 0 before row 1: the right half wins if it has a
    // match, or if both halves miss. Locate that leaf with the last match,
    // then search its unvisited prefix for the first match within the leaf.
    let last = long.iter().rposition(&equal);
    let pivot = last.unwrap_or(long.len() - 1);
    let (mut start, mut len) = (0, long.len());
    while !trace_fits(1, len) {
        let middle = len / 2;
        if pivot - start < middle {
            len = middle;
        } else {
            start += middle;
            len -= middle;
        }
    }
    match last {
        Some(i) => (start + long[start..i].iter().position(equal).unwrap_or(i - start), true),
        None => (start, false),
    }
}

/// Score and recover one global path. Never called for a locations or distance task.
/// `prefix` prepends already matched bytes without shifting the recovered output.
pub(crate) fn find_global<O: Output>(
    query: &[u8], target: &[u8], pairs: &[EqualityPair], limit: Option<usize>, prefix: usize,
) -> Result<Option<(usize, O)>, AlignError> {
    let (m, n) = (query.len(), target.len());
    let buffer = |d| {
        let mut output = O::buffer(prefix + m, prefix + n, d)?;
        if prefix != 0 {
            output.append_run(0, EditOperation::Match, prefix)?;
        }
        Ok::<O, AlignError>(output)
    };
    if m.min(n) == 1 {
        let split = O::LEGACY_SPLIT_TIES && m == 1 && !trace_fits(1, n);
        let single = OneByte::find(query, target, pairs, split);
        let score = single.score();
        if limit.is_some_and(|k| score > k) {
            return Ok(None);
        }
        if score == 0 {
            return Ok(Some((0, O::single(EditOperation::Match, prefix + m))));
        }
        let mut output = buffer(score)?;
        for Run { op, len } in single.runs() {
            output.append_run(0, op, len)?;
        }
        return Ok(Some((score, output)));
    }
    // Distance zero has no DP replay to save. Literal equality suffices even
    // with extra equality pairs; differing bytes must still go through Myers.
    if query == target {
        return Ok(Some((0, O::single(EditOperation::Match, prefix + m))));
    }
    let profile = Profile::new(Sequence::new(query), pairs)?;
    let seq = Sequence::new(target);
    if trace_fits(profile.blocks, n) {
        if let Some(k) = limit {
            // A concrete diagonal path (plus the length gap) certifies a hit.
            // A failed certificate says nothing about an alignment with indels
            // or extra equalities: score first, then reconstruct *every* hit.
            if k == 0 || (k < m.max(n) && diagonal_cost(query, target) > k) {
                let Some(d) = distance_global(&profile, seq, k, &mut Vec::new()) else { return Ok(None) };
                let mut output = buffer(d)?;
                append_alignment(query, target, pairs, d, Some(profile), &mut output)?;
                return Ok(Some((d, output)));
            }
        }
        let mut trace = Trace::new(profile.blocks, n)?;
        let distance = search_distance(m.max(n), limit, |k| {
            if m.abs_diff(n) > k {
                return None;
            }
            // Failed attempts overwrite this storage, never append to an old pass.
            trace.blocks.clear();
            trace.columns.clear();
            if profile.blocks == 1 {
                let score = trace_single(&profile, seq, &mut trace) as usize;
                (score <= k).then_some(score)
            } else {
                global(&profile, seq, k, n, Some(&mut trace))?.distance(m, k)
            }
        });
        let Some(d) = distance else { return Ok(None) };
        drop(profile);
        if d == 0 {
            return Ok(Some((0, O::single(EditOperation::Match, prefix + m))));
        }
        let mut output = buffer(d)?;
        append_traceback(m, n, d, &trace, &mut output)?;
        return Ok(Some((d, output)));
    }

    let mut cuts = cuts::Cuts::new(m, n);
    let distance = search_distance(m.max(n), limit, |k| cuts.score(&profile, seq, k));
    let Some(d) = distance else { return Ok(None) };
    if d == 0 {
        return Ok(Some((0, O::single(EditOperation::Match, prefix + m))));
    }
    let mut output = buffer(d)?;
    if cuts.is_empty() {
        append_alignment(query, target, pairs, d, Some(profile), &mut output)?;
    } else {
        drop(profile);
        cuts.resolve(query, target, pairs, d)?;
        cuts.append(query, target, pairs, d, (0, 0, 0), &mut output)?;
    }
    Ok(Some((d, output)))
}

fn diagonal_cost(query: &[u8], target: &[u8]) -> usize {
    query.len().abs_diff(target.len()) + query.iter().zip(target).filter(|(a, b)| a != b).count()
}

/// Reconstructs an optimal global path for the selected query/target span.
///
/// `distance` is the exact distance already found by the distance kernel. All
/// recursive children retain that invariant; no unbounded scalar DP is needed.
/// `profile`, when present, is the caller's forward profile of `query` with `pairs`.
pub(crate) fn reconstruct<O: Output>(
    query: &[u8],
    target: &[u8],
    pairs: &[EqualityPair],
    distance: usize,
    profile: Option<Profile>,
) -> Result<O, AlignError> {
    let mut output = O::buffer(query.len(), target.len(), distance)?;
    append_alignment(query, target, pairs, distance, profile, &mut output)?;
    Ok(output)
}

fn operation_buffer(m: usize, n: usize, distance: usize) -> Result<Vec<EditOperation>, AlignError> {
    let difference = m.abs_diff(n);
    let excess = distance
        .checked_sub(difference)
        .expect("an optimal distance is at least the length difference");
    // A path has at most floor((m + n + distance) / 2) operations. This
    // equivalent expression avoids overflowing m + n and reserves only once.
    let capacity = m.max(n)
        .checked_add(excess / 2)
        .ok_or(AlignError::InputTooLong)?;
    checked_allocation::<EditOperation>(capacity)?;
    Ok(Vec::with_capacity(capacity))
}

fn checked_allocation<T>(count: usize) -> Result<(), AlignError> {
    match count.checked_mul(size_of::<T>()) {
        Some(bytes) if bytes <= isize::MAX as usize => Ok(()),
        _ => Err(AlignError::InputTooLong),
    }
}

/// `profile`, when present, is the forward profile of exactly `query` and
/// `pairs`. It is consumed here, so it never survives into recursive children.
fn append_alignment<O: Output>(
    query: &[u8],
    target: &[u8],
    pairs: &[EqualityPair],
    distance: usize,
    profile: Option<Profile>,
    output: &mut O,
) -> Result<(), AlignError> {
    if query.is_empty() {
        return output.append_run(0, EditOperation::Delete, target.len());
    }
    if target.is_empty() {
        return output.append_run(0, EditOperation::Insert, query.len());
    }
    if distance == 0 {
        debug_assert_eq!(query.len(), target.len());
        return output.append_run(0, EditOperation::Match, query.len());
    }

    let blocks = (query.len() - 1) / WORD_BITS + 1;
    // A one-column matrix cannot be divided further. Its trace is linear in
    // the query even when it exceeds the usual small-trace memory threshold.
    if target.len() == 1 || trace_fits(blocks, target.len()) {
        let profile = match profile {
            Some(profile) => profile,
            None => Profile::new(Sequence::new(query), pairs)?,
        };
        debug_assert_eq!(profile.len, query.len());
        let mut trace = Trace::new(blocks, target.len())?;
        let score = if blocks == 1 {
            Some(trace_single(&profile, Sequence::new(target), &mut trace))
        } else {
            global(
                &profile,
                Sequence::new(target),
                distance,
                target.len(),
                Some(&mut trace),
            )
            .expect("a known optimal path must remain in the global band")
            .value(query.len())
        };
        debug_assert_eq!(
            score,
            Some(isize::try_from(distance).map_err(|_| AlignError::InputTooLong)?)
        );
        drop(profile);
        return append_traceback(query.len(), target.len(), distance, &trace, output);
    }

    let middle = target.len() / 2;
    let (row, left_distance) = split(query, target, pairs, distance, profile)?;
    // No profiles or columns survive into either child.
    append_alignment(
        &query[..row],
        &target[..middle],
        pairs,
        left_distance,
        None,
        output,
    )?;
    append_alignment(
        &query[row..],
        &target[middle..],
        pairs,
        distance - left_distance,
        None,
        output,
    )
}

fn trace_fits(blocks: usize, columns: usize) -> bool {
    trace_fits_budget(blocks, columns, TRACE_BYTES)
}

fn trace_fits_budget(blocks: usize, columns: usize, budget: usize) -> bool {
    // Trace stores compact active blocks, but using the full-width upper bound
    // keeps its allocation bounded without depending on the eventual band.
    // Rust's isize scores/ranges make this threshold differ slightly from C.
    size_of::<crate::kernel::Block>()
        .checked_mul(blocks)
        .and_then(|bytes| bytes.checked_add(size_of::<(usize, std::ops::Range<usize>)>()))
        .and_then(|bytes| bytes.checked_mul(columns))
        .is_some_and(|bytes| bytes < budget)
}

#[inline(never)]
fn append_traceback<O: Output>(
    mut row: usize,
    mut column: usize,
    distance: usize,
    trace: &Trace,
    output: &mut O,
) -> Result<(), AlignError> {
    let mut score = isize::try_from(distance).map_err(|_| AlignError::InputTooLong)?;
    let (first, tail) = output.begin_leaf();
    while row != 0 && column != 0 {
        // The current cell is retained. Its positive vertical delta tests an
        // insertion without decoding another score. Other predecessors may
        // be absent from the band. Preserve insertion/deletion/diagonal order.
        let lower_score = score.checked_sub(1);
        if lower_score.is_some() && trace.can_insert(row, column) {
            output.append_step(first, EditOperation::Insert)?;
            row -= 1;
            score -= 1;
        } else {
            let (left, diagonal) = trace.scores(row, column - 1);
            if lower_score.is_some() && left == lower_score {
                output.append_step(first, EditOperation::Delete)?;
                column -= 1;
                score -= 1;
            } else {
                let diagonal =
                    diagonal.expect("an optimal traceback must have a retained predecessor");
                debug_assert!(diagonal == score || Some(diagonal) == lower_score);
                output.append_step(first, if diagonal == score {
                    EditOperation::Match
                } else {
                    EditOperation::Mismatch
                })?;
                row -= 1;
                column -= 1;
                score = diagonal;
            }
        }
    }
    debug_assert_eq!(score as usize, row.max(column));
    output.append_run(first, EditOperation::Insert, row)?;
    output.append_run(first, EditOperation::Delete, column)?;
    // Only this child's appended suffix is backward.
    output.finish_leaf(first, tail)
}

fn split(
    query: &[u8],
    target: &[u8],
    pairs: &[EqualityPair],
    distance: usize,
    forward: Option<Profile>,
) -> Result<(usize, usize), AlignError> {
    let profile = match forward {
        Some(profile) => profile,
        None => Profile::new(Sequence::new(query), pairs)?,
    };
    let middle = target.len() / 2;
    // Keep the full target length: pruning depends on the far endpoint.
    let left = global(&profile, Sequence::new(target), distance, middle, None)
        .expect("a known optimal path must cross the forward split column");
    drop(profile);
    let reverse_profile = Profile::new(Sequence::reversed(query), pairs)?;
    let right = global(
        &reverse_profile,
        Sequence::reversed(target),
        distance,
        target.len() - middle,
        None,
    )
    .expect("a known optimal path must cross the reverse split column");
    drop(reverse_profile);

    Ok(split_rows(query.len(), target.len(), middle, distance, left.band(), right.band())
        .next().expect("an optimum crosses the retained split columns"))
}

fn split_rows<'a>(
    m: usize, n: usize, middle: usize, distance: usize, left: Band<'a>, right: Band<'a>,
) -> impl Iterator<Item = (usize, usize)> + 'a {
    let rows = |band: Band<'_>| {
        (band.first * WORD_BITS + 1, ((band.first + band.blocks.len()) * WORD_BITS).min(m))
    };
    let (lf, ll) = rows(left);
    let (rf, rl) = rows(right);
    let first = lf.max(m - rl).max(1);
    let last = ll.min(m - rf).min(m - 1);
    // Original Hirschberg tie order: increasing interior rows, then 0, then m.
    (first..=last).chain([0, m]).filter_map(move |row| {
        let f = if row == 0 { middle } else { left.value(row)? as usize };
        let r = if row == m { n - middle } else { right.value(m - row)? as usize };
        (f <= distance && r == distance - f).then_some((row, f))
    })
}

