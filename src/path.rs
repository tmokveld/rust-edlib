use std::mem::size_of;

use crate::aligner::{AlignError, EditOperation, EqualityPair};
use crate::kernel::{Column, Profile, Sequence, Trace, global, trace_single};

const WORD_BITS: usize = u64::BITS as usize;
const TRACE_BYTES: usize = 1024 * 1024;

/// Reconstructs an optimal global path for the selected query/target span.
///
/// `distance` is the exact distance already found by the distance kernel. All
/// recursive children retain that invariant; no unbounded scalar DP is needed.
/// `profile` is the caller's forward profile of `query` with `pairs`.
pub(crate) fn reconstruct(
    query: &[u8],
    target: &[u8],
    pairs: &[EqualityPair],
    distance: usize,
    profile: Profile,
) -> Result<Vec<EditOperation>, AlignError> {
    let difference = query.len().abs_diff(target.len());
    let excess = distance
        .checked_sub(difference)
        .expect("an optimal distance is at least the length difference");
    // A path has at most floor((m + n + distance) / 2) operations. This
    // equivalent expression avoids overflowing m + n and reserves only once.
    let capacity = query
        .len()
        .max(target.len())
        .checked_add(excess / 2)
        .ok_or(AlignError::InputTooLong)?;
    checked_allocation::<EditOperation>(capacity)?;
    let mut output = Vec::with_capacity(capacity);
    append_alignment(query, target, pairs, distance, Some(profile), &mut output)?;
    Ok(output)
}

fn checked_allocation<T>(count: usize) -> Result<(), AlignError> {
    match count.checked_mul(size_of::<T>()) {
        Some(bytes) if bytes <= isize::MAX as usize => Ok(()),
        _ => Err(AlignError::InputTooLong),
    }
}

/// `profile`, when present, is the forward profile of exactly `query` and
/// `pairs`. It is consumed here, so it never survives into recursive children.
fn append_alignment(
    query: &[u8],
    target: &[u8],
    pairs: &[EqualityPair],
    distance: usize,
    profile: Option<Profile>,
    output: &mut Vec<EditOperation>,
) -> Result<(), AlignError> {
    if query.is_empty() {
        return append_run(output, EditOperation::Delete, target.len());
    }
    if target.is_empty() {
        return append_run(output, EditOperation::Insert, query.len());
    }
    if distance == 0 {
        debug_assert_eq!(query.len(), target.len());
        return append_run(output, EditOperation::Match, query.len());
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
    let (row, left_distance, right_distance) =
        split(query, target, pairs, distance, middle, profile)?;
    // split() owns all profiles and columns, so none survive into recursion.
    // Both target halves are strictly shorter: stack depth is O(log n).
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
        right_distance,
        None,
        output,
    )
}

fn trace_fits(blocks: usize, columns: usize) -> bool {
    // Trace stores compact active blocks, but using the full-width upper bound
    // keeps its allocation bounded without depending on the eventual band.
    // Rust's isize scores/ranges make this threshold differ slightly from C.
    size_of::<crate::kernel::Block>()
        .checked_mul(blocks)
        .and_then(|bytes| bytes.checked_add(size_of::<(usize, std::ops::Range<usize>)>()))
        .and_then(|bytes| bytes.checked_mul(columns))
        .is_some_and(|bytes| bytes < TRACE_BYTES)
}

fn append_run(
    output: &mut Vec<EditOperation>,
    operation: EditOperation,
    count: usize,
) -> Result<(), AlignError> {
    let end = output
        .len()
        .checked_add(count)
        .ok_or(AlignError::InputTooLong)?;
    debug_assert!(end <= output.capacity());
    output.resize(end, operation);
    Ok(())
}

#[inline(never)]
fn append_traceback(
    mut row: usize,
    mut column: usize,
    distance: usize,
    trace: &Trace,
    output: &mut Vec<EditOperation>,
) -> Result<(), AlignError> {
    let start = output.len();
    let mut score = isize::try_from(distance).map_err(|_| AlignError::InputTooLong)?;
    while row != 0 && column != 0 {
        // The current cell is retained. Its positive vertical delta tests an
        // insertion without decoding another score. Other predecessors may
        // be absent from the band. Preserve insertion/deletion/diagonal order.
        let lower_score = score.checked_sub(1);
        if lower_score.is_some() && trace.can_insert(row, column) {
            output.push(EditOperation::Insert);
            row -= 1;
            score -= 1;
        } else {
            let (left, diagonal) = trace.scores(row, column - 1);
            if lower_score.is_some() && left == lower_score {
                output.push(EditOperation::Delete);
                column -= 1;
                score -= 1;
            } else {
                let diagonal =
                    diagonal.expect("an optimal traceback must have a retained predecessor");
                debug_assert!(diagonal == score || Some(diagonal) == lower_score);
                output.push(if diagonal == score {
                    EditOperation::Match
                } else {
                    EditOperation::Mismatch
                });
                row -= 1;
                column -= 1;
                score = diagonal;
            }
        }
    }
    debug_assert_eq!(score as usize, row.max(column));
    append_run(output, EditOperation::Insert, row)?;
    append_run(output, EditOperation::Delete, column)?;
    // Only this child's appended suffix is backward; already emitted siblings
    // remain in forward order. No intermediate alignment vectors are allocated.
    output[start..].reverse();
    Ok(())
}

fn split(
    query: &[u8],
    target: &[u8],
    pairs: &[EqualityPair],
    distance: usize,
    middle: usize,
    forward: Option<Profile>,
) -> Result<(usize, usize, usize), AlignError> {
    let profile = match forward {
        Some(profile) => profile,
        None => Profile::new(Sequence::new(query), pairs)?,
    };
    debug_assert_eq!(profile.len, query.len());
    // Keep the FULL target length in both calls. Truncating its borrowed view
    // at the split would change the endpoint-aware global band's geometry.
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

    let length = query.len();
    let (left_first, left_last) = column_rows(&left, length);
    let (right_first, right_last) = column_rows(&right, length);
    let first = left_first.max(length - right_last).max(1);
    let last = left_last.min(length - right_first).min(length - 1);
    let wanted = isize::try_from(distance).map_err(|_| AlignError::InputTooLong)?;
    // Native Hirschberg visits interior rows in increasing query order before
    // either boundary. Intersect the bands first instead of scanning all rows.
    for row in first..=last {
        if let (Some(forward), Some(reverse)) = (left.value(row), right.value(length - row)) {
            if forward <= wanted && reverse == wanted - forward {
                return Ok((row, forward as usize, reverse as usize));
            }
        }
    }
    if let Some(reverse) = right.value(length) {
        if middle <= distance && reverse as usize == distance - middle {
            return Ok((0, middle, reverse as usize));
        }
    }
    let remaining = target.len() - middle;
    if let Some(forward) = left.value(length) {
        if remaining <= distance && forward as usize == distance - remaining {
            return Ok((length, forward as usize, remaining));
        }
    }
    unreachable!("a known optimal path must cross a retained split cell")
}

fn column_rows(column: &Column, length: usize) -> (usize, usize) {
    // first/last are absolute block indices. The final block may have padding;
    // saturating the padded endpoint is harmless before clipping to length.
    let first = column.first * WORD_BITS + 1;
    let last = column
        .last
        .saturating_add(1)
        .saturating_mul(WORD_BITS)
        .min(length);
    (first, last)
}
