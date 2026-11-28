//! Score-only global Myers passes with a single word of state per column.

use crate::kernel::horizontal;

pub(super) const WIDTH: usize = u64::BITS as usize;
const STEP: usize = 8;
pub(super) const BAND: usize = WIDTH - STEP;

struct Column {
    positive: u64,
    negative: u64,
}

impl Column {
    const INITIAL: Self = Self { positive: u64::MAX, negative: 0 };

    #[inline]
    fn advance(&mut self, eq: u64) -> (u64, u64) {
        let xv = eq | self.negative;
        let (ph, mh) = horizontal(eq, self.positive, self.negative);
        // The top boundary advances by insertion. Complementing Ph before
        // shifting encodes that positive low bit without an extra OR.
        let not_ph = !ph << 1;
        self.positive = (mh << 1) | (not_ph & !xv);
        self.negative = xv & !not_ph;
        (ph, mh)
    }

    #[inline]
    fn value(&self, top: usize, rows: usize) -> usize {
        let used = u64::MAX.checked_shr((WIDTH - rows) as u32).unwrap_or(0);
        top + (self.positive & used).count_ones() as usize
            - (self.negative & used).count_ones() as usize
    }
}

#[inline]
fn masks(query: &[u8]) -> [u64; 256] {
    let mut masks = [0u64; 256];
    // A run of equal bytes would otherwise serialize eight read/modify/writes
    // to the same mask. Mixed groups still use ordinary byte indexing.
    let mut chunks = query.chunks_exact(8);
    let mut base = 0;
    for chunk in chunks.by_ref() {
        let bytes: &[u8; 8] = chunk.try_into().unwrap();
        if *bytes == [bytes[0]; 8] {
            masks[bytes[0] as usize] |= 0xff << base;
        } else {
            for (i, &byte) in bytes.iter().enumerate() {
                masks[byte as usize] |= 1 << (base + i);
            }
        }
        base += 8;
    }
    for (i, &byte) in chunks.remainder().iter().enumerate() {
        masks[byte as usize] |= 1 << (base + i);
    }
    masks
}

/// Every global path crosses the middle target column at some query row i.
/// Minimize prefix(i) + suffix(rows-i). The two half-column recurrences
/// are independent, cutting the serial dependency through the target in half.
pub(super) fn short(query: &[u8], target: &[u8]) -> usize {
    debug_assert!(!query.is_empty() && query.len() <= WIDTH && query.len() <= target.len());
    if let [byte] = query {
        return target.len() - usize::from(target.contains(byte));
    }
    let masks = masks(query);
    let rows = query.len();
    // On tiny targets, joining two columns costs more than advancing one.
    if target.len() <= 12 {
        let mut column = Column::INITIAL;
        for &byte in target {
            column.advance(masks[byte as usize]);
        }
        return column.value(target.len(), rows);
    }
    let (prefix, suffix) = target.split_at(target.len() / 2);
    let reverse = |bits: u64| bits.reverse_bits() >> (WIDTH - rows);
    let (mut left, mut right) = (Column::INITIAL, Column::INITIAL);
    let mut back = suffix.iter().rev();
    for (&a, &b) in prefix.iter().zip(back.by_ref()) {
        left.advance(masks[a as usize]);
        right.advance(reverse(masks[b as usize]));
    }
    if let Some(&b) = back.next() {
        right.advance(reverse(masks[b as usize]));
    }

    let cost = |i| left.value(prefix.len(), i) + right.value(suffix.len(), rows - i);
    let gap = target.len() - rows;
    let slack = (cost(rows / 2) - gap) / 2;
    // The midpoint split supplies a feasible script. Length gaps alone
    // exclude rows outside this interval from improving it.
    let lo = prefix.len().saturating_sub(gap + slack);
    let hi = rows.min(prefix.len() + slack);
    let mut score = cost(lo);
    let mut best = score;
    let (mut lp, mut lm) = (left.positive >> lo, left.negative >> lo);
    let (mut rp, mut rm) = (reverse(right.positive) >> lo, reverse(right.negative) >> lo);
    for _ in lo..hi {
        // Moving the split down adds a left row and removes a right row.
        score = score + ((lp & 1) + (rm & 1)) as usize - ((lm & 1) + (rp & 1)) as usize;
        best = best.min(score);
        (lp, lm, rp, rm) = (lp >> 1, lm >> 1, rp >> 1, rm >> 1);
    }
    best
}

/// Unlike the split pass, this tracks the last row so a tight miss can stop
/// early. The remaining columns can lower that row's score by at most one each.
pub(super) fn short_at_most(query: &[u8], target: &[u8], k: usize) -> Option<usize> {
    if k >= target.len() {
        let score = short(query, target);
        return (score <= k).then_some(score);
    }
    let masks = masks(query);
    let mut col = Column::INITIAL;
    let high = 1 << (query.len() - 1);
    let mut lower = query.len() as isize - target.len() as isize;
    for &byte in target {
        let (ph, mh) = col.advance(masks[byte as usize]);
        lower += 1 + isize::from(ph & high != 0) - isize::from(mh & high != 0);
        if lower > k as isize {
            return None;
        }
    }
    Some(lower as usize)
}

/// A moving window of 64 query rows suffices for all paths of cost <= BAND:
/// on diagonal z = j-i any path needs at least |z| + |gap-z| edits to pass
/// through it and still reach the endpoint. Above BAND the window only
/// supplies an upper bound, so do not return that value as a score.
pub(super) fn sliding(query: &[u8], target: &[u8]) -> Option<usize> {
    debug_assert!(query.len() > WIDTH && query.len() <= target.len());
    let gap = target.len() - query.len();
    debug_assert!(gap <= BAND);
    let mut masks = masks(&query[..WIDTH]);
    let mut col = Column::INITIAL;
    // Finish each tile when the lowest viable diagonal reaches its bottom
    // row. After the slide, its top boundary is above every viable row of
    // the next column. No path costing <= BAND needs the discarded rows.
    let below = (BAND - gap) / 2;
    let (mut base, mut start, mut top) = (0, 0, 0);
    loop {
        let remaining = query.len() - base;
        let end = if remaining == WIDTH { target.len() } else { base + WIDTH - below };
        let rotation = (base % WIDTH) as u32;
        for &byte in &target[start..end] {
            col.advance(masks[byte as usize].rotate_right(rotation));
        }
        top += end - start;
        if remaining == WIDTH {
            let score = col.value(top, WIDTH);
            return (score <= BAND).then_some(score);
        }
        let shift = STEP.min(remaining - WIDTH);
        top = col.value(top, shift);
        // Carry the overlapping rows forward; initialize new bottom rows by
        // deletion from the old bottom. Mask slots wrap in absolute query order.
        col.positive = !(!col.positive >> shift);
        col.negative >>= shift;
        let mut bit = 1u64.rotate_left(rotation);
        for (&old, &new) in query[base..base + shift].iter()
            .zip(&query[base + WIDTH..base + WIDTH + shift])
        {
            masks[old as usize] &= !bit;
            masks[new as usize] |= bit;
            bit = bit.rotate_left(1);
        }
        base += shift;
        start = end;
    }
}
