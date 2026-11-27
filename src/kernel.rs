use crate::aligner::{AlignError, EqualityPair};
use std::{mem::size_of, ops::Range};

pub(crate) fn check_buffer<T>(len: usize) -> Result<(), AlignError> {
    len.checked_mul(size_of::<T>())
        .filter(|&n| n <= isize::MAX as usize)
        .map(|_| ())
        .ok_or(AlignError::InputTooLong)
}

/// Initial budget above the mandatory gap for unbounded score-only searches.
pub(crate) fn initial_slack(gap: usize, overlap: usize) -> usize {
    (gap.min(overlap) / 4).max(64)
}

/// Reversed traversal without copying input bytes.
#[derive(Clone, Copy)]
pub(crate) struct Sequence<'a> {
    bytes: &'a [u8],
    reverse: bool,
}
impl<'a> Sequence<'a> {
    pub fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            reverse: false,
        }
    }
    pub fn reversed(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            reverse: true,
        }
    }
    pub fn len(self) -> usize {
        self.bytes.len()
    }
    #[inline]
    pub fn at(self, i: usize) -> u8 {
        self.bytes[if self.reverse {
            self.bytes.len() - 1 - i
        } else {
            i
        }]
    }
}

/// Only query symbols and equality partners need rows. Unknown target symbols
/// share row zero. No target scan, transformed input, or 256 x m initialization.
pub(crate) struct Profile {
    pub len: usize,
    pub blocks: usize,
    // Row IDs need nine bits, irrespective of query length. A compact index
    // avoids moving/clearing 2 KiB per call; resolve its word offset once per
    // target column, outside the block loop.
    indices: [u16; 256],
    words: Vec<u64>,
}

impl Profile {
    // Keep construction in the caller: out-of-line Result returns copied the
    // profile index twice in short alignments (see results/profile-short).
    #[inline(always)]
    pub fn new(query: Sequence<'_>, pairs: &[EqualityPair]) -> Result<Self, AlignError> {
        let blocks = query.len().div_ceil(64);
        let mut indices = [0u16; 256];
        let mut rows = 1;
        let mut add = |byte: u8| {
            if indices[byte as usize] == 0 {
                indices[byte as usize] = rows as u16;
                rows += 1;
            }
        };
        // Precheck makes the per-symbol offset multiplication safe.
        check_buffer::<u64>(blocks.checked_mul(257).ok_or(AlignError::InputTooLong)?)?;
        for i in 0..query.len() {
            add(query.at(i));
        }
        for pair in pairs {
            add(pair.first);
            add(pair.second);
        }
        let mut words = vec![0; rows * blocks];
        // Row IDs below `rows` need `bits` bits. Tiny alphabets use bit planes
        // in registers; others scatter. The scatter's read-modify-write chains
        // serialize through memory whenever a symbol repeats closely.
        let bits = usize::BITS - (rows - 1).leading_zeros();
        match bits {
            1 => Self::planes::<1>(&mut words, &indices, blocks, rows, query),
            2 => Self::planes::<2>(&mut words, &indices, blocks, rows, query),
            3 => Self::planes::<3>(&mut words, &indices, blocks, rows, query),
            _ => {
                for i in 0..query.len() {
                    words[indices[query.at(i) as usize] as usize * blocks + i / 64] |=
                        1u64 << (i % 64);
                }
            }
        }
        if !pairs.is_empty() {
            Self::add_equalities(&mut words, &indices, blocks, query, pairs);
        }
        // A suffix of wildcard padding cannot influence any real query row.
        let used = query.len() % 64;
        if used != 0 {
            let padding = u64::MAX << used;
            for row in words.chunks_exact_mut(blocks) {
                row[blocks - 1] |= padding;
            }
        }
        Ok(Self {
            len: query.len(),
            blocks,
            indices,
            words,
        })
    }

    /// Bit-sliced construction for row IDs of at most `P` bits. Plane `j`
    /// collects bit `j` of every position's row ID in one 64-position block;
    /// a row's word is the intersection of planes (or complements) matching
    /// its ID. Every row ID is nonzero, so positions past the query end, whose
    /// planes are zero, never enter a row. Exactly the scatter loop's words.
    #[inline(always)]
    fn planes<const P: usize>(
        words: &mut [u64],
        indices: &[u16; 256],
        blocks: usize,
        rows: usize,
        query: Sequence<'_>,
    ) {
        for b in 0..blocks {
            let start = b * 64;
            let mut plane = [0u64; P];
            for i in start..query.len().min(start + 64) {
                let id = u64::from(indices[query.at(i) as usize]);
                for (j, word) in plane.iter_mut().enumerate() {
                    *word |= ((id >> j) & 1) << (i - start);
                }
            }
            for row in 1..rows {
                let mut mask = u64::MAX;
                for (j, &word) in plane.iter().enumerate() {
                    mask &= if (row >> j) & 1 != 0 { word } else { !word };
                }
                words[row * blocks + b] = mask;
            }
        }
    }

    // Kept out of the common path so ordinary equality needs no large table.
    fn add_equalities(
        words: &mut [u64],
        indices: &[u16; 256],
        blocks: usize,
        query: Sequence<'_>,
        pairs: &[EqualityPair],
    ) {
        let mut neighbors = [[0u64; 4]; 256];
        for pair in pairs {
            let (a, b) = (pair.first as usize, pair.second as usize);
            neighbors[a][b / 64] |= 1 << (b % 64);
            neighbors[b][a / 64] |= 1 << (a % 64);
        }
        // Read the original query, never an already-extended profile: OR-ing
        // extended rows would accidentally make equalities transitive.
        for i in 0..query.len() {
            for (group, &bits) in neighbors[query.at(i) as usize].iter().enumerate() {
                let mut bits = bits;
                while bits != 0 {
                    let symbol = group * 64 + bits.trailing_zeros() as usize;
                    words[indices[symbol] as usize * blocks + i / 64] |= 1 << (i % 64);
                    bits &= bits - 1;
                }
            }
        }
    }

    #[inline]
    fn row(&self, symbol: u8) -> &[u64] {
        let offset = self.indices[symbol as usize] as usize * self.blocks;
        &self.words[offset..offset + self.blocks]
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Block {
    pub p: u64,
    pub m: u64,
    pub score: isize,
}

impl Block {
    pub(crate) fn new(score: isize) -> Self {
        Self { p: u64::MAX, m: 0, score }
    }

    #[inline]
    pub(crate) fn advance(&mut self, eq: u64, carry: Deltas, high: usize) -> Deltas {
        let xv = eq | self.m;
        let (ph, mh) = horizontal(eq | carry.m, self.p, self.m);
        let next = Deltas { p: (ph >> high) & 1, m: (mh >> high) & 1 };
        let ph = (ph << 1) | carry.p;
        self.p = (mh << 1) | carry.m | !(xv | ph);
        self.m = ph & xv;
        self.score += next.value();
        next
    }

    /// Vertical differences encode D[r,c]-D[r-1,c]. Bits above the requested
    /// row must be subtracted from the bottom score; popcount avoids 64 steps.
    #[inline]
    pub fn value(self, row: usize) -> isize {
        let mask = if row == 63 { 0 } else { u64::MAX << (row + 1) };
        self.score - (self.p & mask).count_ones() as isize + (self.m & mask).count_ones() as isize
    }
    pub(crate) fn all_larger(self, k: isize, width: usize) -> bool {
        let mut value = self.score;
        for i in (0..width).rev() {
            if value <= k {
                return false;
            }
            value -= ((self.p >> i) & 1) as isize;
            value += ((self.m >> i) & 1) as isize;
        }
        true
    }
}

/// Split +/- horizontal carries avoid decoding/re-encoding between words.
#[derive(Clone, Copy)]
pub(crate) struct Deltas {
    p: u64,
    pub(crate) m: u64,
}

impl Deltas {
    // D[0,c] = c when the leading target end is not free.
    pub(crate) const ANCHORED: Self = Self { p: 1, m: 0 };

    #[inline]
    pub(crate) fn value(self) -> isize {
        self.p as isize - self.m as isize
    }
}

/// Myers' horizontal deltas without materializing `Xh = (sum ^ pv) | eq`,
/// where `sum = (eq & pv) + pv`. Exact bitwise identities give
/// `Xh | pv = sum | eq | pv` and `pv & Xh = pv & (eq | !sum)`, removing the
/// XOR and one OR from the loop-carried dependency chain. Returns `(Ph, Mh)`.
#[inline(always)]
fn horizontal(eq: u64, pv: u64, mv: u64) -> (u64, u64) {
    // Overflow is the discarded carry beyond this word, not a Rust overflow.
    let sum = (eq & pv).wrapping_add(pv);
    (mv | !(sum | (eq | pv)), pv & (eq | !sum))
}

/// Advance complete 64-row blocks. In Myers' multiword sum `(Eq & Pv) + Pv`,
/// the carry out of one half is exactly its high `Mh` bit. A `u128` addition
/// delivers it directly to the next half; bits 63 and 127 still give the
/// two separate bottom-score deltas needed by pruning and traceback.
#[inline]
pub(crate) fn advance_words(blocks: &mut [Block], row: &[u64], carry: &mut Deltas) -> isize {
    debug_assert_eq!(blocks.len(), row.len());
    let join = |a: u64, b: u64| u128::from(a) | (u128::from(b) << 64);
    let mut pairs = blocks.chunks_exact_mut(2);
    let mut eqs = row.chunks_exact(2);
    let mut last = 0;
    for (b, e) in pairs.by_ref().zip(eqs.by_ref()) {
        let pv = join(b[0].p, b[1].p);
        let mv = join(b[0].m, b[1].m);
        let eq = join(e[0], e[1]);
        let xv = eq | mv;
        let eq = eq | u128::from(carry.m);
        let sum = (eq & pv).wrapping_add(pv);
        let ph = mv | !(sum | (eq | pv));
        let mh = pv & (eq | !sum);
        b[0].score += (ph as u64 >> 63) as isize - (mh as u64 >> 63) as isize;
        let next = Deltas { p: (ph >> 127) as u64, m: (mh >> 127) as u64 };
        b[1].score += next.value();
        last = b[1].score;
        let ph = (ph << 1) | u128::from(carry.p);
        let pv = (mh << 1) | u128::from(carry.m) | !(xv | ph);
        let mv = ph & xv;
        (b[0].p, b[1].p) = (pv as u64, (pv >> 64) as u64);
        (b[0].m, b[1].m) = (mv as u64, (mv >> 64) as u64);
        *carry = next;
    }
    if let ([block], [eq]) = (pairs.into_remainder(), eqs.remainder()) {
        *carry = block.advance(*eq, *carry, 63);
        last = block.score;
    }
    last
}

#[inline]
fn advance_single<const FREE: bool>(
    eq: u64,
    high: u64,
    pv: &mut u64,
    mv: &mut u64,
    score: &mut isize,
) {
    let xv = eq | *mv;
    let (ph, mh) = horizontal(eq, *pv, *mv);
    *score += isize::from(ph & high != 0) - isize::from(mh & high != 0);
    let ph = (ph << 1) | u64::from(!FREE);
    *pv = (mh << 1) | !(xv | ph);
    *mv = ph & xv;
}

pub(crate) struct Column {
    pub blocks: Vec<Block>,
    pub first: usize,
    pub last: usize,
}

impl Column {
    pub fn value(&self, row: usize) -> Option<isize> {
        let b = row.checked_sub(1)? / 64;
        (self.first..=self.last)
            .contains(&b)
            .then(|| self.blocks[b].value((row - 1) % 64))
    }
}

/// Compact saved active blocks. Only small subproblems create a Trace; larger
/// ones retain two columns and recurse. Capacity is checked and reserved once.
pub(crate) struct Trace {
    pub blocks: Vec<Block>,
    pub columns: Vec<(usize, Range<usize>)>,
}

impl Trace {
    pub fn new(max_blocks: usize, columns: usize) -> Result<Self, AlignError> {
        let count = max_blocks
            .checked_mul(columns)
            .ok_or(AlignError::InputTooLong)?;
        check_buffer::<Block>(count)?;
        check_buffer::<(usize, Range<usize>)>(columns)?;
        Ok(Self {
            blocks: Vec::with_capacity(count),
            columns: Vec::with_capacity(columns),
        })
    }

    /// Scores at a positive row and the row above it in one saved column.
    /// Within a block, decode once and obtain the neighbor from its delta.
    pub fn scores(&self, row: usize, col: usize) -> (Option<isize>, Option<isize>) {
        if col == 0 {
            return (Some(row as isize), Some(row as isize - 1));
        }
        let (first, range) = &self.columns[col - 1];
        let block = |b: usize| {
            (*first..*first + range.len())
                .contains(&b)
                .then(|| &self.blocks[range.start + b - first])
        };
        let b = (row - 1) / 64;
        let bit = (row - 1) % 64;
        let current = block(b);
        let value = current.map(|block| block.value(bit));
        let above = if row == 1 {
            Some(col as isize)
        } else if bit == 0 {
            // The previous block can be absent even when the current one is
            // retained, or present when the current one is absent.
            block(b - 1).map(|block| block.score)
        } else {
            current.map(|block| {
                value.unwrap() - ((block.p >> bit) & 1) as isize
                    + ((block.m >> bit) & 1) as isize
            })
        };
        (value, above)
    }

    /// On a retained traceback cell, a positive vertical delta means the cell
    /// above costs one less. Avoid decoding that score with two popcounts.
    pub fn can_insert(&self, row: usize, col: usize) -> bool {
        let (first, range) = &self.columns[col - 1];
        let local = (row - 1) / 64 - first;
        // At a block boundary the predecessor may have been pruned. Row zero
        // is the implicit matrix boundary and is always available.
        if row > 1 && row % 64 == 1 && local == 0 {
            return false;
        }
        self.blocks[range.start + local].p & (1 << ((row - 1) % 64)) != 0
    }

    fn push(&mut self, blocks: &[Block], first: usize, end: usize) {
        let start = self.blocks.len();
        self.blocks.extend_from_slice(&blocks[first..end]);
        self.columns.push((first, start..self.blocks.len()));
    }
}

/// The traced `global` for a one-block query, with the exact optimal distance
/// as its bound. With one block, any band shrink in `global` empties the band
/// and fails, which a known optimal path rules out. So every column retains
/// exactly block zero, and band arithmetic, strong reductions and slice copies
/// cannot change the trace. This records the identical blocks and column
/// ranges and returns the final score.
pub(crate) fn trace_single(profile: &Profile, target: Sequence<'_>, trace: &mut Trace) -> isize {
    debug_assert_eq!(profile.blocks, 1);
    let mut block = Block::new(64);
    for c in 0..target.len() {
        block.advance(profile.row(target.at(c))[0], Deltas::ANCHORED, 63);
        let start = trace.blocks.len();
        trace.blocks.push(block);
        trace.columns.push((0, start..start + 1));
    }
    block.value(profile.len - 1)
}

pub(crate) fn distance_global(profile: &Profile, target: Sequence<'_>, k: usize) -> Option<usize> {
    if profile.len.abs_diff(target.len()) > k {
        return None;
    }

    if profile.len <= 64 {
        let (mut pv, mut mv, mut score) = (u64::MAX, 0, profile.len as isize);
        let high = 1 << (profile.len - 1);
        for c in 0..target.len() {
            advance_single::<false>(
                profile.row(target.at(c))[0],
                high,
                &mut pv,
                &mut mv,
                &mut score,
            );
            // A remaining column can decrease the actual last row by at most 1.
            if score - (target.len() - c - 1) as isize > k as isize {
                return None;
            }
        }
        return (score <= k as isize).then_some(score as usize);
    }

    let column = global(profile, target, k, target.len(), None)?;
    column
        .value(profile.len)
        .filter(|&s| s <= k as isize)
        .map(|s| s as usize)
}

/// NW band restricted by both accumulated cost and mandatory remaining length
/// imbalance. `stop` may be an interior column for Hirschberg, but `target.len`
/// remains the full dimension: pruning must not mistake the midpoint for an end.
pub(crate) fn global(
    profile: &Profile,
    target: Sequence<'_>,
    bound: usize,
    stop: usize,
    mut trace: Option<&mut Trace>,
) -> Option<Column> {
    let (m, n, count) = (profile.len as isize, target.len() as isize, profile.blocks);
    if profile.len.abs_diff(target.len()) > bound {
        return None;
    }

    let mut k = bound.min(profile.len.max(target.len())) as isize;
    let padding = count as isize * 64 - m;
    let mut first = 0;
    let mut end = count.min(((k.min((k + m - n) / 2) + 1) as usize).div_ceil(64));
    let mut blocks: Vec<Block> = (0..count)
        .map(|b| Block::new((b as isize + 1) * 64))
        .collect();

    for c in 0..stop {
        let row = profile.row(target.at(c));
        let mut carry = Deltas::ANCHORED;
        // The last block's score stays in a register: the band checks below
        // reload a score from memory only after the band actually shrinks.
        let mut last_score = advance_words(&mut blocks[first..end], &row[first..end], &mut carry);
        let c = c as isize;
        k = k.min(
            last_score
                + (n - c - 1).max(m - end as isize * 64)
                + if end == count { padding } else { 0 },
        );

        // Whole-block bounds deliberately allow up to 63 extra rows. Every
        // 2048 columns, exact edge-cell bounds recover that slack.
        if end < count && end as isize * 64 - 1 <= k - last_score + 126 - n + c + m {
            let mut block = Block::new(last_score + 64 - carry.value());
            block.advance(row[end], carry, 63);
            last_score = block.score;
            blocks[end] = block;
            end += 1;
        }
        while end > first {
            if last_score < k + 64 && end as isize * 64 - 1 <= k - last_score + 127 - n + c + m {
                break;
            }
            end -= 1;
            if end > first {
                last_score = blocks[end - 1].score;
            }
        }
        while first < end {
            let score = blocks[first].score;
            if score < k + 64 && (first as isize + 1) * 64 - 1 >= score - k - n + m + c {
                break;
            }
            first += 1;
        }
        if c % 2048 == 0 {
            while end > first {
                let b = end - 1;
                let cells = if end == count {
                    64 - padding as usize
                } else {
                    64
                };
                let viable = (0..cells).rev().any(|i| {
                    let s = blocks[b].value(i);
                    s <= k && (b * 64 + i) as isize <= k - s - n + c + m + 1
                });
                if viable {
                    break;
                }
                end -= 1;
            }
            while first < end {
                let cells = if first + 1 == count {
                    64 - padding as usize
                } else {
                    64
                };
                let viable = (0..cells).rev().any(|i| {
                    let s = blocks[first].value(i);
                    s <= k && (first * 64 + i) as isize >= s - k - n + c + m
                });
                if viable {
                    break;
                }
                first += 1;
            }
        }
        if first == end {
            return None;
        }
        if let Some(saved) = trace.as_mut() {
            saved.push(&blocks, first, end);
        }
    }

    Some(Column {
        blocks,
        first,
        last: end - 1,
    })
}

struct SemiBest<'a> {
    k: isize,
    distance: Option<usize>,
    ends: Option<&'a mut Vec<isize>>,
}

impl SemiBest<'_> {
    /// Record an admissible score. Only a score-only search can finish at zero.
    #[inline]
    fn report(&mut self, score: isize, end: isize) -> bool {
        if score > self.k {
            return false;
        }
        self.distance = Some(score as usize);
        if let Some(ends) = &mut self.ends {
            if score < self.k {
                ends.clear();
            }
            ends.push(end);
            self.k = score;
            false
        } else if score == 0 {
            true
        } else {
            self.k = score - 1; // Ties cannot improve a score-only result.
            false
        }
    }
}

fn semi_single<const FREE: bool>(
    profile: &Profile,
    target: impl Iterator<Item = u8>,
    best: &mut SemiBest<'_>,
) {
    let (mut pv, mut mv, mut score) = (u64::MAX, 0, profile.len as isize);
    let high = 1 << (profile.len - 1);
    for (c, symbol) in target.enumerate() {
        if !FREE && c >= profile.len + best.k as usize {
            break;
        }
        advance_single::<FREE>(profile.row(symbol)[0], high, &mut pv, &mut mv, &mut score);
        if best.report(score, c as isize) {
            break;
        }
    }
}

pub(crate) fn semiglobal(
    profile: &Profile,
    target: impl Iterator<Item = u8>,
    bound: usize,
    free: bool,
    mut ends: Option<&mut Vec<isize>>,
    blocks: &mut Vec<Block>,
) -> Option<usize> {
    if let Some(ends) = &mut ends {
        ends.clear();
    }
    let mut best = SemiBest { k: bound.min(profile.len) as isize, distance: None, ends };
    // Edlib omits boundary zero for a whole-word query and a nonempty target.
    // That reporting convention does not exclude the empty substring as a score.
    if profile.len % 64 != 0 || best.ends.is_none() {
        best.report(profile.len as isize, -1);
    }
    if profile.len <= 64 {
        if free {
            semi_single::<true>(profile, target, &mut best);
        } else {
            semi_single::<false>(profile, target, &mut best);
        }
    } else {
        semi_blocks(profile, target, free, &mut best, blocks);
    }
    best.distance
}

fn semi_blocks(
    profile: &Profile,
    target: impl Iterator<Item = u8>,
    free: bool,
    best: &mut SemiBest<'_>,
    blocks: &mut Vec<Block>,
) {
    let count = profile.blocks;
    let full = profile.len / 64;
    let tail = profile.len % 64;
    let width = |b| if b < full { 64 } else { tail };
    let mut first = 0;
    let mut end = count.min((best.k as usize + 1).div_ceil(64));
    blocks.clear();
    blocks.reserve_exact(count);
    blocks.extend((1..=end).map(|b| Block::new((b * 64).min(profile.len) as isize)));
    for (c, symbol) in target.enumerate() {
        let k = best.k;
        let row = profile.row(symbol);
        let mut carry = Deltas { p: u64::from(!free), m: 0 };
        let full_end = end.min(full);
        let bottom = advance_words(&mut blocks[first..full_end], &row[first..full_end], &mut carry);
        if end > full {
            // Track the last real row: a zero is final on this column, not
            // several columns later when it reaches the bottom of the padding.
            blocks[full].advance(row[full], carry, tail - 1);
        }
        if end < count && bottom - carry.value() <= k && (row[end] & 1 != 0 || carry.m != 0) {
            let w = width(end);
            let mut block = Block::new(bottom - carry.value() + w as isize);
            block.advance(row[end], carry, w - 1);
            blocks.push(block);
            end += 1;
        } else {
            while end > first && blocks[end - 1].score >= k + 64 {
                end -= 1;
            }
        }
        if c % 2048 == 0 {
            while end > first && blocks[end - 1].all_larger(k, width(end - 1)) {
                end -= 1;
            }
        }
        if free {
            // A match can start at the next column even if this one missed.
            end = end.max(1);
        } else {
            while first < end && blocks[first].score >= k + 64 {
                first += 1;
            }
            if c % 2048 == 0 {
                while first < end && blocks[first].all_larger(k, width(first)) {
                    first += 1;
                }
            }
        }
        if end <= first {
            break;
        }
        blocks.truncate(end);
        if end == count && best.report(blocks[end - 1].score, c as isize) {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        EqualityPair, Profile, Sequence, Trace, distance_global, global, horizontal, semiglobal,
        trace_single,
    };

    #[test]
    fn score_only_infix_stops_on_the_first_exact_end() {
        let target = vec![b'a'; 1_000_000];
        for len in [1, 64, 65, 129] {
            let profile = Profile::new(Sequence::new(&target[..len]), &[]).unwrap();
            for bound in [0, usize::MAX] {
                let mut columns = 0;
                let counted = target.iter().copied().inspect(|_| columns += 1);
                assert_eq!(semiglobal(&profile, counted, bound, true, None, &mut Vec::new()), Some(0));
                assert_eq!(columns, len, "query length {len}, cutoff {bound}");
            }
        }
    }

    /// The scatter construction, kept only as the bit-plane test oracle.
    fn scatter_words(profile: &Profile, query: Sequence<'_>) -> Vec<u64> {
        let rows = profile.words.len() / profile.blocks;
        let mut words = vec![0; rows * profile.blocks];
        for i in 0..query.len() {
            words[profile.indices[query.at(i) as usize] as usize * profile.blocks + i / 64] |=
                1 << (i % 64);
        }
        words
    }

    #[test]
    fn bit_plane_profiles_equal_scatter_profiles() {
        let mut state = 0x9e37_79b9_u32;
        let mut next = |bound: u32| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            state % bound
        };
        let pairs = [EqualityPair::new(b'A', b'z'), EqualityPair::new(b'q', b'C')];
        // 1..=9 distinct symbols span plane counts 1, 2 and 3 and the scatter fallback.
        for symbols in 1..=9u32 {
            for len in [1, 2, 31, 63, 64, 65, 127, 128, 129, 200, 513] {
                let query: Vec<u8> = (0..len).map(|_| b"ACGTNacgt"[next(symbols) as usize]).collect();
                for reverse in [false, true] {
                    for extra in [&pairs[..0], &pairs[..]] {
                        let sequence = if reverse { Sequence::reversed(&query) } else { Sequence::new(&query) };
                        let profile = Profile::new(sequence, extra).unwrap();
                        let mut expected = scatter_words(&profile, sequence);
                        if !extra.is_empty() {
                            Profile::add_equalities(&mut expected, &profile.indices, profile.blocks, sequence, extra);
                        }
                        if len % 64 != 0 {
                            for row in expected.chunks_exact_mut(profile.blocks) {
                                row[profile.blocks - 1] |= u64::MAX << (len % 64);
                            }
                        }
                        assert_eq!(profile.words, expected, "{symbols} symbols, length {len}, reverse {reverse}");
                    }
                }
            }
        }
    }

    /// Myers' original formulation, kept only as the test oracle.
    fn reference(eq: u64, pv: u64, mv: u64) -> (u64, u64) {
        let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;
        (mv | !(xh | pv), pv & xh)
    }

    #[test]
    fn shortened_horizontal_deltas_equal_myers_expressions() {
        let edges = [0, 1, u64::MAX, u64::MAX - 1, 1 << 63, (1 << 63) - 1, 0x5555_5555_5555_5555];
        for &eq in &edges {
            for &pv in &edges {
                for &mv in &edges {
                    assert_eq!(horizontal(eq, pv, mv), reference(eq, pv, mv));
                }
            }
        }
        // SplitMix64: deterministic words, including carries across every bit.
        let mut state = 0x5eed_ed1b_4285_fee6_u64;
        let mut next = || {
            state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        };
        for _ in 0..1_000_000 {
            let (eq, pv) = (next(), next());
            // Valid Myers states never set a bit in both Pv and Mv.
            let mv = next() & !pv;
            assert_eq!(horizontal(eq, pv, mv), reference(eq, pv, mv));
        }
    }

    #[test]
    fn single_block_traces_equal_banded_global_traces() {
        let mut state = 0x7ace_5eed_u32;
        let mut next = |bound: u32| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            state % bound
        };
        let pairs = [EqualityPair::new(b'A', b'C')];
        // Exact optimal bounds, including distance 0 and length-imbalanced
        // targets, every one-block query length, padding and equalities.
        for symbols in [1, 2, 4, 9] {
            for len in 1..=64 {
                for extra in [&pairs[..0], &pairs[..]] {
                    for _ in 0..4 {
                        let query: Vec<u8> = (0..len)
                            .map(|_| b"ACGTNacgt"[next(symbols) as usize])
                            .collect();
                        let n = 1 + next(3 * len as u32) as usize;
                        let target: Vec<u8> = (0..n)
                            .map(|_| b"ACGTNacgt"[next(symbols) as usize])
                            .collect();
                        let profile = Profile::new(Sequence::new(&query), extra).unwrap();
                        let distance =
                            distance_global(&profile, Sequence::new(&target), len.max(n)).unwrap();
                        let mut expected = Trace::new(1, n).unwrap();
                        let column = global(
                            &profile,
                            Sequence::new(&target),
                            distance,
                            n,
                            Some(&mut expected),
                        )
                        .unwrap();
                        let mut actual = Trace::new(1, n).unwrap();
                        let score = trace_single(&profile, Sequence::new(&target), &mut actual);
                        assert_eq!(Some(score), column.value(len));
                        assert_eq!(score, distance as isize);
                        assert_eq!(actual.columns, expected.columns);
                        let words = |t: &Trace| {
                            t.blocks
                                .iter()
                                .map(|b| (b.p, b.m, b.score))
                                .collect::<Vec<_>>()
                        };
                        assert_eq!(words(&actual), words(&expected), "{query:?} {target:?}");
                    }
                }
            }
        }
    }
}
