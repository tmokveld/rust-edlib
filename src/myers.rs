//! Score-only anchored Myers search.

use crate::AlignError;
use crate::kernel::{Block, Deltas, advance_words, initial_slack};

/// The caller has trimmed the shared pinned end and cropped the free end.
/// `lower` is a proven lower bound; `limit` is the inclusive cutoff.
pub(crate) fn distance<const REVERSE: bool>(
    query: &[u8],
    target: &[u8],
    lower: usize,
    limit: usize,
) -> Result<Option<usize>, AlignError> {
    let m = query.len();
    debug_assert!(lower <= limit && limit <= m);
    if lower == m {
        return Ok(Some(m));
    }
    if m <= 64 {
        return Ok(single::<REVERSE>(query, target, lower, limit));
    }

    // Even an inactive edge reconstructed via deletions costs at most row +
    // column. Leave room for a whole word in the signed pruning expressions.
    m.checked_add(target.len())
        .and_then(|sum| sum.checked_add(64))
        .filter(|&sum| sum <= isize::MAX as usize)
        .ok_or(AlignError::InputTooLong)?;
    let profile = Profile::new::<REVERSE>(query)?;
    let mut blocks = Vec::new();
    blocks
        .try_reserve_exact(profile.count)
        .map_err(|_| AlignError::InputTooLong)?;
    // Only use the empty flank as an incumbent when it is within the cutoff.
    let empty = (limit == m).then_some(m);
    let upper = limit.min(m - 1);
    let gap = m.saturating_sub(target.len());
    let mut bound = if limit < m {
        limit
    } else {
        (gap + initial_slack(gap, target.len())).max(lower).min(upper)
    };
    loop {
        let reach = match banded::<REVERSE>(&profile, target, lower, bound, &mut blocks) {
            Pass::Hit(score) => return Ok(Some(score)),
            Pass::Miss(_) if bound == upper => return Ok(empty),
            Pass::Miss(reach) => reach,
        };
        // The failed pass's query extent predicts a useful budget, not a
        // lower bound. Start again from the pinned end at the new budget.
        let forecast = bound.saturating_mul(m) / reach;
        let grown = bound + (bound - gap);
        bound = grown.max(forecast.min(bound.saturating_mul(16))).min(upper);
    }
}

#[inline]
fn bytes<const REVERSE: bool>(slice: &[u8]) -> impl Iterator<Item = u8> + '_ {
    let mut iter = slice.iter().copied();
    std::iter::from_fn(move || if REVERSE { iter.next_back() } else { iter.next() })
}

struct Profile {
    len: usize,
    count: usize,
    // Only query symbols get their own row; unknown target bytes share row zero.
    indices: [u16; 256],
    words: Vec<u64>,
}

impl Profile {
    fn new<const REVERSE: bool>(query: &[u8]) -> Result<Self, AlignError> {
        let count = query.len().div_ceil(64);
        let mut indices = [0; 256];
        let mut rows = 1;
        for &symbol in query {
            let id = &mut indices[symbol as usize];
            if *id == 0 {
                *id = rows as u16;
                rows += 1;
            }
        }
        let size = count.checked_mul(rows).ok_or(AlignError::InputTooLong)?;
        let mut words = Vec::new();
        words
            .try_reserve_exact(size)
            .map_err(|_| AlignError::InputTooLong)?;
        words.resize(size, 0);
        // For small alphabets, collect bit planes in registers rather than
        // repeatedly read-modify-writing a symbol's equality word in memory.
        match rows {
            2 => Self::planes::<REVERSE, 1>(&mut words, &indices, count, rows, query),
            3..=4 => Self::planes::<REVERSE, 2>(&mut words, &indices, count, rows, query),
            5..=8 => Self::planes::<REVERSE, 3>(&mut words, &indices, count, rows, query),
            _ => {
                for (i, symbol) in bytes::<REVERSE>(query).enumerate() {
                    let row = indices[symbol as usize] as usize;
                    words[row * count + i / 64] |= 1 << (i % 64);
                }
            }
        }
        Ok(Self {
            len: query.len(),
            count,
            indices,
            words,
        })
    }

    fn planes<const REVERSE: bool, const P: usize>(
        words: &mut [u64],
        indices: &[u16; 256],
        count: usize,
        rows: usize,
        query: &[u8],
    ) {
        let mut query = bytes::<REVERSE>(query);
        for block in 0..count {
            let mut planes = [0u64; P];
            for (i, symbol) in query.by_ref().take(64).enumerate() {
                let id = u64::from(indices[symbol as usize]);
                for (p, plane) in planes.iter_mut().enumerate() {
                    *plane |= ((id >> p) & 1) << i;
                }
            }
            // Unused bits match no nonzero ID, so the final word needs no mask.
            for row in 1..rows {
                let mut mask = u64::MAX;
                for (p, &plane) in planes.iter().enumerate() {
                    mask &= if (row >> p) & 1 != 0 { plane } else { !plane };
                }
                words[row * count + block] = mask;
            }
        }
    }

    #[inline]
    fn row(&self, symbol: u8) -> &[u64] {
        let start = self.indices[symbol as usize] as usize * self.count;
        &self.words[start..start + self.count]
    }
}

fn single<const REVERSE: bool>(
    query: &[u8], target: &[u8], lower: usize, limit: usize,
) -> Option<usize> {
    let mut row = [0u64; 256];
    for (i, symbol) in bytes::<REVERSE>(query).enumerate() {
        row[symbol as usize] |= 1 << i;
    }
    let m = query.len();
    let mut best = m.min(limit + 1); // No hit yet if the empty flank exceeds the cutoff.
    let mut block = Block::new(m as isize);
    for (c, symbol) in bytes::<REVERSE>(target).enumerate() {
        if c + 1 >= m + best {
            break; // The length gap alone rules out a strict improvement.
        }
        block.advance(row[symbol as usize], Deltas::ANCHORED, m - 1);
        best = best.min(block.score as usize);
        // Each remaining column can reduce the last-row score by at most one.
        if best == lower || block.score as usize >= best + (target.len() - c - 1) {
            break;
        }
    }
    (best <= limit).then_some(best)
}

enum Pass {
    Hit(usize),
    // Query extent of the last nonempty band, not a reachable endpoint.
    Miss(usize),
}

fn banded<const REVERSE: bool>(
    profile: &Profile,
    target: &[u8],
    lower: usize,
    bound: usize,
    blocks: &mut Vec<Block>,
) -> Pass {
    let m = profile.len;
    let count = profile.count;
    let last = count - 1;
    // Track the last real row's score directly, rather than padding with
    // wildcards and delaying/recovering its scores from later columns.
    let tail = m - last * 64;
    let width = |b| if b == last { tail } else { 64 };
    let mut first = 0;
    let mut k = bound as isize;
    let mut found = false;
    blocks.clear();
    blocks.extend(
        (1..=count.min((bound + 1).div_ceil(64))).map(|b| Block::new((b * 64).min(m) as isize)),
    );
    for (c, symbol) in bytes::<REVERSE>(target).enumerate() {
        if c >= m + k as usize {
            break; // Prefix length c + 1 exceeds m + k.
        }
        let row = profile.row(symbol);
        let mut end = blocks.len();
        let full_end = end.min(last);
        let mut carry = Deltas::ANCHORED;
        let last_score =
            advance_words(&mut blocks[first..full_end], &row[first..full_end], &mut carry);
        if end == count {
            blocks[last].advance(row[last], carry, tail - 1);
        }
        if end < count
            && last_score - carry.value() <= k
            && (row[end] & 1 != 0 || carry.m != 0)
        {
            // The inactive next word starts as deletions from the previous
            // column's bottom cell. A match or -1 carry can open that word.
            let previous = last_score - carry.value();
            let w = width(end);
            let mut block = Block::new(previous + w as isize);
            block.advance(row[end], carry, w - 1);
            blocks.push(block);
            end += 1;
        } else {
            while end > first && blocks[end - 1].score >= k + 64 {
                end -= 1;
            }
        }
        while first < end && blocks[first].score >= k + 64 {
            first += 1;
        }
        // Cheap whole-word pruning allows slack; periodically inspect the
        // vertical deltas to drop an edge word with no viable real cells.
        if c % 2048 == 0 {
            while end > first && blocks[end - 1].all_larger(k, width(end - 1)) {
                end -= 1;
            }
            while first < end && blocks[first].all_larger(k, width(first)) {
                first += 1;
            }
        }
        if first == end {
            break;
        }
        blocks.truncate(end);
        if end == count && blocks[last].score <= k {
            k = blocks[last].score;
            found = true;
            if k as usize == lower {
                break;
            }
        }
    }
    if found {
        Pass::Hit(k as usize)
    } else {
        Pass::Miss((blocks.len() * 64).min(m))
    }
}
