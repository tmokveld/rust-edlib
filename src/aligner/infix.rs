use super::AlignError;

/// Literal infix search after the empty-input and length-gap checks.
/// Returns the start of a match, without a query profile or a list of ends.
pub(super) fn find_exact(query: &[u8], target: &[u8]) -> Result<Option<usize>, AlignError> {
    let m = query.len();
    debug_assert!(m > 0 && m <= target.len());
    if let [byte] = query {
        return Ok(target.iter().position(|b| b == byte));
    }
    let Some((head, tail)) = query.split_first_chunk::<8>() else {
        return Ok(target.windows(m).position(|window| window == query));
    };

    // Charge a full window for each failed verification. Once another one
    // would exceed a target-length budget, finish in linear time with KMP.
    let mut budget = target.len();
    for (start, window) in target.windows(m).enumerate() {
        if window.first_chunk::<8>() != Some(head) {
            continue;
        }
        if &window[8..] == tail {
            return Ok(Some(start));
        }
        budget -= m;
        if budget < m {
            // Only starts through `start` were ruled out, not overlapping windows.
            let next = start + 1;
            return Ok(kmp(query, &target[next..])?.map(|offset| next + offset));
        }
    }
    Ok(None)
}

fn kmp(query: &[u8], target: &[u8]) -> Result<Option<usize>, AlignError> {
    if target.len() < query.len() {
        return Ok(None);
    }
    // borders[n-1] is the longest proper prefix that is a suffix of query[..n].
    let advance = |mut n: usize, byte: u8, borders: &[usize]| {
        while n > 0 && query[n] != byte {
            n = borders[n - 1];
        }
        n + usize::from(query[n] == byte)
    };
    let mut borders = Vec::new();
    borders.try_reserve_exact(query.len()).map_err(|_| AlignError::InputTooLong)?;
    borders.push(0);
    let mut matched = 0;
    for &byte in &query[1..] {
        matched = advance(matched, byte, &borders);
        borders.push(matched);
    }

    matched = 0;
    for (i, &byte) in target.iter().enumerate() {
        matched = advance(matched, byte, &borders);
        if matched == query.len() {
            return Ok(Some(i + 1 - matched));
        }
    }
    Ok(None)
}
