#![allow(dead_code)]

use rust_edlib::aligner::{
    AlignConfig, Alignment, AlignmentMode, AlignmentTask, CigarFormat, EqualityPair,
};
use std::ops::Range;

pub const SEED: u64 = 0x5eed_ed1b_4285_fee6;
pub const MODES: [AlignmentMode; 3] = [
    AlignmentMode::Global,
    AlignmentMode::Prefix,
    AlignmentMode::Infix,
];
pub const TASKS: [AlignmentTask; 3] = [
    AlignmentTask::Distance,
    AlignmentTask::Locations,
    AlignmentTask::Path,
];
pub const BOUNDARY_LENGTHS: [usize; 12] = [0, 1, 2, 63, 64, 65, 127, 128, 129, 191, 256, 513];

#[derive(Clone, Debug)]
pub struct Expected {
    pub distance: usize,
    pub ends: Vec<isize>,
    pub starts: Vec<usize>,
}

pub fn equal(a: u8, b: u8, pairs: &[EqualityPair]) -> bool {
    a == b
        || pairs
            .iter()
            .any(|p| (a == p.first() && b == p.second()) || (a == p.second() && b == p.first()))
}

/// Independent scalar Wagner-Fischer DP. Each cell carries its earliest optimal
/// start, so locations are checked without trusting a bit-vector traceback.
/// Space is O(target length); time is O(query length * target length).
pub fn quadratic(
    query: &[u8],
    target: &[u8],
    mode: AlignmentMode,
    pairs: &[EqualityPair],
) -> Expected {
    let infix = mode == AlignmentMode::Infix;
    let mut row: Vec<(usize, usize)> = (0..=target.len())
        .map(|j| if infix { (0, j) } else { (j, 0) })
        .collect();
    for (i, &q) in query.iter().enumerate() {
        let mut diagonal = row[0];
        row[0] = (i + 1, 0);
        for (j, &t) in target.iter().enumerate() {
            let above = row[j + 1];
            row[j + 1] = (diagonal.0 + usize::from(!equal(q, t, pairs)), diagonal.1)
                .min((above.0 + 1, above.1))
                .min((row[j].0 + 1, row[j].1));
            diagonal = above;
        }
    }
    let distance = if mode == AlignmentMode::Global {
        row[target.len()].0
    } else {
        row.iter().map(|cell| cell.0).min().unwrap()
    };
    let offsets: Vec<usize> = if mode == AlignmentMode::Global {
        vec![target.len()]
    } else {
        row.iter()
            .enumerate()
            .filter_map(|(j, cell)| {
                // Native edlib only observes column zero through the padding in a
                // partial word. This documented reporting convention is separate
                // from the mathematical DP optimum, which always includes it.
                let omitted =
                    j == 0 && !query.is_empty() && query.len() % 64 == 0 && !target.is_empty();
                (cell.0 == distance && !omitted).then_some(j)
            })
            .collect()
    };
    Expected {
        distance,
        ends: offsets.iter().map(|&j| j as isize - 1).collect(),
        starts: offsets.iter().map(|&j| row[j].1).collect(),
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Snapshot {
    pub distance: Option<usize>,
    pub ends: Vec<isize>,
    pub starts: Option<Vec<usize>>,
    pub path: Option<Vec<u8>>,
    pub standard: Option<String>,
    pub extended: Option<String>,
}

impl Snapshot {
    pub fn capture(result: &Alignment) -> Self {
        let ends = result.end_locations().to_vec();
        let offsets: Vec<usize> = ends.iter().map(|&end| (end + 1) as usize).collect();
        assert_eq!(result.end_offsets().len(), offsets.len());
        assert_eq!(result.end_offsets().collect::<Vec<_>>(), offsets);
        let starts = result.start_locations().map(<[_]>::to_vec);
        let spans = starts.as_ref().map(|starts| {
            starts
                .iter()
                .zip(&offsets)
                .map(|(&start, &end)| start..end)
                .collect::<Vec<_>>()
        });
        assert_eq!(
            result.locations().map(|spans| spans.len()),
            spans.as_ref().map(Vec::len)
        );
        assert_eq!(result.locations().map(Iterator::collect::<Vec<_>>), spans);
        Self {
            distance: result.edit_distance(),
            ends,
            starts,
            path: result
                .path()
                .map(|path| path.iter().map(|&op| op as u8).collect()),
            standard: result.cigar(CigarFormat::Standard).map(|c| c.to_string()),
            extended: result.cigar(CigarFormat::Extended).map(|c| c.to_string()),
        }
    }
}

pub fn limits(distance: usize) -> Vec<Option<usize>> {
    let mut limits = vec![
        None,
        Some(0),
        Some(distance.saturating_sub(1)),
        Some(distance),
        Some(distance.saturating_add(1)),
        Some(usize::MAX),
    ];
    limits.dedup();
    limits
}

pub fn assert_snapshot(
    query: &[u8],
    target: &[u8],
    config: &AlignConfig<'_>,
    expected: &Expected,
    result: &Snapshot,
) {
    let context = || {
        format!(
            "seed={SEED:#x} mode={:?} task={:?} limit={:?} pairs={:?} query={query:?} target={target:?}",
            config.mode,
            config.task,
            config.max_distance,
            config
                .additional_equalities
                .iter()
                .map(|p| (p.first(), p.second()))
                .collect::<Vec<_>>()
        )
    };
    let found = config.max_distance.is_none_or(|k| expected.distance <= k);
    assert_eq!(
        result.distance,
        found.then_some(expected.distance),
        "{}",
        context()
    );
    if !found {
        assert!(result.ends.is_empty(), "{} {result:?}", context());
        assert!(result.starts.is_none(), "{} {result:?}", context());
        assert!(result.path.is_none(), "{} {result:?}", context());
        assert!(
            result.standard.is_none() && result.extended.is_none(),
            "{}",
            context()
        );
        return;
    }
    assert_eq!(result.ends, expected.ends, "{}", context());
    if config.task == AlignmentTask::Distance {
        assert!(result.starts.is_none(), "{}", context());
    } else {
        assert_eq!(
            result.starts.as_deref(),
            Some(expected.starts.as_slice()),
            "{}",
            context()
        );
    }
    if config.task == AlignmentTask::Path {
        let path = result
            .path
            .as_ref()
            .unwrap_or_else(|| panic!("missing computed path: {}", context()));
        let span = expected.starts[0]..(expected.ends[0] + 1) as usize;
        replay(
            query,
            target,
            config.additional_equalities,
            span,
            expected.distance,
            path,
        )
        .unwrap_or_else(|error| panic!("{error}: {} path={path:?}", context()));
        check_cigar(result.standard.as_deref().unwrap(), path, false)
            .unwrap_or_else(|error| panic!("{error}: {}", context()));
        check_cigar(result.extended.as_deref().unwrap(), path, true)
            .unwrap_or_else(|error| panic!("{error}: {}", context()));
    } else {
        assert!(result.path.is_none(), "{}", context());
        assert!(
            result.standard.is_none() && result.extended.is_none(),
            "{}",
            context()
        );
    }
}

/// Accept numeric operations to replay both implementations independently.
pub fn replay(
    query: &[u8],
    target: &[u8],
    pairs: &[EqualityPair],
    span: Range<usize>,
    distance: usize,
    path: &[u8],
) -> Result<(), String> {
    if span.start > span.end || span.end > target.len() {
        return Err(format!("invalid span {span:?}"));
    }
    let (mut q, mut t, mut cost) = (0, span.start, 0);
    for (index, &op) in path.iter().enumerate() {
        match op {
            0 | 3 => {
                let (&qb, &tb) = query
                    .get(q)
                    .zip(target.get(t).filter(|_| t < span.end))
                    .ok_or_else(|| format!("diagonal out of bounds at operation {index}"))?;
                if equal(qb, tb, pairs) != (op == 0) {
                    return Err(format!(
                        "wrong match/mismatch at operation {index}: {qb} versus {tb}"
                    ));
                }
                q += 1;
                t += 1;
                cost += usize::from(op == 3);
            }
            1 if q < query.len() => {
                q += 1;
                cost += 1;
            }
            2 if t < span.end => {
                t += 1;
                cost += 1;
            }
            _ => return Err(format!("invalid operation {op} at index {index}")),
        }
    }
    if (q, t, cost) != (query.len(), span.end, distance) {
        return Err(format!(
            "path consumed/cost {:?}; expected {:?}",
            (q, t, cost),
            (query.len(), span.end, distance)
        ));
    }
    Ok(())
}

/// Parse CIGAR rather than reproducing its formatter. Check every expanded
/// operation and forbid zero-length or adjacent identical runs.
pub fn check_cigar(cigar: &str, path: &[u8], extended: bool) -> Result<(), String> {
    let (mut count, mut index, mut previous, mut digits) = (0usize, 0usize, None, false);
    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            if !digits && byte == b'0' {
                return Err("CIGAR has a leading zero".into());
            }
            count = count
                .checked_mul(10)
                .and_then(|n| n.checked_add((byte - b'0') as usize))
                .ok_or("CIGAR count overflow")?;
            digits = true;
            continue;
        }
        if !digits
            || count == 0
            || previous == Some(byte)
            || count > path.len().saturating_sub(index)
        {
            return Err(format!("invalid CIGAR run in {cigar:?}"));
        }
        for &op in &path[index..index + count] {
            let valid = match (byte, extended, op) {
                (b'I', _, 1) | (b'D', _, 2) | (b'=', true, 0) | (b'X', true, 3) => true,
                (b'M', false, 0 | 3) => true,
                _ => false,
            };
            if !valid {
                return Err(format!("CIGAR {cigar:?} disagrees with operation {index}"));
            }
        }
        index += count;
        previous = Some(byte);
        count = 0;
        digits = false;
    }
    if digits || index != path.len() {
        return Err(format!("incomplete CIGAR {cigar:?}"));
    }
    Ok(())
}

#[derive(Clone)]
pub struct Rng(u64);
impl Rng {
    pub fn new(seed: u64) -> Self {
        assert_ne!(seed, 0);
        Self(seed)
    }
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    pub fn bytes(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.next() as u8).collect()
    }
}

pub fn binary_strings(max_len: usize) -> Vec<Vec<u8>> {
    (0..=max_len)
        .flat_map(|len| {
            (0..1usize << len).map(move |bits| (0..len).map(|i| ((bits >> i) & 1) as u8).collect())
        })
        .collect()
}

/// Mutations include every edit orientation; inputs and failures are reproducible
/// without depending on a platform RNG or UTF-8 conversion.
pub fn mutated(rng: &mut Rng, query: &[u8]) -> Vec<u8> {
    let mut target = Vec::new();
    for &byte in query {
        match rng.next() % 13 {
            0 => {}
            1 => target.extend([byte, rng.next() as u8]),
            2 => target.push(byte.wrapping_add(1)),
            _ => target.push(byte),
        }
    }
    target
}

/// Closed-form fixture: each B requires an edit against an all-A query. The
/// flanking A bytes make the full target the unique optimal end and start.
/// This avoids quadratic reference work in strong-band/Hirschberg stress cases.
pub fn known_substitutions(len: usize, differences: usize) -> (Vec<u8>, Vec<u8>, Expected) {
    assert!(differences + 2 <= len);
    let query = vec![b'A'; len];
    let mut target = query.clone();
    let begin = (len - differences) / 2;
    target[begin..begin + differences].fill(b'B');
    (
        query,
        target,
        Expected {
            distance: differences,
            ends: vec![len as isize - 1],
            starts: vec![0],
        },
    )
}

/// A supersequence realizes the length-difference lower bound exactly, giving
/// an independent large global-indel oracle in either query/target direction.
pub fn known_insertions(len: usize, stride: usize) -> (Vec<u8>, Vec<u8>, Expected) {
    assert!(stride > 0);
    let mut rng = Rng::new(SEED ^ 0x1ade1);
    let query: Vec<u8> = rng.bytes(len).into_iter().map(|byte| byte % 251).collect();
    let mut target = Vec::with_capacity(len + len / stride + 1);
    for (index, &byte) in query.iter().enumerate() {
        target.push(byte);
        if index % stride == stride / 2 {
            target.push(255);
        }
    }
    let expected = Expected {
        distance: target.len() - query.len(),
        ends: vec![target.len() as isize - 1],
        starts: vec![0],
    };
    (query, target, expected)
}
