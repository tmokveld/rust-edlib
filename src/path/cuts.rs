//! Retain forward bands at the columns Hirschberg would split. One reverse
//! sweep can certify cuts for many descendants without changing a single tie.

use super::{Output, TRACE_BYTES, WORD_BITS, append_alignment, split, split_rows, trace_fits, trace_fits_budget};
use crate::aligner::{AlignError, EqualityPair};
use crate::kernel::{Band, Block, Profile, Sequence, distance_global, global_recording};

struct Site {
    column: usize,
    state: State,
}

enum State {
    Unknown,
    Band { first: usize, blocks: Box<[Block]> },
    Cut { row: usize, score: usize },
}

pub(super) struct Cuts(Vec<Site>);

impl Cuts {
    pub fn new(m: usize, n: usize) -> Self {
        if m <= WORD_BITS || n == 1 { return Self(Vec::new()); }
        // Estimate the depth at which rectangles fit. Actual query splits can
        // be uneven: an absent site just leaves an ordinary Hirschberg split.
        let mut leaves = 1;
        let max_sites = TRACE_BYTES / (size_of::<Site>() + size_of::<Block>());
        while leaves <= n / 2 && 2 * leaves <= max_sites
            && !trace_fits(m.div_ceil(leaves).div_ceil(WORD_BITS), n.div_ceil(leaves))
        {
            leaves *= 2;
        }
        let mut sites = Vec::with_capacity(leaves - 1);
        fn plan(sites: &mut Vec<Site>, start: usize, len: usize, leaves: usize) {
            if leaves == 1 { return; }
            let half = len / 2;
            plan(sites, start, half, leaves / 2);
            sites.push(Site { column: start + half, state: State::Unknown });
            plan(sites, start + half, len - half, leaves / 2);
        }
        plan(&mut sites, 0, n, leaves);
        Self(sites)
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn score(&mut self, profile: &Profile, target: Sequence<'_>, k: usize) -> Option<usize> {
        if self.0.is_empty() {
            return distance_global(profile, target, k, &mut Vec::new());
        }
        let mut room = (TRACE_BYTES - self.metadata_bytes()) / size_of::<Block>();
        let mut i = 0;
        let mut next = self.0[0].column;
        let mut exhausted = false;
        let end = global_recording(profile, target, k, target.len(), |col, blocks, first, end| {
            if col != next { return; }
            if end - first > room {
                exhausted = true;
                next = usize::MAX;
                return; // No more checkpoints; keep scoring to the real endpoint.
            }
            self.0[i].state = State::Band { first, blocks: blocks[first..end].into() };
            room -= end - first;
            i += 1;
            next = self.0.get(i).map_or(usize::MAX, |site| site.column);
        });
        let d = end.and_then(|column| column.distance(profile.len, k));
        if exhausted {
            self.0 = Vec::new();
        } else if d.is_none() {
            for site in &mut self.0[..i] { site.state = State::Unknown; }
        }
        d
    }

    pub fn resolve(
        &mut self, query: &[u8], target: &[u8], pairs: &[EqualityPair], d: usize,
    ) -> Result<(), AlignError> {
        let Some(first) = self.0.first() else { return Ok(()) };
        let (m, n) = (query.len(), target.len());
        let stop = n - first.column;
        let profile = Profile::new(Sequence::reversed(query), pairs)?;
        let mut it = self.0.iter_mut().rev();
        let mut site = it.next();
        global_recording(&profile, Sequence::reversed(target), d, stop, |col, blocks, first, end| {
            let Some(s) = site.as_deref_mut() else { return };
            if n - col != s.column { return; }
            let State::Band { first: saved_first, blocks: saved } = &s.state else { unreachable!() };
            let left = Band { first: *saved_first, blocks: saved };
            let right = Band { first, blocks: &blocks[first..end] };
            s.state = {
                let mut rows = split_rows(m, n, s.column, d, left, right);
                let (row, score) = rows.next().expect("an optimum crosses every checkpoint");
                // At the root use the original tie rule. Elsewhere the ancestor
                // endpoints aren't known yet. Every local optimum extends to a
                // global optimum, so only a unique global crossing suffices.
                if s.column == n / 2 || rows.next().is_none() {
                    State::Cut { row, score }
                } else {
                    State::Unknown
                }
            };
            site = it.next();
        }).expect("the optimal reverse band must reach the first checkpoint");
        Ok(())
    }

    fn metadata_bytes(&self) -> usize {
        self.0.capacity() * size_of::<Site>()
    }

    pub fn append<O: Output>(
        &mut self, query: &[u8], target: &[u8], pairs: &[EqualityPair], d: usize,
        origin: (usize, usize, usize), output: &mut O,
    ) -> Result<(), AlignError> {
        let (m, n) = (query.len(), target.len());
        let blocks = m.div_ceil(WORD_BITS);
        let needs_trace = m != 0 && n != 0 && d != 0;
        if self.is_empty() || !needs_trace || n == 1 || trace_fits(blocks, n) {
            // Keep the original leaf/split boundary (it determines ties). If
            // its rectangle needs the remaining budget, relinquish unused cuts.
            if !self.is_empty() && needs_trace && !trace_fits_budget(blocks, n, TRACE_BYTES - self.metadata_bytes()) {
                self.0 = Vec::new();
            }
            return append_alignment(query, target, pairs, d, None, output);
        }
        let (q0, t0, d0) = origin;
        let middle = n / 2;
        let cut = self.0.binary_search_by_key(&(t0 + middle), |site| site.column).ok()
            .and_then(|i| match self.0[i].state {
                State::Cut { row, score } => Some((row - q0, score - d0)),
                _ => None,
            });
        let (row, left) = match cut {
            Some(cut) => cut,
            None => split(query, target, pairs, d, None)?,
        };
        self.append(&query[..row], &target[..middle], pairs, left, origin, output)?;
        self.append(&query[row..], &target[middle..], pairs, d - left,
            (q0 + row, t0 + middle, d0 + left), output)
    }
}
