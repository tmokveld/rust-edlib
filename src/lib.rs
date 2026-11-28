//! Native, safe Rust Levenshtein alignment using banded Myers bit vectors.
//!
//! Use [`distance`], [`prefix_distance`], [`suffix_distance`], or [`infix_distance`]
//! for a score alone. Their bounded variants, [`distance_at_most`],
//! [`prefix_distance_at_most`], [`suffix_distance_at_most`], and
//! [`infix_distance_at_most`], return `Ok(Some(d))` when the minimum score `d`
//! is at most the cutoff, or `Ok(None)` otherwise. An infix score-only search
//! stops at the first exact match; [`align`] collects optimal ends even for
//! [`AlignmentTask::Distance`].
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aligner;
mod kernel;
mod myers;
mod path;
mod word;

pub use aligner::{
    AlignConfig, AlignError, Alignment, AlignmentMode, AlignmentTask, Cigar, CigarFormat,
    EditOperation, EqualityPair, GlobalProbe, Run, RunAlignment, align, align_runs, distance,
    distance_at_most, distance_with_hint, infix_distance, infix_distance_at_most, prefix_distance,
    prefix_distance_at_most, suffix_distance, suffix_distance_at_most,
};
