//! Native, safe Rust Levenshtein alignment using banded Myers bit vectors.
//!
//! Use [`distance`], [`prefix_distance`], or [`suffix_distance`] for a score
//! alone, or [`align`] for a configurable alignment.
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aligner;
mod kernel;
mod myers;
mod path;

pub use aligner::{
    AlignConfig, AlignError, Alignment, AlignmentMode, AlignmentTask, Cigar, CigarFormat,
    EditOperation, EqualityPair, align, distance, prefix_distance, suffix_distance,
};
