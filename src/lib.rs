//! Native, safe Rust Levenshtein alignment using banded Myers bit vectors.
//!
//! Use [`distance`] for a global edit score or [`align`] for a configurable alignment.
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aligner;
mod kernel;
mod path;

pub use aligner::{
    AlignConfig, AlignError, Alignment, AlignmentMode, AlignmentTask, Cigar, CigarFormat,
    EditOperation, EqualityPair, align, distance,
};
