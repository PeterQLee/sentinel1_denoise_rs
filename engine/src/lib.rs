#![feature(get_mut_unchecked)]
//! Backbone engine for Sentinel-1 noise floor removal.
//! Main interfaces are within the [interface] and [postprocess] modules

extern crate ndarray;

extern crate lapack;
extern crate blas;
extern crate openblas_src;
extern crate libc;
pub mod parse;

pub mod apply;
pub mod read_from_archive;
pub mod interface;
#[macro_use]
pub mod prep_lp;
pub mod est_lp;
pub mod estimate;
pub mod postprocess;
