#![feature(get_mut_unchecked)]
#[macro_use(s)]
extern crate ndarray;
//extern crate ndarray_linalg;
extern crate lapack;
extern crate blas;
extern crate openblas_src;
extern crate libc;
pub mod parse;

pub mod apply;
pub mod read_from_archive;
pub mod py_interface;
#[macro_use]
pub mod prep_lp;
pub mod est_lp;
pub mod estimate;
