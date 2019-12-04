extern crate denoise_engine;
use denoise_engine::parse::NoiseField;
use denoise_engine::estimate::*;
use std::io::prelude::*;
use std::fs;
use ndarray::Array2;
use hdf5;
use hdf5::types::Array;

#[test]
fn test_k_est() {
    let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_vol.hdf5", "r").unwrap();
    let x = handle.dataset("x").unwrap().read_2d().unwrap();
    let y = handle.dataset("y").unwrap().read_2d().unwrap();


    let w_ = handle.dataset("w").unwrap().read_1d().unwrap(); //TODO: reshape into slice
    let w = w_.into_vec();
    let swath_bounds_ =  handle.dataset("swath_bounds").unwrap().read_1d().unwrap(); //TODO: reshape into nested slices
    let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
    let lambda2_ = 1.0;
    let mu = 1.7899;
    let gamma = 2.0;
    
    let kvals = estimate_k_values(x.as_view(), y.as_view().
                                  &w, swath_bounds,
                                  mu, gamma, lambda_, lambda2_);
    
}
