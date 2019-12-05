extern crate denoise_engine;
extern crate openblas_src;
use denoise_engine::parse::{NoiseField, SwathElem};
use denoise_engine::estimate::*;
use std::io::prelude::*;
use std::fs;
use ndarray::{Array1, Array2, Array3, ArrayD};
use hdf5;
use hdf5::types::Array;

#[test]
fn test_k_est() {
    let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_arrs.hdf5", "r").unwrap();
    let mut x = handle.dataset("x").unwrap().read_2d().unwrap();
    let mut y = handle.dataset("y").unwrap().read_2d().unwrap();

    for a in x.iter_mut() {
        *a = *a/10000.0;
    }

    for a in y.iter_mut() {
        *a = *a/10000.0;
    }

    let w_:Array1<i64> = handle.dataset("w").unwrap().read_1d().unwrap(); //TODO: reshape into slice
    let mut w:Vec<usize> = Vec::new();
    for i in 0..5{
        w.push(w_[i] as usize);
    }
    let swath_bounds_:ArrayD<f64> =  handle.dataset("swath_bounds").unwrap().read_dyn().unwrap(); //TODO: reshape into nested slices
    let mut swath_bounds:Vec<&[SwathElem]> = Vec::new();

    let mut sb1:Vec<SwathElem> = Vec::new();
    let mut sb2:Vec<SwathElem> = Vec::new();
    let mut sb3:Vec<SwathElem>= Vec::new();
    let mut sb4:Vec<SwathElem> = Vec::new();
    let mut sb5:Vec<SwathElem> = Vec::new();

    for i in 0..4{
        sb1.push(SwathElem{ fa:swath_bounds_[[0,i,0]] as usize,
                                 la:swath_bounds_[[0,i,1]] as usize,
                                 fr:swath_bounds_[[0,i,2]] as usize,
                                 lr:swath_bounds_[[0,i,3]] as usize});
        
        sb2.push(SwathElem{ fa:swath_bounds_[[1,i,0]] as usize,
                                 la:swath_bounds_[[1,i,1]] as usize,
                                 fr:swath_bounds_[[1,i,2]] as usize,
                                 lr:swath_bounds_[[1,i,3]] as usize});

        sb3.push(SwathElem{ fa:swath_bounds_[[2,i,0]] as usize,
                                 la:swath_bounds_[[2,i,1]] as usize,
                                 fr:swath_bounds_[[2,i,2]] as usize,
                                 lr:swath_bounds_[[2,i,3]] as usize});
        
        sb4.push(SwathElem{ fa:swath_bounds_[[3,i,0]] as usize,
                                 la:swath_bounds_[[3,i,1]] as usize,
                                 fr:swath_bounds_[[3,i,2]] as usize,
                                 lr:swath_bounds_[[3,i,3]] as usize});
        
        
        sb5.push(SwathElem{ fa:swath_bounds_[[4,i,0]] as usize,
                                 la:swath_bounds_[[4,i,1]] as usize,
                                 fr:swath_bounds_[[4,i,2]] as usize,
                                 lr:swath_bounds_[[4,i,3]] as usize});

    }

    swath_bounds.push(&sb1);
    swath_bounds.push(&sb2);
    swath_bounds.push(&sb3);
    swath_bounds.push(&sb4);
    swath_bounds.push(&sb5);

    
    let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
    let lambda2_ = 1.0;
    let mu = 1.7899;
    let gamma = 2.0;
    
    let kvals = estimate_k_values(x.view(), y.view(),
                                  &w, &swath_bounds,
                                  mu, gamma, lambda_, lambda2_);
    println!("{:?}",kvals);
}
