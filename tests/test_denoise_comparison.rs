/* extern crate denoise_engine;
use denoise_engine::parse::{NoiseField, SwathElem};
use denoise_engine::apply::apply_swath_scale;
use denoise_engine::estimate::*;
use std::io::prelude::*;
use std::fs;
use ndarray::{Array1, Array2, Array3, ArrayD};
use hdf5;
use hdf5::types::Array;

#[test]
fn apply_comparison() {
    let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_arrs.hdf5", "r").unwrap();
    let mut x = handle.dataset("x").unwrap().read_2d().unwrap();
    let mut y = handle.dataset("y").unwrap().read_2d().unwrap();
    let mut k = handle.dataset("k").unwrap().read_1d().unwrap();    

    let swath_bounds_:ArrayD<f64> =  handle.dataset("swath_bounds").unwrap().read_dyn().unwrap(); 
    let mut swath_bounds:Vec<&[SwathElem]> = Vec::new();

    let mut sb1:Vec<SwathElem> = Vec::new();
    let mut sb2:Vec<SwathElem> = Vec::new();
    let mut sb3:Vec<SwathElem> = Vec::new();
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

    apply_swath_scale(x.view_mut(), y.view(), k.view(), &swath_bounds);

    
    let mut denoised:Array2<f32> = handle.dataset("denoised").unwrap().read_2d().unwrap();
    //println!("{}",x[(100,100)]);
    for i in 0..9992 {
        for j in 0..10400 {
            let a = denoised[(i,j)] as f64;
            let b = x[(i,j)];
            if (a-b).abs() >= 0.2 {
                println!("{} {} : {} {}", i,j,a,b);

                println!("around python\n{} {} {}\n {} {} {}", denoised[(i,j-1)], denoised[(i,j)], denoised[(i,j+1)], denoised[(i+1,j-1)], denoised[(i+1,j)], denoised[(i+1,j+1)]);
                println!("around rust\n{} {} {}\n {} {} {}", x[(i,j-1)], x[(i,j)], x[(i,j+1)], x[(i+1,j-1)], x[(i+1,j)], x[(i+1,j+1)]);

                assert!( (a-b).abs() < 0.2);
            }
        }
    }
}
*/
