extern crate denoise_engine;
use denoise_engine::parse::{NoiseField, SwathElem};
use std::fs;
use ndarray::{Array1, Array2, Array3, ArrayD};
use hdf5;
use hdf5::types::Array;
//use itertools::Itertools;



#[test]
fn subswath_bounds() {
    let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_arrs.hdf5", "r").unwrap();
    let xmlfile = fs::read_to_string("/mnt/D2/Data/Sentinel/test/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.SAFE/annotation/s1a-ew-grd-hv-20180902t164932-20180902t165032-023522-028faa-002.xml").unwrap();

    let (swath_bounds, w) = SwathElem::new(&xmlfile);
    let w_:Array1<i64> = handle.dataset("w").unwrap().read_1d().unwrap(); //TODO: reshape into slice

    let swath_bounds_:ArrayD<f64> = handle.dataset("swath_bounds").unwrap().read_dyn().unwrap(); //TODO: reshape into nested slices

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

    for (x,y) in sb1.iter().zip(swath_bounds[0].iter()) {
        assert!((x.fa == y.fa && x.fr == y.fr && x.la == y.la && x.lr == y.lr));
    }

    for (x,y) in sb2.iter().zip(swath_bounds[1].iter()) {
        assert!((x.fa == y.fa && x.fr == y.fr && x.la == y.la && x.lr == y.lr));
    }

    for (x,y) in sb3.iter().zip(swath_bounds[2].iter()) {
        assert!((x.fa == y.fa && x.fr == y.fr && x.la == y.la && x.lr == y.lr));
    }
    for (x,y) in sb4.iter().zip(swath_bounds[3].iter()) {
        assert!((x.fa == y.fa && x.fr == y.fr && x.la == y.la && x.lr == y.lr));
    }

    for (x,y) in sb5.iter().zip(swath_bounds[4].iter()) {
        assert!((x.fa == y.fa && x.fr == y.fr && x.la == y.la && x.lr == y.lr));
    }

    for (x,y) in w.iter().zip(w_.iter()) {
        assert!((*x as usize) == (*y as usize));
    }
}

