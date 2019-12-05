
/*extern crate denoise_engine;
use denoise_engine::parse::NoiseField;
use std::io::prelude::*;
use std::fs;
use ndarray::Array2;
use hdf5;
use hdf5::types::Array;

#[test]
fn test_noise_field() {
    let xmlfile = fs::read_to_string("/mnt/D2/Data/Sentinel/test/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.SAFE/annotation/calibration/noise-s1a-ew-grd-hv-20180902t164932-20180902t165032-023522-028faa-002.xml").unwrap();
    let m = NoiseField::new(&xmlfile);
    //panic!("");adf
    //let ot = hdf5::File::open("tests/outarr.hdf5", "w").unwrap();
    //ods.write(&m.data).unwrap();
    //let ods = ot.new_dataset::<f64>().create("az", (9992, 10400)).unwrap();
    
    let handle = hdf5::File::open("tests/testarr.hdf5","r").unwrap();
    let gt = handle.dataset("full").unwrap().read_2d().unwrap();
    println!("{} {}", m.data[(0,0)], gt[(0,0)]);
    assert!(m.data.all_close(&gt, 1e-4));
}
*/
