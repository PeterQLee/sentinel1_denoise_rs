extern crate denoise_engine;
use denoise_engine::parse::NoiseField;
use std::io::prelude::*;
use std::fs;
fn main() {
    let xmlfile = fs::read_to_string("/mnt/D2/Data/Sentinel/test/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.SAFE/annotation/calibration/noise-s1a-ew-grd-hv-20180902t164932-20180902t165032-023522-028faa-002.xml").unwrap();
    let m = NoiseField::new(&xmlfile);
    panic!("");
}
