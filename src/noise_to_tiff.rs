extern crate denoise_engine;
use denoise_engine::parse::NoiseField;
use std::io::prelude::*;
use std::fs::File;
use std::env;
use zip::read::{ZipArchive};


fn get_xml_and_tiff_buffers(zipfile:&mut ZipArchive<File>) {
    for i in 0..zipfile.len() {
        let entry = zipfile.by_index(i).unwrap();
        let name = entry.name();
        // Use name to determine if this entry is a feasible one...
    }
}

fn main() {
    if env::args().len() < 3 {
        println!("./noise_to_tiff [path_to_zip_file] [path_to_output_file]");
    }

    let args:Vec<_> = env::args().collect();
    let zip_path = &args[0];
    let tiff_path = &args[1];

    let handle = File::open(zip_path).unwrap();
    let zipfile = ZipArchive::new(handle).unwrap();



    
}
