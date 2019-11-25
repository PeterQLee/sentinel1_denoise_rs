extern crate denoise_engine;
use denoise_engine::parse::NoiseField;
use std::io::prelude::*;
use std::fs;
use std::env;
use read::ZipArchive;

fn main() {
    if env::args().len() < 3 {
        println!("./noise_to_tiff [path_to_zip_file] [path_to_output_file]");
    }

    let args:Vec<_> = env::args().collect();
    let zip_path = args[0];
    let tiff_path = args[1];

    let handle = File::open(
    let zipfile = zip::ZipArchive::new(reader).unwrap();
}
