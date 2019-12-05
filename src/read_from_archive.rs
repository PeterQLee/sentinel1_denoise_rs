/// Module for reading from sentinel-1 zip archive.
use crate::parse::NoiseField;
use zip::read::{ZipArchive};
use std::io::prelude::*;
use std::fs::File;

/// Get the dateid prefix for the relevant files.
///ex: S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.SAFE-report-20180902T190505.pdf
/// Returns the uppercase id and the lower case id
fn get_id_prefix(token:&str) -> (String, String){
    // Takes the pdf name as input
    
}

/// Gets the original 16-bit image and noise field from a given path.
pub fn get_data_from_zip_path(path:&str) -> () {
    let file_h = fs::File::open(path);
    match file_h {
        Ok(f) => {
            let ziparch = zip:ZipArchive(f);

            // TODO: find ids and stuff
            
        }
        None => {
            println!("Cannot open zipfile {}", path);
        }
    }

}
