extern crate denoise_engine;
use denoise_engine::parse::{NoiseField, SwathElem};
use denoise_engine::read_from_archive::get_data_from_zip_path;

#[test]
fn test_load_sample(){
    let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    let zipval = get_data_from_zip_path(path);
    
    match zipval {
        Some((swath_bounds, w, noisefield, x16)) => {
            println!("Load successful");
        }

        None => {
            panic!("File parsed incorrectly");
        }
    }
}

/// Test the entire pipeline, comparing the final result to the compiled hdarr
#[test]
fn test_full_moty() {
}
    // let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_arrs.hdf5", "r").unwrap();
    // let mut x_gt = handle.dataset("x").unwrap().read_2d().unwrap();
    // let mut y_gt = handle.dataset("y").unwrap().read_2d().unwrap();
    // let mut k_gt = handle.dataset("k").unwrap().read_1d().unwrap();
