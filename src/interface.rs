
use crate::parse::{NoiseField, SwathElem};
use crate::read_from_archive::get_data_from_zip_path;
use crate::apply::{apply_swath_scale, prep_measurement, restore_scale};
use crate::estimate::*;
extern crate libc;
use ndarray::{Array1, Array2};
//use std::mem;
//use libc;
use std::ptr;

#[repr(C)] pub struct OutArr {data:*mut libc::c_double, rows:libc::c_int, cols:libc::c_int}

///Entry point for denoising zip file from C api
#[no_mangle]
pub extern fn denoise_zip(path:*const u8, pathlen:libc::c_int) -> OutArr{
    
    let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
    let lambda2_ = 1.0;
    let mu = 1.7899;
    let gamma = 2.0;
    // parse the path arguments
    let zip_path = unsafe {std::slice::from_raw_parts(path, (pathlen as usize)*std::mem::size_of::<u8>())};
    let zipval = get_data_from_zip_path(std::str::from_utf8(zip_path).unwrap());

    
    match zipval {
        Some((swath_bounds, w, mut noisefield, x16)) => {
            let mut x_:Option<Array2<f64>> = None;
            {
                let mut y = noisefield.data.view_mut();
                x_ = Some(prep_measurement(x16.view(), y));
            }
            let mut x = x_.unwrap();

            let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
            
            let k = estimate_k_values(x.view(), noisefield.data.view(), &w, &sb, mu, gamma, lambda_, lambda2_);
            
            {
                let mut y = noisefield.data.view_mut();
                restore_scale(x.view_mut(), y);
            }
            apply_swath_scale(x.view_mut(), noisefield.data.view(), k.view(), &sb);

            let (rows, cols) = x.dim();
            assert!(x.is_standard_layout());
            let result = OutArr {
                data:x.as_mut_ptr(),
                rows:x.shape()[0] as libc::c_int,
                cols:x.shape()[1] as libc::c_int
            };

            std::mem::forget(x);

            return result;
        }
        None => {
            println!("File parsed incorrectly");
            return OutArr{data:ptr::null::<libc::c_double>() as *mut libc::c_double, rows:0, cols:0};
        }
    }       
}

