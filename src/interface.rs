extern crate denoise_engine;
use denoise_engine::parse::{NoiseField, SwathElem};
use denoise_engine::read_from_archive::get_data_from_zip_path;
use denoise_engine::apply::{apply_swath_scale, prep_measurement, restore_scale};
use denoise_engine::estimate::*;
use ndarray::{Array1, Array2};

#[repr(C)] pub struct OutArr {data:*libc::c_double, rows:libc::c_int, cols:libc::c_int}


#[no mangle]
pub extern fn denoise_zip(path:*const u8, pathlen:libc::c_int) -> OutArr{
    
    let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
    let lambda2_ = 1.0;
    let mu = 1.7899;
    let gamma = 2.0;
    let zipval = get_data_from_zip_path(zip_path);

    
    match zipval {
        Some((swath_bounds, w, mut noisefield, x16)) => {
            println!("Load successful");
            let mut x_:Option<Array2<f64>> = None;
            {
                let mut y = noisefield.data.view_mut();
                x_ = Some(prep_measurement(x16.view(), y));
            }
            let mut x = x_.unwrap();

            let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
            
            let k = estimate_k_values(x.view(), noisefield.data.view(), &w, &sb, mu, gamma, lambda_, lambda2_);
            println!("k={:?}", k);
            
            {
                let mut y = noisefield.data.view_mut();
                restore_scale(x.view_mut(), y);
            }
            apply_swath_scale(x.view_mut(), noisefield.data.view(), k.view(), &sb);

            let (rows, cols) = x.dims();
            assert!(x.is_standard_layout());
            //x.as_mut_ptr() ? however don't know if this will satisfy borrow checker
            
        }
        None => {
            panic!("File parsed incorrectly");
        }
    }       
}

