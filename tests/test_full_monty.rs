/*extern crate denoise_engine;
use denoise_engine::parse::{NoiseField, SwathElem};
use denoise_engine::read_from_archive::get_data_from_zip_path;
use denoise_engine::apply::{apply_swath_scale, prep_measurement, restore_scale};
use denoise_engine::estimate::*;

use ndarray::{Array1, Array2, Array3, ArrayD};
use ndarray::Zip;
use hdf5;
use hdf5::types::Array;

#[test]
fn test_load_sample(){
    let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    let zipval = get_data_from_zip_path(path);

    let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_arrs.hdf5", "r").unwrap();
    let mut x_gt:Array2<f64> = handle.dataset("x").unwrap().read_2d().unwrap();
    match zipval {
        Some((swath_bounds, w, noisefield, x16)) => {
            println!("Load successful");
            println!("{} {} {} {}", x16.shape()[0], x16.shape()[1], x_gt.shape()[0], x_gt.shape()[1]);
            assert!(x16.shape()[0] == x_gt.shape()[0] && x16.shape()[1]==x_gt.shape()[1]);
        }

        None => {
            panic!("File parsed incorrectly");
        }
    }
}

/// Test the entire pipeline, comparing the final result to the compiled hdarr
#[test]
fn test_full_monty() {
    let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_arrs.hdf5", "r").unwrap();
    let mut x_gt:Array2<f64> = handle.dataset("x").unwrap().read_2d().unwrap();
    let mut y_gt:Array2<f64> = handle.dataset("y").unwrap().read_2d().unwrap();
    let mut denoised:Array2<f32> = handle.dataset("denoised").unwrap().read_2d().unwrap();

    let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
    let lambda2_ = 1.0;
    let mu = 1.7899;
    let gamma = 2.0;
    
    let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    let zipval = get_data_from_zip_path(path);


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
            


            for i in 0..9992 {
                for j in 0..10400 {
                    let a = denoised[(i,j)] as f64;
                    let b = x[(i,j)];
                    if (a-b).abs() >= 40.0 {
                        println!("{} {} : {} {}", i,j,a,b);
                        
                        println!("around python\n{} {} {}\n {} {} {}", denoised[(i,j-1)], denoised[(i,j)], denoised[(i,j+1)], denoised[(i+1,j-1)], denoised[(i+1,j)], denoised[(i+1,j+1)]);
                        println!("around rust\n{} {} {}\n {} {} {}", x[(i,j-1)], x[(i,j)], x[(i,j+1)], x[(i+1,j-1)], x[(i+1,j)], x[(i+1,j+1)]);
                        
                        assert!( (a-b).abs() < 0.2);
                    }
                }
            }

            
        }

        None => {
            panic!("File parsed incorrectly");
        }
    }
    

}

#[test]
fn test_use_k(){
    let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_arrs.hdf5", "r").unwrap();
    let mut x_gt:Array2<f64> = handle.dataset("x").unwrap().read_2d().unwrap();
    let mut y_gt:Array2<f64> = handle.dataset("y").unwrap().read_2d().unwrap();
    let mut denoised:Array2<f32> = handle.dataset("denoised").unwrap().read_2d().unwrap();

    let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
    let lambda2_ = 1.0;
    let mu = 1.7899;
    let gamma = 2.0;
    
    let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    let zipval = get_data_from_zip_path(path);


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
            
            //let k = estimate_k_values(x.view(), noisefield.data.view(), &w, &sb, mu, gamma, lambda_, lambda2_);
            let k = handle.dataset("k").unwrap().read_1d().unwrap();
            println!("k={:?}", k);
            
            {
                let mut y = noisefield.data.view_mut();
                restore_scale(x.view_mut(), y);
            }
            apply_swath_scale(x.view_mut(), noisefield.data.view(), k.view(), &sb);
            


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

        None => {
            panic!("File parsed incorrectly");
        }
    }
    

}



#[test]
fn test_load_x(){
    let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    let zipval = get_data_from_zip_path(path);

    let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_arrs.hdf5", "r").unwrap();
    let mut x_gt:Array2<f64> = handle.dataset("x").unwrap().read_2d().unwrap();

    match zipval {
        Some((swath_bounds, w, noisefield, x16)) => {
            println!("Load successful");
            println!("{} {} {} {}", x16.shape()[0], x16.shape()[1], x_gt.shape()[0], x_gt.shape()[1]);
            assert!(x16.shape()[0] == x_gt.shape()[0] && x16.shape()[1]==x_gt.shape()[1]);
            let mut x:Array2<f64> = Array2::zeros((x16.shape()[0], x16.shape()[1]));
            
            Zip::from(&mut x)
                .and(x16.view())
                .apply(|a,b| *a = (*b as f64)*(*b as f64));


            for i in 0..9992 {
                for j in 0..10400 {
                    let a = x_gt[(i,j)] as f64;
                    let b = x[(i,j)];
                    if (a-b).abs() >= 0.2 {
                        println!("{} {} : {} {}", i,j,a,b);
                        
                        println!("around python\n{} {} {}\n {} {} {}", x_gt[(i,j-1)], x_gt[(i,j)], x_gt[(i,j+1)], x_gt[(i+1,j-1)], x_gt[(i+1,j)], x_gt[(i+1,j+1)]);
                        println!("around rust\n{} {} {}\n {} {} {}", x[(i,j-1)], x[(i,j)], x[(i,j+1)], x[(i+1,j-1)], x[(i+1,j)], x[(i+1,j+1)]);
                        
                        assert!( (a-b).abs() < 0.2);
                    }
                }
            }

        }

        None => {
            panic!("File parsed incorrectly");
        }
    }
}



#[test]
fn test_load_y(){
    let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    let zipval = get_data_from_zip_path(path);

    let handle = hdf5::File::open("/mnt/D2/Data/Sentinel/test_arrs.hdf5", "r").unwrap();
    let mut y_gt:Array2<f64> = handle.dataset("y").unwrap().read_2d().unwrap();

    match zipval {
        Some((swath_bounds, w, noisefield, x16)) => {
            let y = noisefield.data.view();

            for i in 0..9992 {
                for j in 0..10400 {
                    let a = y_gt[(i,j)] as f64;
                    let b = y[(i,j)];
                    if (a-b).abs() >= 0.2 {
                        println!("{} {} : {} {}", i,j,a,b);
                        
                        println!("around python\n{} {} {}\n {} {} {}", y_gt[(i,j-1)], y_gt[(i,j)], y_gt[(i,j+1)], y_gt[(i+1,j-1)], y_gt[(i+1,j)], y_gt[(i+1,j+1)]);
                        println!("around rust\n{} {} {}\n {} {} {}", y[(i,j-1)], y[(i,j)], y[(i,j+1)], y[(i+1,j-1)], y[(i+1,j)], y[(i+1,j+1)]);
                        
                        assert!( (a-b).abs() < 0.2);
                    }
                }
            }

        }

        None => {
            panic!("File parsed incorrectly");
        }
    }
}
*/
