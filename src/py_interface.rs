use crate::parse::{SwathElem};
use crate::read_from_archive::{get_data_from_zip_path, SentinelArchiveOutput};
use crate::apply::{apply_swath_scale, prep_measurement, restore_scale, convert_to_f64_f32, convert_to_u16_f32};
use crate::estimate::*;
extern crate libc;
use ndarray::{Array2, arr1};
use std::ptr;
use numpy::{PyArray, PyArray1, PyArray2};
use pyo3::prelude::{Py, pymodule,  PyModule, PyResult, Python};
use pyo3::wrap_pyfunction;
use pyo3::exceptions;
//extern crate openblas_src;
extern crate lapack_src;
//extern crate lapacke;


#[repr(C)] pub struct OutArr {crosspol:*mut libc::c_float,
                              copol:*mut libc::c_float,
                              rows:libc::c_int,
                              cols:libc::c_int}

///Entry point for denoising zip file from C api
#[no_mangle]
pub extern fn denoise_zip(path:*const u8, pathlen:libc::c_int) -> OutArr{
    
    let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
    let lambda2_ = 1.0;
    let mu = 1.7899;
    let gamma = 2.0;
    // parse the path arguments
    let zip_path = unsafe {std::slice::from_raw_parts(path, (pathlen as usize)*std::mem::size_of::<u8>())};
    let zipval = get_data_from_zip_path(std::str::from_utf8(zip_path).unwrap(), true);

    
    match zipval {
        Ok(archout) => {
            match archout {
                SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, _lpargs) => {
                    let mut x_:Option<Array2<f64>> = None;
                    {
                        let mut y = noisefield.data.view_mut();
                        x_ = Some(prep_measurement(x16.view(), y));
                    }
                    let mut x = x_.unwrap();

                    let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
                    
                    let k = estimate_k_values(x.view(), noisefield.data.view(), &w, &sb, mu, gamma, lambda_, lambda2_);
                    
                    //{
                     //   let mut y = noisefield.data.view_mut();
                        //restore_scale(x.view_mut(), y);
                    //}
                    apply_swath_scale(x.view_mut(), noisefield.data.view(), k.view(), &sb);

                    restore_scale(x.view_mut());

                    let (rows, cols) = x.dim();
                    assert!(x.is_standard_layout());
                    let mut crosspol = convert_to_f64_f32(x.view());
                    let mut copol = convert_to_u16_f32(co16.view());
                    
                    let result = OutArr {
                        crosspol:crosspol.as_mut_ptr(),
                        copol:copol.as_mut_ptr(),
                        rows:x.shape()[0] as libc::c_int,
                        cols:x.shape()[1] as libc::c_int
                    };

                    std::mem::forget(crosspol);
                    std::mem::forget(copol);

                    return result;
                },
            
                _ => { println!("File parsed incorrectly");
                       return OutArr{crosspol:ptr::null::<libc::c_float>() as *mut libc::c_float,
                                     copol:ptr::null::<libc::c_float>() as *mut libc::c_float,
                                     rows:0, cols:0};}
            }
        },
        Err(e) => {
            println!("File parsed incorrectly: {}",e);
            return OutArr{crosspol:ptr::null::<libc::c_float>() as *mut libc::c_float,
                          copol:ptr::null::<libc::c_float>() as *mut libc::c_float,
                          rows:0, cols:0};
        }
    }       
}


/// Noise floor removal engine for Sentinel-1
/// Two types of noise removal methods are available.
/// Note that the engine currently only removes noise floor from cross-pol images.
///
/// 1. A linear noise floor removal method that rescales the ESA provided noise floor
///    that is provided in each Sentinel-1 product. This is the application of the method in
///    P.Q. Lee, L. Xu, D.A. Clausi. 2020. Sentinel-1 additive noise removal from cross-polarization extra-wide TOPSAR with dynamic least-squares. https://doi.org/10.1016/j.rse.2020.111982
///    Currently only available in GRD EW mode images.
///    These methods are prefaced by linear_[...]
///
/// 2. A non-linear noise floor removal method that computes the noise floor as a power function of
///    the antenna pattern. Parameters are estimated with linear programming. This is the application
///    of the method in [].
///    Can be applied to both EW and IW GRD mode images.
///    These methods are prefaced by lp_[...]
#[pymodule]
fn denoise_engine(_py: Python, m:&PyModule) -> PyResult<()> {
    /// Applies the lstsquares estimation method to retrieve scaling parameters, k
    /// and 
    #[pyfn(m, "get_dualpol_data")]
    fn linear_get_dualpol_data(__py:Python, zippath:&str) -> PyResult<(Py<PyArray2<f64>>,Py<PyArray2<u16>>, Py<PyArray1<f64>>)> {
        let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
        let lambda2_ = 1.0;
        let mu = 1.7899;
        let gamma = 2.0;
        let zipval = get_data_from_zip_path(zippath, true);
        match zipval {
            Ok(archout) => {
                match archout {
                    SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, _lpargs) => {
                        let mut x_:Option<Array2<f64>> = None;
                        {
                            let mut y = noisefield.data.view_mut();
                            x_ = Some(prep_measurement(x16.view(), y));
                        }
                        let mut x = x_.unwrap();
                        
                        let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
                        
                        let k = estimate_k_values(x.view(), noisefield.data.view(), &w, &sb, mu, gamma, lambda_, lambda2_);
                        
                        apply_swath_scale(x.view_mut(), noisefield.data.view(), k.view(), &sb);
                        
                        restore_scale(x.view_mut());


                        let py_cross = PyArray::from_array(__py,&x).to_owned();
                        let py_co = PyArray::from_array(__py,&co16).to_owned();
                        let py_k = PyArray::from_array(__py, &k).to_owned();
                        return Ok((py_cross, py_co, py_k));

                    },
                    _ => {}

                }
                
            },
            Err(e) => {println!("{}",e);}
        }

        return exceptions::ValueError.into();
        //return Err(exceptions::ValueError("bad parsing path"));
    }


    #[pyfn(m, "get_noscale_data")]
    fn get_noscale_data(__py:Python, zippath:&str) -> PyResult<(Py<PyArray2<f64>>,Py<PyArray2<u16>>, Py<PyArray1<f64>>)> {
        let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
        let lambda2_ = 1.0;
        let mu = 1.7899;
        let gamma = 2.0;
        let zipval = get_data_from_zip_path(zippath, true);
        match zipval {
            Ok(archout) => {
                match archout {
                    SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, _lpargs) => {
                        let mut x_:Option<Array2<f64>> = None;
                        {
                            let mut y = noisefield.data.view_mut();
                            x_ = Some(prep_measurement(x16.view(), y));
                        }
                        let mut x = x_.unwrap();
                        
                        let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
                        
                        //let k = estimate_k_values(x.view(), noisefield.data.view(), &w, &sb, mu, gamma, lambda_, lambda2_);
                        let k = arr1(&[1.0,1.0,1.0,1.0,1.0]);
                        
                        apply_swath_scale(x.view_mut(), noisefield.data.view(), k.view(), &sb);
                        
                        restore_scale(x.view_mut());


                        let py_cross = PyArray::from_array(__py,&x).to_owned();
                        let py_co = PyArray::from_array(__py,&co16).to_owned();
                        let py_k = PyArray::from_array(__py, &k).to_owned();
                        return Ok((py_cross, py_co, py_k));

                    },
                    _ => {}

                }
                
            },
           Err(e) => {println!("{}",e);}
        }

        return exceptions::ValueError.into();
        //return Err(exceptions::ValueError("bad parsing path"));
    }

    #[pyfn(m, "get_noise_data")]
    fn get_noise_data(__py:Python, zippath:&str) -> PyResult<Py<PyArray2<f64>>> {
        let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
        let lambda2_ = 1.0;
        let mu = 1.7899;
        let gamma = 2.0;
        let zipval = get_data_from_zip_path(zippath, true);
        match zipval {
            Ok(archout) => {
                match archout {
                    SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, _lpargs) => {
                        
                        let mut y = noisefield.data.view_mut();

                        let py_noise = PyArray::from_array(__py,&y).to_owned();
                        return Ok(py_noise);

                    },
                    _ => {}

                }
                
            },
            Err(e) => {println!("{}",e);}
        }

        return exceptions::ValueError.into();
        //return Err(exceptions::ValueError("bad parsing path"));
    }
    
    #[pyfn(m, "get_raw_crosspol")]
    /// Deprecated: Returns the un denoised crosspol
    fn get_raw_crosspol(__py:Python, zippath:&str) -> PyResult<Py<PyArray2<f64>>> {
        let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
        let lambda2_ = 1.0;
        let mu = 1.7899;
        let gamma = 2.0;
        let zipval = get_data_from_zip_path(zippath, true);
        match zipval {
            Ok(archout) => {
                match archout {
                    SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, _lpargs) => {
                        
                        let mut x_:Option<Array2<f64>> = None;
                        {
                            let mut y = noisefield.data.view_mut();
                            x_ = Some(prep_measurement(x16.view(), y));
                        }
                        let mut x = x_.unwrap();


                        let py_cross = PyArray::from_array(__py,&x).to_owned();
                        return Ok(py_cross);

                    },
                    _ => {}

                }
                
            },
            Err(e) => {println!("{}",e);}
        }

        return exceptions::ValueError.into();
        //return Err(exceptions::ValueError("bad parsing path"));
    }

    #[pyfn(m, "get_raw_data")]
    ///Returns the original cross pol, co pol, and noise field.
    fn get_raw_data(__py:Python, zippath:&str) -> PyResult<(Py<PyArray2<u16>>, Py<PyArray2<u16>>, Py<PyArray2<f64>>)> {
        let lambda_ = &[0.1,0.1,6.75124,2.78253,10.0]; //convert to array
        let lambda2_ = 1.0;
        let mu = 1.7899;
        let gamma = 2.0;
        let zipval = get_data_from_zip_path(zippath, true);
        match zipval {
            Ok(archout) => {
                match archout {
                    SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, _lpargs) => {                 
                        let y = noisefield.data.view_mut();
                        
                        let py_cross = PyArray::from_array(__py,&x16).to_owned();
                        let py_co = PyArray::from_array(__py,&co16).to_owned();
                        let py_y = PyArray::from_array(__py, &y).to_owned();
                        return Ok((py_cross, py_co, py_y));

                    },
                    _ => {}

                }
                
            },
            Err(e) => {println!("{}",e);}
        }

        return exceptions::ValueError.into();
        //return Err(exceptions::ValueError("bad parsing path"));

    }


    m.add_wrapped(wrap_pyfunction!(get_noscale_data))?;
    m.add_wrapped(wrap_pyfunction!(get_dualpol_data))?;
    m.add_wrapped(wrap_pyfunction!(get_noise_data))?;
    m.add_wrapped(wrap_pyfunction!(get_raw_crosspol))?;
    Ok(())
}
