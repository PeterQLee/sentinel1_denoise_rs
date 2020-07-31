
use crate::parse::{SwathElem,  HyperParams, LinearConfig};
use crate::interface;
extern crate libc;
// use ndarray::{Array2, arr1};
// use std::ptr;
use numpy::{PyArray, PyArray1, PyArray2};
use pyo3::prelude::{Py, pymodule,  PyModule, PyResult, Python};
use pyo3::{wrap_pyfunction, exceptions};
use pyo3::types::PyList;

use std::sync::Arc;
//extern crate openblas_src;
extern crate lapack_src;
//extern crate lapacke;
// use crate::prep_lp::*;
/*
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
                    let y = noisefield.data.view_mut();
                    let mut x = prep_measurement(x16.view(), y);


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
                        rows:rows as libc::c_int,
                        cols:cols as libc::c_int
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
*/

/// Noise floor removal engine for Sentinel-1
/// Two types of noise removal methods are available.
/// Note that the engine currently only removes noise floor from cross-pol images.
///
/// 1. A linear noise floor removal method that rescales the ESA provided noise floor
///    that is provided in each Sentinel-1 product. This is the application of the method in
///    P.Q. Lee, L. Xu, D.A. Clausi. 2020. Sentinel-1 additive noise removal from cross-polarization extra-wide TOPSAR with dynamic least-squares. https://doi.org/10.1016/j.rse.2020.111982
///    Currently only available in GRD EW mode images.
///    These methods are prefaced by linear_[...].
///
/// 2. A non-linear noise floor removal method that computes the noise floor as a power function of
///    the antenna pattern. Parameters are estimated with linear programming. This is the application
///    of the method in [].
///    Can be applied to both EW and IW GRD mode images.
///    These methods are prefaced by lp_[...].
///
/// Unless stated otherwise, all methods will return the measurements in intensity (square) units.
/// This is in case one wants to perform more significant analysis on the image directly after noise floor removal.
/// Typically one wants to do analysis in amplitude (linear) units. An issue creeping up from this is
/// the potential of applying the square root to potentially negative measurements.
/// To accomodate this, we provide postprocessing methods to handle this aspect
/// with the methods prefaced by post_[...].
#[pymodule]
fn denoise_engine(_py: Python, m:&PyModule) -> PyResult<()> {
    /// Applies the lstsquares estimation method to retrieve scaling parameters, k,
    /// rescales the noise floor, y, and subtracts it from the image, x.
    /// Returns the values back in square intensity units.
    /// Parameters:
    ///
    /// zippath: str
    ///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
    ///
    ///
    /// Returns:
    /// (cross, co, k)
    ///
    /// cross: ndarray(2)
    ///    array holding processed cross-polarized measurements (square)
    ///
    /// co: ndarray(2)
    ///    array holding unprocessed co-polarized measurements (linear)
    ///
    /// k: ndarray(2)
    ///    array of the estimated linear parameters
    #[pyfn(m, "linear_get_dualpol_data")]
    fn linear_get_dualpol_data(__py:Python, zippath:&str) -> PyResult<(Py<PyArray2<f64>>,Py<PyArray2<u16>>, Py<PyArray1<f64>>)> {
	match interface::linear_get_dualpol_data(zippath, &LinearConfig::default()) {
	    Ok((x, co16, k)) => {
		let py_cross = PyArray::from_array(__py,&x).to_owned();
		let py_co = PyArray::from_array(__py,&co16).to_owned();
		let py_k = PyArray::from_array(__py, &k).to_owned();
		return Ok((py_cross, py_co, py_k));
	    },
	    Err(e) => {
		println!("{}",e);
		return exceptions::ValueError.into();
	    }
	}
    }


    #[pyfn(m, "linear_get_customscale_data")]
    /// Applies the linear scaling method using custom user provided scales, k.
    /// Returns the values back in square intensity units.
    ///
    /// Parameters:
    ///
    /// zippath: str
    ///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
    /// k: ndarray(1)
    ///     One dimensional array in 64-bit float that indicates the linear scaling parameters
    ///     to apply to each subswath. Length of array must equal five.
    ///     To apply the standard ESA noise removal, make k = np.array([1.0,1.0,1.0,1.0,1.0])
    ///
    ///
    /// Returns:
    /// (cross, co)
    ///
    /// cross: ndarray(2)
    ///    array holding processed cross-polarized measurements (square)
    ///
    /// co: ndarray(2)
    ///    array holding unprocessed co-polarized measurements (linear)
    ///
    fn linear_get_customscale_data(__py:Python, zippath:&str, py_k:&PyArray1<f64>) -> PyResult<(Py<PyArray2<f64>>,Py<PyArray2<u16>>)> {
	let k = py_k.as_array();
	match interface::linear_get_customscale_data(zippath, k, &LinearConfig::default()) {
	    Ok((x, co16)) => {
                let py_cross = PyArray::from_array(__py,&x).to_owned();
                let py_co = PyArray::from_array(__py,&co16).to_owned();
                return Ok((py_cross, py_co));

            },
	    Err(e) => {
		println!("{}",e);
		return exceptions::ValueError.into();
	    }
	}
    }
    

    #[pyfn(m, "linear_get_raw_data")]
    /// Returns the original cross pol, co pol, and noise field from the archive.
    /// Note that arrays are in linear units.
    ///
    /// Parameters:
    ///
    /// zippath: str
    ///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
    ///
    /// Returns:
    /// (cross, co, y)
    ///
    /// cross: ndarray(2)
    ///    array holding processed cross-polarized measurements (linear)
    /// co: ndarray(2)
    ///    array holding unprocessed co-polarized measurements (linear)
    /// y: ndarray(2)
    ///    array holding noise field (linear)
    fn linear_get_raw_data(__py:Python, zippath:&str) -> PyResult<(Py<PyArray2<u16>>, Py<PyArray2<u16>>, Py<PyArray2<f64>>)> {
	match interface::linear_get_raw_data(zippath) {
	    Ok((x, co16, y)) => {
                let py_cross = PyArray::from_array(__py,&x).to_owned();
                let py_co = PyArray::from_array(__py,&co16).to_owned();
		let py_y = PyArray::from_array(__py,&y).to_owned();
                return Ok((py_cross, py_co, py_y));

            },
	    Err(e) => {
		println!("{}",e);
		return exceptions::ValueError.into();
	    }
	}

    }

    /// Applies the linear programming method to restimate a noise floor based
    /// on the characteristics of the original image and the 
    /// Returns the values back in square intensity units.
    /// Parameters:
    ///
    /// zippath: str
    ///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
    /// lstsq_rescale: bool
    ///     Indicate whether you want to apply the least squares method from
    ///     linear_get_dualpol_data to get the baseline minimum offset values for
    ///     the method. Ignored if product type is 
    ///     true for applying the method
    ///     false to just use the default ESA noise floor for this.
    ///
    /// Returns:
    /// (cross, co, m, v)
    ///
    /// cross: ndarray(2)
    ///    array holding processed cross-polarized measurements (square)
    /// co: ndarray(2)
    ///    array holding unprocessed co-polarized measurements (linear)
    /// m: ndarray(1)
    ///    Array of slope / exponent parameters estimated.
    /// b: ndarray(1)
    ///    Array of intercept parameters estimated
    #[pyfn(m, "lp_get_dualpol_data")]
    fn lp_get_dualpol_data<'p>(__py:Python<'p>, zippath:&str, lstsq_rescale:bool)
			   -> PyResult<(Py<PyArray2<f64>>,
					Py<PyArray2<u16>>,
					&'p PyList,
					&'p PyList
    )>{
       match interface::lp_get_dualpol_data(zippath, lstsq_rescale, &LinearConfig::default()) {
	   Ok((xv, co16, params)) => {

	       let xout = Arc::try_unwrap(xv).expect("Could not unwrap");
	       let xview = xout.to_ndarray();
               let py_cross = PyArray::from_array(__py, &xview).to_owned();
               let py_co = PyArray::from_array(__py, &co16).to_owned();

	       
	       // convert lp params to vectors and outputs
	       let m:&PyList = PyList::new(__py, params.iter()
					   .map(|i| PyList::new(__py, i.iter()
								.map(|j| j.m))));
	       let b:&PyList = PyList::new(__py, params.iter()
			    .map(|i| PyList::new(__py, i.iter()
						 .map(|j| j.b))));
	       
	       
	       return Ok((py_cross, py_co, m, b));
	       
            },
	    Err(e) => {
		println!("{}",e);
		return exceptions::ValueError.into();
	    }
	}
    }

    m.add_wrapped(wrap_pyfunction!(linear_get_dualpol_data))?;
    m.add_wrapped(wrap_pyfunction!(linear_get_customscale_data))?;
    m.add_wrapped(wrap_pyfunction!(linear_get_raw_data))?;
    m.add_wrapped(wrap_pyfunction!(lp_get_dualpol_data))?;
    Ok(())
}
