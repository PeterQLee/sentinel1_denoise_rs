//! Python interface to methods.
/*
MIT License

Copyright (c) 2020 Peter Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

use s1_noisefloor_engine::parse::{LinearConfig, HyperParams};
use s1_noisefloor_engine::interface;
use s1_noisefloor_engine::postprocess;
use s1_noisefloor_engine::prep_lp;


use numpy::{PyArray, PyArray1, PyArray2};
use pyo3::prelude::{*};
use pyo3::{wrap_pyfunction, exceptions};
use pyo3::types::{PyList, PyString, PyAny};


use std::sync::Arc;



/// Noise floor removal engine for Sentinel-1
/// Two types of noise removal methods are available.
/// Note that the engine currently only removes noise floor from cross-pol images.
///
/// 1. A linear noise floor removal method that rescales the ESA provided noise floor
///    that is provided in each Sentinel-1 product. This is the application of the method in
///    P. Q. Lee, L. Xu, D. A. Clausi. 2020. Sentinel-1 additive noise removal from cross-polarization extra-wide TOPSAR with dynamic least-squares. Remote Sensing of Environment. 248. https://doi.org/10.1016/j.rse.2020.111982
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
fn s1_noisefloor(_py: Python, m:&PyModule) -> PyResult<()> {
    /// Applies the lstsquares estimation method to retrieve scaling parameters, k,
    /// rescales the noise floor, y, and subtracts it from the image, x.
    /// Returns the values back in square intensity units.
    ///
    /// Parameters:
    ///
    /// archpath: str
    ///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
    /// config_path: str or None
    ///     Optional path to config file. If None (or non-string) will use default configuration.
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
    /// k: ndarray(1)
    ///    array of the estimated linear parameters
    #[pyfn(m, "linear_get_dualpol_data")]
    fn linear_get_dualpol_data(__py:Python, archpath:&str, config_path:&PyAny) -> PyResult<(Py<PyArray2<f64>>,Py<PyArray2<u16>>, Py<PyArray1<f64>>)> {
	let cfgpath:PyResult<_> = PyString::from_object(
	    config_path,
	    "utf-8",
	    "");
	let lin_param = match cfgpath {
	    Ok(pth) => {
		let s = pth.to_string_lossy();
		match LinearConfig::parse_config(&s) {
		    Ok(d) => d,
		    Err(e) => {
			eprintln!("Could not parse config {}",e);
			return exceptions::ValueError.into();
		    }
		}
	    }
	    Err(_e) => {
		println!("Could not parse config path or was not provided. Using default.");
		LinearConfig::default()
	    }
	};
	match interface::linear_get_dualpol_data(archpath, &lin_param) {
	    Ok((x, co16, k)) => {
		let py_cross = PyArray::from_array(__py,&x).to_owned();
		let py_co = PyArray::from_array(__py,&co16).to_owned();
		let py_k = PyArray::from_array(__py, &k).to_owned();
		return Ok((py_cross, py_co, py_k));
	    },
	    Err(e) => {
		eprintln!("{}",e);
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
    /// archpath: str
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
    fn linear_get_customscale_data(__py:Python, archpath:&str, py_k:&PyArray1<f64>) -> PyResult<(Py<PyArray2<f64>>,Py<PyArray2<u16>>)> {
	let k = py_k.as_array();

	match interface::linear_get_customscale_data(archpath, k) {
	    Ok((x, co16)) => {
                let py_cross = PyArray::from_array(__py,&x).to_owned();
                let py_co = PyArray::from_array(__py,&co16).to_owned();
                return Ok((py_cross, py_co));

            },
	    Err(e) => {
		eprintln!("{}",e);
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
    /// archpath: str
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
    fn linear_get_raw_data(__py:Python, archpath:&str) -> PyResult<(Py<PyArray2<u16>>, Py<PyArray2<u16>>, Py<PyArray2<f64>>)> {
        match interface::linear_get_raw_data(archpath) {
            Ok((x, co16, y)) => {
                let py_cross = PyArray::from_array(__py,&x).to_owned();
                let py_co = PyArray::from_array(__py,&co16).to_owned();
                let py_y = PyArray::from_array(__py,&y).to_owned();
                return Ok((py_cross, py_co, py_y));

            },
            Err(e) => {
                eprintln!("{}",e);
                return exceptions::ValueError.into();
            }
        }

    }

    /// Applies the linear programming method to restimate a noise floor based
    /// on the characteristics of the original image and the 
    /// Returns the values back in square intensity units.
    ///
    /// Parameters:
    ///
    /// archpath: str
    ///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
    /// lstsq_rescale: bool
    ///     Indicate whether you want to apply the least squares method from
    ///     linear_get_dualpol_data to get the baseline minimum offset values for
    ///     the method. Ignored if product type is 
    ///     true for applying the method
    ///     false to just use the default ESA noise floor for this.
    /// config_path: str or None
    ///     Optional path to config file. If None (or non-string) will use default configuration.
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
    /// k: ndarray(1)
    ///    Array of linear scales computed by prior to applying the LP.
    #[pyfn(m, "lp_get_dualpol_data")]
    fn lp_get_dualpol_data<'p>(__py:Python<'p>,
                               archpath:&str,
                               lstsq_rescale:bool,
                               config_path:&PyAny)
                           -> PyResult<(Py<PyArray2<f64>>,
                                        Py<PyArray2<u16>>,
                                        &'p PyList,
                                        &'p PyList,
                                        Py<PyArray1<f64>>
                               )>
    {
        let cfgpath:PyResult<String> = config_path.extract();
        let (lin_param, lp_param) = match cfgpath {
            Ok(s) => {
                (match LinearConfig::parse_config(&s) {
                    Ok(d) => d,
                    Err(e) => {
                        println!("Error parsing config {}",e);
                        return exceptions::ValueError.into();
                    }
                }, match HyperParams::parse_config(&s) {
                    Ok(d) => d,
                    Err(e) => {
                        println!("Error parsing config {}",e);
                        return exceptions::ValueError.into();
                    }
                })
            }
            Err(_e) => {
                println!("Could not parse config path or was not provided. Using default.");
                (LinearConfig::default(), HyperParams::default())
            }
        };

            
       match interface::lp_get_dualpol_data(archpath, lstsq_rescale, &lin_param, lp_param) {
           Ok((xv, co16, params, k)) => {

               let xout = Arc::try_unwrap(xv).expect("Could not unwrap");
               let xview = xout.to_ndarray();
               let py_cross = PyArray::from_array(__py, &xview).to_owned();
               let py_co = PyArray::from_array(__py, &co16).to_owned();
               let py_k = PyArray::from_array(__py, &k).to_owned();

               let m:&PyList = PyList::new(__py, params.iter()
                                           .map(|i| i.iter()
                                                .map(|j| j.m))
                                           .flatten().collect::<Vec<f64>>());
               let b:&PyList = PyList::new(__py, params.iter()
                                           .map(|i| i.iter()
                                                .map(|j| j.b))
                                           .flatten().collect::<Vec<f64>>());

               
               return Ok((py_cross, py_co, m, b, py_k));
               
            },
            Err(e) => {
                eprintln!("{}",e);
                return exceptions::ValueError.into();
            }
        }
    }
    ///  Applies the power function noise floor obtained from linear programming,
    ///   with parameters given by the user.
    /// Returns the values back in square intensity units.
    /// 
    /// Parameters:
    ///
    /// archpath: str
    ///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
    /// m: list
    ///     List of slope parameters
    /// b: list
    ///     List of intercept parameters.
    /// config_path: str or None
    ///     Optional path to config file. If None (or non-string) will use default configuration.
    #[pyfn(m, "lp_get_customscale_data")]
    fn lp_get_customscale_data<'p>(__py:Python<'p>, archpath:&str,
                                   m:Vec<f64>,
                                   b:Vec<f64>,
                                   config_path:&PyAny)
                            -> PyResult<(Py<PyArray2<f64>>,
                                        Py<PyArray2<u16>>)>
    {
        let cfgpath:PyResult<String> = config_path.extract();
        let lp_param = match cfgpath {
            Ok(s) => {
                match HyperParams::parse_config(&s) {
                    Ok(d) => d,
                    Err(e) => {
                        println!("Error parsing config {}",e);
                        return exceptions::ValueError.into();
                    }
                }
            },
            Err(_e) => {
                println!("Could not parse config path or was not provided. Using default.");
                HyperParams::default()
            }
        };
        if m.len() != b.len() {
            eprintln!("Invalid number of parameters: len(m)!=len(b)");
            return exceptions::ValueError.into();

        }
        let num_subswaths:usize = match m.len() {
            12 => 5,
            6 => 3,
            _ => {
                eprintln!("Invalid number of parameters");
                return exceptions::ValueError.into();
            }
        };
        match interface::lp_get_customscale_data(archpath, &m, &b, num_subswaths, lp_param) {
            Ok((xv, co16)) => {
                let xout = Arc::try_unwrap(xv).expect("Could not unwrap");
                let xview = xout.to_ndarray();
                let py_cross = PyArray::from_array(__py, &xview).to_owned();
                let py_co = PyArray::from_array(__py, &co16).to_owned();
                return Ok((py_cross, py_co))
            },
        
            Err(e) => {
                eprintln!("An error occurred. No output written: {}",e);                
            }
        }
            
        return exceptions::ValueError.into();
    }

    /// Applies multilooking to the input image, sets negative values to 0, and
    /// square roots the output values.
    ///
    /// Parameters:
    /// py_x: Input array for multilooking (square units)
    /// row_factor: integer amount to mean reduce row by.
    /// col_factor: integer amount to mean reduce col by.
    /// num_cores: number of multithreading cores to use.
    ///
    /// Output:
    /// x : multilooked image (linear units)
    #[pyfn(m, "post_multilook_and_floor")]
    fn post_multilook_and_floor(__py:Python, py_x:&PyArray2<f64>,
                                row_factor:u64,
                                col_factor:u64,
                                num_cores:u64) -> PyResult<Py<PyArray2<f64>>> {
        let x = py_x.to_owned_array();
        let xp = Arc::new(prep_lp::TwoDArray::from_ndarray(x));
        let o = postprocess::multilook_and_floor(xp,
                                                 row_factor as usize,
                                                 col_factor as usize,
                                                 num_cores as usize);
        let op = o.to_ndarray();
        let py_o = PyArray::from_array(__py, &op).to_owned();
        Ok(py_o)
        
    }


    m.add_wrapped(wrap_pyfunction!(linear_get_dualpol_data))?;
    m.add_wrapped(wrap_pyfunction!(linear_get_customscale_data))?;
    m.add_wrapped(wrap_pyfunction!(linear_get_raw_data))?;
    m.add_wrapped(wrap_pyfunction!(lp_get_dualpol_data))?;
    m.add_wrapped(wrap_pyfunction!(lp_get_customscale_data))?;

    m.add_wrapped(wrap_pyfunction!(post_multilook_and_floor))?;
    Ok(())
}
