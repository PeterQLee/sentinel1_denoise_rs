///Entry point for denoising zip file from C api
/// Under failure, assume that all elements are null.
use std::ptr;
use s1_noisefloor_engine::parse::{LinearConfig, HyperParams};
use s1_noisefloor_engine::interface;
use s1_noisefloor_engine::postprocess;
use s1_noisefloor_engine::prep_lp;
use ndarray;
use std::sync::Arc;

#[repr(C)]
/// Struct for holding results from Linear results
/// crosspol: pointer to crosspol in f64 (square)
/// copol: pointer to u16
/// k: pointer to linear scales
/// rows: # rows in image
/// cols: # cols in image
pub struct LinearResult {crosspol:*mut f64,
			 copol:*mut u16,
			 k:*mut f64,
                         rows:usize,
                         cols:usize,
                         n_subswaths:usize
}
macro_rules! null_linear {
    () => {
        LinearResult{crosspol:ptr::null::<f64>() as *mut f64,
                     copol:ptr::null::<u16>() as *mut u16,
                     k:ptr::null::<f64>() as *mut f64,
                     rows:0, cols:0, n_subswaths:0}
    }
}


#[repr(C)]
/// struct for holding raw results.
/// crosspol: pointer to crosspol in f64 (square)
/// copol: pointer to u16
/// y: pointer to noise field.
/// rows: # rows in image
/// cols: # cols in image
pub struct RawResult {
    crosspol:*mut u16,
    copol:*mut u16,
    y:*mut f64,
    rows:usize,
    cols:usize
}
macro_rules! null_raw {
    () => {
        RawResult{crosspol:ptr::null::<u16>() as *mut u16,
                  copol:ptr::null::<u16>() as *mut u16,
                  y:ptr::null::<f64>() as *mut f64,
                  rows:0, cols:0}
    }

}

#[repr(C)]
/// Struct for holding results from Linear results
/// crosspol: pointer to crosspol in f64 (square)
/// copol: pointer to u16
/// k: pointer to noise field.
/// m: pointer to slope params
/// b: pointer to intercept params
/// rows: # rows in image
/// cols: # cols in image
/// n_subswaths: number of subswaths in image.
/// plen: number of elements in m.
pub struct LPResult {crosspol:*mut f64,
                     copol:*mut u16,
                     k:*mut f64,
                     m:*mut f64,
                     b:*mut f64,
                     rows:usize,
                     cols:usize,
                     n_subswaths:usize,
                     plen:usize
}
macro_rules! null_lp {
    () => { LPResult{crosspol:ptr::null::<f64>() as *mut f64,
                     copol:ptr::null::<u16>() as *mut u16,
                     k:ptr::null::<f64>() as *mut f64,
                     m:ptr::null::<f64>() as *mut f64,
                     b:ptr::null::<f64>() as *mut f64,
                     rows:0,
                     cols:0,
                     n_subswaths:0,
                     plen:0}
    }
}

#[repr(C)]
/// Simple array structure for multilook processing
/// Assumes Row major order, contiguous.
/// data: pointer to data
/// rows: number of rows in image
/// cols: number of columns in image
pub struct SimpleArray2D {
    data:*mut f64,
    rows:usize,
    cols:usize,
}


#[no_mangle]
/// Applies the lstsquares estimation method to retrieve scaling parameters, k,
/// rescales the noise floor, y, and subtracts it from the image, x.
/// Returns the values back in square intensity units.
/// Parameters:
///
/// path: char*
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// pathlen: size_t
/// configpath: char*
///     Path to the config ini file.
/// configpathlen: size_t
///
/// Returns:
/// Struct of LinearResult. Upon error, all entries will be null.
pub extern fn linear_get_dualpol_data(path:*const u8, pathlen:usize,
                                      configpath:*const u8, configpathlen:usize)
                                      -> LinearResult
{
    // parse the path arguments
    if pathlen == 0 {
        return null_linear!();
    }
    let archpath = std::str::from_utf8(unsafe {std::slice::from_raw_parts(path, (pathlen)*std::mem::size_of::<u8>())}).unwrap();
    let lin_param = match configpathlen {
        0 => LinearConfig::default(),
        _ => {
            let s = std::str::from_utf8(unsafe {std::slice::from_raw_parts(configpath, (configpathlen)*std::mem::size_of::<u8>())}).unwrap();
            match LinearConfig::parse_config(s) {
                Ok(d) => d,
                Err(e) => {
                    eprintln!("Could not parse config {}",e);
                    return null_linear!();
                }
            }
        }
    };

    
    match interface::linear_get_dualpol_data(archpath, &lin_param) {
        Ok((x, co16, k)) => {
            // TODO check x and co16 are contigious
            let nss = k.len();
            let rows = x.nrows();
            let cols = x.ncols();

            // Prevent over-freeing. The calling C-function owns these now.
            let mut x_d = std::mem::ManuallyDrop::new(x.into_raw_vec());
            let mut co16_d = std::mem::ManuallyDrop::new(co16.into_raw_vec());
            let mut k_d = std::mem::ManuallyDrop::new(k.into_raw_vec());
            LinearResult {
                crosspol:x_d.as_mut_ptr(),
                copol:co16_d.as_mut_ptr(),
                k:k_d.as_mut_ptr(),
                rows:rows,
                cols:cols,
                n_subswaths:nss,
            }

        },
        Err(e) => {
            eprintln!("{}",e);
            return LinearResult{crosspol:ptr::null::<libc::c_float>() as *mut f64,
                                copol:ptr::null::<libc::c_float>() as *mut u16,
                                k:ptr::null::<libc::c_float>() as *mut f64,
                                rows:0, cols:0, n_subswaths:0};
        }
    }
}

#[no_mangle]
/// Applies the linear scaling method using custom user provided scales, k.
/// Returns the values back in square intensity units.
///
/// Parameters:
///
/// path: char*
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// pathlen: size_t
/// 
/// kp: *double
///     One dimensional array in 64-bit float that indicates the linear scaling parameters
///     to apply to each subswath. Length of array must equal five.
///     To apply the standard ESA noise removal
/// klen: size_t
///
///
/// Returns:
/// LinearResult
pub extern fn linear_get_customscale_data(path:*const u8, pathlen:usize,
                                          kp:*mut f64, klen:usize) -> LinearResult{
    if pathlen == 0 {
        return null_linear!();
    }
    let archpath = std::str::from_utf8(unsafe {std::slice::from_raw_parts(path, (pathlen)*std::mem::size_of::<u8>())}).unwrap();
    let k = unsafe {ndarray::ArrayView::from_shape_ptr(ndarray::Dim([klen]), kp)};


    match interface::linear_get_customscale_data(archpath, k) {
        Ok((x,co16)) => {
            let rows = x.nrows();
            let cols = x.ncols();
            let mut x_d = std::mem::ManuallyDrop::new(x.into_raw_vec());
            let mut co16_d = std::mem::ManuallyDrop::new(co16.into_raw_vec());
            LinearResult {
                crosspol:x_d.as_mut_ptr(),
                copol:co16_d.as_mut_ptr(),
                k:kp,
                rows:rows,
                cols:cols,
                n_subswaths:klen,
            }
        },
        Err(e) => {
            eprintln!("{}",e);
            return null_linear!();
        }
    }
}

#[no_mangle]
/// Returns the original cross pol, co pol, and noise field from the archive.
/// Note that arrays are in linear units.
///
/// Parameters:
///
/// path: char *
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// pathlen: size_t
///
/// Returns:
/// RawResult
pub extern fn linear_get_raw_data(path:*const u8, pathlen:usize) -> RawResult
{
    if pathlen == 0 {
        return null_raw!();
    }
    let archpath = std::str::from_utf8(unsafe {std::slice::from_raw_parts(path, (pathlen)*std::mem::size_of::<u8>())}).unwrap();

    match interface::linear_get_raw_data(archpath) {
        Ok((x, co16, y)) => {

            let rows = x.nrows();
            let cols = x.ncols();
            let mut x_d = std::mem::ManuallyDrop::new(x.into_raw_vec());
            let mut co16_d = std::mem::ManuallyDrop::new(co16.into_raw_vec());
            let mut y_d = std::mem::ManuallyDrop::new(y.into_raw_vec());
            
            RawResult {
                crosspol:x_d.as_mut_ptr(),
                copol:co16_d.as_mut_ptr(),
                y:y_d.as_mut_ptr(),
                rows:rows,
                cols:cols,
            }
        },
        Err(e) => {
            eprintln!("{}",e);
            return null_raw!();
        }
    }
    
}



#[no_mangle]
/// Applies the linear programming method to restimate a noise floor based
/// on the characteristics of the original image and the 
/// Returns the values back in square intensity units.
///
/// Parameters:
///
/// path: char*
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// pathlen: size_t
/// 
/// lstsq_rescale: unsigned char
///     Indicate whether you want to apply the least squares method from
///     linear_get_dualpol_data to get the baseline minimum offset values for
///     the method. Ignored if product type is IW
///     >0 for applying the method
///     =0 to just use the default ESA noise floor for this.
/// configpath: char *
/// configpathlen: size_t
///     set to 0 if there is no config file and you want to use default settings.
///
/// Returns:
/// LPResult
pub extern fn lp_get_dualpol_data(path:*const u8,
                       pathlen:usize,
                       lstsq_rescale_:u8,
                       configpath:*const u8,
                       configpathlen:usize) -> LPResult
{
    if pathlen == 0 {
        return null_lp!();
    }
    let archpath = std::str::from_utf8(unsafe {std::slice::from_raw_parts(path, (pathlen)*std::mem::size_of::<u8>())}).unwrap();
    let lstsq_rescale:bool = lstsq_rescale_ > 0;
    
    let (lin_param, lp_param) = match configpathlen {
        0 => {println!("Could not parse config path or was not provided. Using default.");
              (LinearConfig::default(), HyperParams::default())},
        
        _ => {
            let s = std::str::from_utf8(unsafe {std::slice::from_raw_parts(configpath, (configpathlen)*std::mem::size_of::<u8>())}).unwrap();
            (match LinearConfig::parse_config(&s) {
                Ok(d) => d,
                Err(e) => {
                    eprintln!("Error parsing config {}",e);
                    return null_lp!();
                }
            }, match HyperParams::parse_config(&s) {
                Ok(d) => d,
                Err(e) => {
                    eprintln!("Error parsing config {}",e);
                    return null_lp!();
                }
            })
        }
    };

    
    match interface::lp_get_dualpol_data(archpath, lstsq_rescale, &lin_param, lp_param) {
        Ok((xv, co16, params, k)) => {

            let xout = Arc::try_unwrap(xv).expect("Could not unwrap");
            //let mut xview = xout.to_ndarray();
            let m:Vec<f64> = params.iter().map(|i| i.iter()
                                               .map(|j| j.m))
                .flatten().collect();
            let b:Vec<f64> = params.iter().map(|i| i.iter()
                                               .map(|j| j.b))
                .flatten().collect();
            
            let num_ss = params.len();
            let plen = m.len();

            let rows = co16.nrows();
            let cols = co16.ncols();

            // Manually drop memory so no double frees occur.
            let mut xout_d = std::mem::ManuallyDrop::new(xout);
            let mut m_d = std::mem::ManuallyDrop::new(m);
            let mut b_d = std::mem::ManuallyDrop::new(b);
            let mut k_d = std::mem::ManuallyDrop::new(k.into_raw_vec());
            let mut co16_d = std::mem::ManuallyDrop::new(co16.into_raw_vec());
            let ptr = xout_d.to_raw_pointer();

            
            LPResult {
                crosspol:ptr,
                copol:co16_d.as_mut_ptr(),
                k:k_d.as_mut_ptr(),
                m:m_d.as_mut_ptr(),
                b:b_d.as_mut_ptr(),
                rows:rows,
                cols:cols,
                n_subswaths:num_ss,
                plen:plen
            }
            
            
            
        },
        Err(e) => {
            eprintln!("{}",e);
            return null_lp!();
        }
    }
}


#[no_mangle]
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
pub extern fn lp_get_customscale_data(path:*const u8,
                           pathlen:usize,
                           m:*mut f64,
                           b:*mut f64,
                           plen:usize,
                           configpath:*const u8,
                           configpathlen:usize) -> LPResult
{
    if pathlen == 0 {
        return null_lp!();
    }
    let archpath = std::str::from_utf8(unsafe {std::slice::from_raw_parts(path, (pathlen)*std::mem::size_of::<u8>())}).unwrap();


    let lp_param = match configpathlen {
        0 => {println!("Could not parse config path or was not provided. Using default.");
              HyperParams::default()},
        _ => {
            let s = std::str::from_utf8(unsafe {std::slice::from_raw_parts(configpath, (configpathlen)*std::mem::size_of::<u8>())}).unwrap();
            match HyperParams::parse_config(&s) {
                Ok(d) => d,
                Err(e) => {
                    eprintln!("Error parsing config {}",e);
                    return null_lp!();
                }
            }
        }
    };


    let num_subswaths:usize = match plen {
        12 => 5,
        6 => 3,
        _ => {
            eprintln!("Invalid number of parameters");
            return null_lp!();
        }
    };
    let m_s = unsafe{std::slice::from_raw_parts(m, (plen)*std::mem::size_of::<f64>())};
    let b_s = unsafe{std::slice::from_raw_parts(b, (plen)*std::mem::size_of::<f64>())};
    match interface::lp_get_customscale_data(archpath, m_s, b_s, num_subswaths, lp_param) {
        Ok((xv, co16)) => {
            let xout = Arc::try_unwrap(xv).expect("Could not unwrap");
            let rows = xout.rows;
            let cols = xout.cols;
            let mut xout_d = std::mem::ManuallyDrop::new(xout);
            let mut co16_d = std::mem::ManuallyDrop::new(co16.into_raw_vec());
            let ptr = xout_d.to_raw_pointer();
            
            LPResult {
                crosspol:ptr,
                copol:co16_d.as_mut_ptr(),
                k:std::ptr::null::<f64>() as *mut f64,
                m:m,
                b:b,
                rows:rows,
                cols:cols,
                n_subswaths:0,
                plen:plen,
            }

        },
        
        Err(e) => {
            eprintln!("An error occurred. No output written: {}",e);
            return null_lp!();
        }
    }
    

}

#[no_mangle]
/// Applies multilooking to the input image, sets negative values to 0, and
/// square roots the output values.
///
/// Parameters:
/// x: Input array for multilooking (square units) - SimpleArray2D
/// row_factor: integer amount to mean reduce row by.
/// col_factor: integer amount to mean reduce col by.
/// num_cores: number of multithreading cores to use.
///
/// Output:
/// x : multilooked image (linear units) - SimpleArray2D
pub extern fn post_multilook_and_floor(x:*mut SimpleArray2D,
			    row_factor:usize,
			    col_factor:usize,
			    num_cores:usize) -> SimpleArray2D {
    let x_s = unsafe{std::slice::from_raw_parts((*x).data, ((*x).rows*(*x).cols)*std::mem::size_of::<f64>())};
    let x_v = x_s.to_vec();
    let xp = Arc::new(prep_lp::TwoDArray::from_vec(x_v,
						   unsafe{(*x).rows}, unsafe{(*x).cols}));
    let o = postprocess::multilook_and_floor(xp,
					     row_factor,
					     col_factor,
					     num_cores );
    let (orow, ocol) = (o.rows, o.cols);

    
    let mut o_d = std::mem::ManuallyDrop::new(o);
    let op = o_d.to_raw_pointer();
    SimpleArray2D {
	data:op,
	rows:orow,
	cols:ocol
    }

	
}
