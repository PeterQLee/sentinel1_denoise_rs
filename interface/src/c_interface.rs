///Entry point for denoising zip file from C api

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
#[repr(C)]
/// Struct for holding results from Linear results
pub struct LinearResult {crosspol:*mut f64,
				    copol:*mut u16,
				    k:*mut f64,
				    rows:libc::c_int,
				    cols:libc::c_int,
				    n_subswaths:libc::c_int
}


#[no_mangle]
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
/// Struct of LinearResult. Upon error, all entries will be null.
pub extern fn linear_get_dualpol_data(path:*const u8, pathlen:libc::c_int) -> LinearResult{
    // parse the path arguments
    let zippath = unsafe {std::slice::from_raw_parts(path, (pathlen as usize)*std::mem::size_of::<u8>())};
    match interface::linear_get_dualpol_data(zippath, &LinearConfig::default()) {
	Ok((x, co16, k)) => {
	    // TODO check x and co16 are contigious
	    let nss = k.len();
	    let result = LinearResult {
		crosspol:x.as_mut_ptr(),
		copol:co16.as_mut_ptr(),
		k:k.as_mut_ptr(),
		rows:rows as libc::c_int,
		cols:cols as libc::c_int,
		n_subswaths:nss,
	    };

	    // Prevent over-freeing. The calling C-function owns these now.
	    std::mem::forget(x);
	    std::mem::forget(co16);
	    std::mem::forget(m);

	    return result;
	},
	Err(e) => {
	    eprintln!("{}",e);
	    return LinearResult{crosspol:ptr::null::<libc::c_float>() as *mut f64,
				copol:ptr::null::<libc::c_float>() as *mut f64,
				k:ptr::null::<libc::c_float>() as *mut f64,
				rows:0, cols:0, n_subswaths:0};
	}
    }
}

// pub extern fn linear_get_customscale_data(path:*const u8, pathlen:libc::c_int, c_k:*mut f64) -> LinearResult {
    
// }
