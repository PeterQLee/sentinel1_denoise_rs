//! Interface for applying the methods in specific cases.
use crate::read_from_archive::{get_data_from_zip_path, SentinelArchiveOutput};
use crate::apply::{apply_swath_scale, prep_measurement,  LpApply};
use crate::parse::{SwathElem,  HyperParams, LinearConfig};
use crate::estimate::*;
use ndarray::{Array2, Array1, ArrayView1, arr1};
use crate::est_lp;

use std::sync::Arc;
use crate::prep_lp::*;


/// Applies the lstsquares estimation method to retrieve scaling parameters, k,
/// rescales the noise floor, y, and subtracts it from the image, x.
/// Returns the values back in square intensity units.
/// Parameters:
///
/// * `archpath` - Path to the zip or directory unpacked from the Sentinel-1 zip archive
///     
/// * `linpar` - Linear parameters to feed into the algorithm.
pub fn linear_get_dualpol_data(archpath:&str, linpar:&LinearConfig) 
                               -> Result<(Array2<f64>, Array2<u16>, Array1<f64>), String>
{
    let archout = get_data_from_zip_path(archpath, true)?;
    
    match archout {
        SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, _lpargs) => {
            let y = noisefield.data.view_mut();
            let mut x = prep_measurement(x16.view(), y);
            
            let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
            
            let k = estimate_k_values(x.view(), noisefield.data.view(), &w, &sb, linpar.mu, linpar.gamma, &linpar.lambda, linpar.lambda2);
            
            apply_swath_scale(x.view_mut(), noisefield.data.view(), k.view(), &sb);
            
            return Ok((x, co16, k));
            
        },
        _ => {}
        
    }

    return Err("Error parsing archive".into());

}

/// Applies the linear scaling method using custom user provided scales, k.
/// Returns the values back in square intensity units.
///
/// Parameters:
///
/// * `archpath` -
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// * `k` -
///     One dimensional array in 64-bit float that indicates the linear scaling parameters
///     to apply to each subswath. Length of array must equal five.
///     To apply the standard ESA noise removal, make k = np.array([1.0,1.0,1.0,1.0,1.0])
pub fn linear_get_customscale_data(archpath:&str, k:ArrayView1<f64>)
                                   -> Result<(Array2<f64>, Array2<u16>), String>
{

    let archout = get_data_from_zip_path(archpath, true)?;
    match archout {
        SentinelArchiveOutput::BothPolOutput(swath_bounds, _w, mut noisefield, x16, co16, _lpargs) => {
            let y = noisefield.data.view_mut();
            let mut x = prep_measurement(x16.view(), y);
            
            let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
            
            apply_swath_scale(x.view_mut(), noisefield.data.view(), k, &sb);
                    
            return Ok((x, co16));
            
        },
        _ => {}
    }
    return Err("Error parsing archive".into());
}

/// Returns the original cross pol, co pol, and noise field from the archive.
/// Note that arrays are in linear units.
///
/// Parameters:
///
/// * `archpath` -
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
pub fn linear_get_raw_data( archpath:&str)
                            -> Result<(Array2<u16>, Array2<u16>, Array2<f64>), String> {
    
    let archout = get_data_from_zip_path(archpath, true)?;
    match archout {
        SentinelArchiveOutput::BothPolOutput(_swath_bounds, _w, noisefield, x16, co16, _lpargs) => {
            let y = noisefield.data;
            return Ok((x16, co16, y));
        },
        _ => {}
    }
    return Err("Error parsing archive".into());
}    

/// Applies the linear programming method to restimate a noise floor based
/// on the characteristics of the original image and the 
/// Returns the values back in square intensity units.
///
/// Parameters:
///
/// * `archpath` -
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// * `lstsq_rescale` - 
///     Indicate whether you want to apply the least squares method from
///     linear_get_dualpol_data to get the baseline minimum offset values for
///     the method. Ignored if product type is 
///     true for applying the method
///     false to just use the default ESA noise floor for this.
/// * `linpar` - linear configuration
/// * `hp` - Hyper parameters to feed into the algorithms

pub fn lp_get_dualpol_data(archpath:&str, lstsq_rescale:bool, linpar:&LinearConfig, hp:HyperParams)
                           -> Result<(Arc<TwoDArray>, Array2<u16>, Vec<Vec<est_lp::lin_params>>, Array1<f64>), String>{
    let archout = get_data_from_zip_path(archpath, true)?;

    match archout {
        SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, lpargs) => {

            let mut mv:Array2<f64>;
            {
                let y_ = noisefield.data.view_mut();
                mv = prep_measurement(x16.view(), y_);
            }
            
            let x = mv.clone();
            
            
            let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
                        
            let k = match lpargs.id.sentmode.as_str() {
                "EW" =>  match lstsq_rescale {
                    true => estimate_k_values(x.view(), noisefield.data.view(), &w, &sb, linpar.mu, linpar.gamma, &linpar.lambda, linpar.lambda2),
                    false => arr1(&[1.0,1.0,1.0,1.0,1.0])
                }
                "IW" => arr1(&[1.0,1.0,1.0]),
                _ => { return Err("Invalid sentinel sensor mode".into()); }
            };
            {
                let nf = noisefield;
                apply_swath_scale(mv.view_mut(), nf.data.view(), k.view(), &sb);
            }
            

            let xv = Arc::new(TwoDArray::from_ndarray(x));
            let hyper:Arc<HyperParams> = Arc::new(hp);
            let mino_list:Vec<Vec<f64>>;
            {
                let base_v = Arc::new(TwoDArray::from_ndarray(mv));

                
                        // compute minimum on based on baseline
                
                mino_list = compute_mino_list(
                    base_v,
                    &lpargs.bt,
                    hyper.clone(),
                    &lpargs.id);
            }
                

            //compute weights
            let w_vec = LpApply::compute_weights_for_affine(
                &xv,
                &lpargs.bt,
                hyper.clone(),
                &lpargs.id);

            let params = select_and_estimate_segments(xv.clone(),
                                                      lpargs.mp_dict,
                                                      &lpargs.bt,
                                                      &lpargs.split_indices,
                                                      mino_list,
                                                      hyper.clone(),
                                                      &lpargs.id);
            
                            
            
            let sb = Arc::new(swath_bounds);
            let sb1 = sb.clone();
            // apply power function
            LpApply::apply_lp_noisefield(xv.clone(),
                                         lpargs.az_noise.clone(),
                                         lpargs.eval_mp_dict,
                                         &lpargs.bt,
                                         &lpargs.eval_split_indices,
                                         sb,
                                         &params,
                                         hyper.clone(),
                                         &lpargs.id);
            // apply affine offsets.
            
            LpApply::apply_affine(xv.clone(),
                                  w_vec,
                                  &lpargs.bt,
                                  sb1,
                                  hyper.clone(),
                                  &lpargs.id);
            
            
            return Ok((xv, co16, params, k));
        },
        _ => {}
    }
    return Err("Error parsing archive".into());
}

///  Applies the power function noise floor obtained from linear programming,
///   with parameters given by the user.
/// Returns the values back in square intensity units.
/// 
/// Parameters:
///
/// * `archpath` -
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// * `m` - 
///     List of slope parameters
/// * `b` - 
///     List of intercept parameters.
/// * `num_subswaths` - Number of subswaths in image
/// * `hp` - Hyper parameters to feed into the algorithms
pub fn lp_get_customscale_data(archpath:&str, m:&[f64],
                               b:&[f64], num_subswaths:usize,
                               hp:HyperParams) -> Result<(Arc<TwoDArray>, Array2<u16>), String> {
    let archout = get_data_from_zip_path(archpath, true)?;
    match archout {
        SentinelArchiveOutput::BothPolOutput(swath_bounds, _w, mut noisefield, x16, co16, lpargs) => {
            let mv:Array2<f64>;
            {
                let y_ = noisefield.data.view_mut();
                mv = prep_measurement(x16.view(), y_);
            }
	    
            let x = mv.clone();
	    let xv = Arc::new(TwoDArray::from_ndarray(x));


	    let hyper:Arc<HyperParams> = Arc::new(hp);

	    // construct params
	    let sb = Arc::new(swath_bounds);
	    let params = params_from_vecs(m, b, num_subswaths);
	    
	    //compute weights
	    let w_vec = LpApply::compute_weights_for_affine(
		&xv,
		&lpargs.bt,
		hyper.clone(),
		&lpargs.id);
	    
	    LpApply::apply_lp_noisefield(xv.clone(),
					 lpargs.az_noise.clone(),
					 lpargs.eval_mp_dict,
					 &lpargs.bt,
					 &lpargs.eval_split_indices,
					 sb.clone(),
					 &params,
					 hyper.clone(),
					 &lpargs.id);
   	    
	    LpApply::apply_affine(xv.clone(),
				  w_vec,
				  &lpargs.bt,
				  sb,
				  hyper.clone(),
				  &lpargs.id);
	    return Ok((xv, co16));

	},
	_ => {}
    }
    return Err("Error parsing archive".into());
}


fn params_from_vecs(m:&[f64], b:&[f64], num_subswaths:usize) -> Vec<Vec<est_lp::lin_params>> {
    match num_subswaths {
	5 => { //EW mode
	    let num_p:&[usize] = &[4,2,2,2,2];
	    let cumu:&[usize] = &[0,4,6,8,10];
	    num_p.iter().zip(cumu.iter())
		.map(|x| m[*x.1..*x.1+*x.0].iter()
		     .zip(b[*x.1..*x.1+*x.0].iter())
		     .map(|y| est_lp::lin_params{m:*y.0, b:*y.1})
		     .collect()).collect()
	},
	3 =>  { //IW mode
	    let num_p:&[usize] = &[2,2,2];
	    let cumu:&[usize] = &[0,2,4];
	    num_p.iter().zip(cumu.iter())
		.map(|x| m[*x.1..*x.1+*x.0].iter()
		     .zip(b[*x.1..*x.1+*x.0].iter())
		     .map(|y| est_lp::lin_params{m:*y.0, b:*y.1})
		     .collect()).collect()
	},
	_ => {panic!("Invalid number of subswaths")}
    }
}
