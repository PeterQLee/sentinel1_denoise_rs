use crate::read_from_archive::{get_data_from_zip_path, SentinelArchiveOutput};
use crate::apply::{apply_swath_scale, prep_measurement, restore_scale, convert_to_f64_f32, convert_to_u16_f32, LpApply};
use crate::parse::{SwathElem,  HyperParams, LinearConfig};
use crate::estimate::*;
use ndarray::{Array2, Array1, ArrayView1, arr1};
use std::ptr;
use numpy::{PyArray, PyArray1, PyArray2};
use crate::est_lp;

use std::sync::Arc;
use crate::prep_lp::*;

pub fn linear_get_dualpol_data(zippath:&str, linpar:&LinearConfig) 
			       -> Result<(Array2<f64>, Array2<u16>, Array1<f64>), String>
{
    let archout = get_data_from_zip_path(zippath, true)?;
    
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


pub fn linear_get_customscale_data(zippath:&str, k:ArrayView1<f64>, linpar:&LinearConfig)
				   -> Result<(Array2<f64>, Array2<u16>), String>
{

    let archout = get_data_from_zip_path(zippath, true)?;
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

pub fn linear_get_raw_data( zippath:&str)
			    -> Result<(Array2<u16>, Array2<u16>, Array2<f64>), String> {
    
    let archout = get_data_from_zip_path(zippath, true)?;
    match archout {
        SentinelArchiveOutput::BothPolOutput(_swath_bounds, _w, noisefield, x16, co16, _lpargs) => {
            let y = noisefield.data;
	    return Ok((x16, co16, y));
	},
	_ => {}
    }
    return Err("Error parsing archive".into());
}    


pub fn lp_get_dualpol_data(zippath:&str, lstsq_rescale:bool, linpar:&LinearConfig)
			   -> Result<(Arc<TwoDArray>, Array2<u16>, Vec<Vec<est_lp::lin_params>>), String>{
    let archout = get_data_from_zip_path(zippath, true)?;

    match archout {
        SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, lpargs) => {

            let mut mv:Array2<f64>;
	    {
		let y_ = noisefield.data.view_mut();
		mv = prep_measurement(x16.view(), y_);
	    }
	    
            let mut x = mv.clone();
	    
            
            let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
                        
            let k = match lpargs.id.sentmode.as_str() {
		"EW" =>  match lstsq_rescale {
		    true => estimate_k_values(x.view(), noisefield.data.view(), &w, &sb, linpar.mu, linpar.gamma, &linpar.lambda, linpar.lambda2),
		    false => arr1(&[1.0,1.0,1.0,1.0,1.0])
		}
		"IW" => arr1(&[1.0,1.0,1.0]),
		_ => { return Err("Invalid sentinel sensor mode".into()); }
	    };

            // TODO: need to test this for IW mode
            apply_swath_scale(mv.view_mut(), noisefield.data.view(), k.view(), &sb);
	    

	    let xv = Arc::new(TwoDArray::from_ndarray(x));
	    let hyper:Arc<HyperParams> = Arc::new(HyperParams::default());
	    let mino_list:Vec<Vec<f64>>;
	    {
		let base_v = Arc::new(TwoDArray::from_ndarray(mv));

		
			// compute minimum on based on baseline
		
		mino_list = compute_mino_list(
		    base_v.clone(),
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
	    println!("prestimate segments");
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
	    
	    
	    return Ok((xv, co16, params));
        },
        _ => {}
    }
    return Err("Error parsing archive".into());
}

pub fn lp_get_customscale_data() {
}
