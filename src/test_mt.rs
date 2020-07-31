extern crate denoise_engine;
use denoise_engine::parse::{NoiseField, SwathElem, TimeRowLut, BurstEntry, RawPattern, HyperParams};
use denoise_engine::prep_lp::*;
use denoise_engine::read_from_archive::{SentinelFormatId,get_data_from_zip_path,  get_id_prefix, SentinelArchiveOutput};
use denoise_engine::apply::{apply_swath_scale, prep_measurement, LpApply};
use hdf5;
use hdf5::types::Array;

use std::fs::File;
use zip::read::{ZipArchive};
use std::io::prelude::*;

use std::sync::Arc;
use ndarray::{Array2, arr1, Ix2};
use std::time::{Instant};
fn main() {
    //let path = "/mnt/D2/Data/Sentinel/ant2_test//S1B_EW_GRDM_1SDH_20200112T183748_20200112T183848_019787_02569C_3399.zip";
    //let path = "/mnt/D2/Data/Sentinel/ant2_test//S1A_EW_GRDM_1SDH_20191224T193600_20191224T193700_030494_037DD6_1F74.zip";
    //let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    let path = "/mnt/D2/Data/Sentinel/ant2_test//S1B_EW_GRDM_1SDH_20191218T193518_20191218T193618_019423_024B12_5B75.zip";
    
    let file_h = File::open(path);

    let zipval = get_data_from_zip_path(path, true);
    let archout = match zipval {
	Ok(ac) => ac,
	Err(e) => panic!(e)};
    /*
    let mut x:Option<Array2<f64>> = None;
    let mut y:Option<Array2<f64>> = None;
    let mut base:Option<Array2<f64>> = None;
     */
    match archout {
	SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16, lpargs) => {

	    let mut y_ = noisefield.data.view_mut();
	    let mut mv = prep_measurement(x16.view(), y_);
	    let x = mv.clone();

	    let y = noisefield.data;
            
	    println!("stmv={}", mv[(0,100)]);
	    println!("stmv={}", mv[(100,0)]);
	    println!("mv={}", mv[(500,500)]);

	    let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
	    let k = arr1(&[1.0,1.0,1.0,1.0,1.0]);
	    apply_swath_scale(mv.view_mut(), y.view(), k.view(), &sb);
	    
	    let base_v = Arc::new(TwoDArray::from_ndarray(mv));
	    let xv = Arc::new(TwoDArray::from_ndarray(x));
	    let yv = Arc::new(TwoDArray::from_ndarray(y));

	    
	    let hyper:Arc<HyperParams> = Arc::new(HyperParams::default());

	    // cmpute minimum on based on ESA
	    let mut mino_list:Vec<Vec<f64>> = compute_mino_list(
		base_v.clone(),
		&lpargs.bt,
		hyper.clone(),
		&lpargs.id);
	    
	    //compute weights
	     let w_vec = LpApply::compute_weights_for_affine(
		 &xv,
		 &lpargs.bt,
		 hyper.clone(),
		 &lpargs.id);

	    println!("w_vec={:?}",w_vec);
	    
	    // computer slope/intercept params
	    let params = select_and_estimate_segments(xv.clone(),
			    			      lpargs.mp_dict,
			    			      &lpargs.bt,
			    			      &lpargs.split_indices,
			    			      mino_list,
			    			      hyper.clone(),
			    			      &lpargs.id);
	    for i in params.iter() {
		println!("Size = {}", i.len());
		for j in i.iter() {
		    println!("m={} b={}", j.m, j.b);
		}
	    }
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
						  
	    println!("Done");
	    // convert x back to ndarray

	    // no more references should exist at this point
	    let xout = Arc::try_unwrap(xv).expect("Could not unwrap");
	    let xview = xout.to_ndarray();

	    let file = hdf5::File::create("/mnt/D2/test_mt841.h5").unwrap();
	    let fout = file.new_dataset::<f64>().create("ab", [xout.rows, xout.cols]).unwrap();

	    fout.write(xview).expect("lol");
	    
	    
	    
	},
	_ => {}
    }
    


}
