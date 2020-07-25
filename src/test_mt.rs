extern crate denoise_engine;
use denoise_engine::parse::{NoiseField, SwathElem, TimeRowLut, BurstEntry, RawPattern, HyperParams};
use denoise_engine::prep_lp::*;
use denoise_engine::read_from_archive::{SentinelFormatId,get_data_from_zip_path,  get_id_prefix, SentinelArchiveOutput};
use denoise_engine::apply::{apply_swath_scale, prep_measurement, restore_scale, LpApply};
use denoise_engine::estimate::*;
use std::fs::File;
use zip::read::{ZipArchive};
use std::io::prelude::*;
use std::io::Cursor;
use std::sync::Arc;
use ndarray::{Array1, Array2, arr1};
use std::time::{Duration, Instant};
fn main() {
    //let path = "/mnt/D2/Data/Sentinel/ant2_test//S1B_EW_GRDM_1SDH_20200112T183748_20200112T183848_019787_02569C_3399.zip";
    let path = "/mnt/D2/Data/Sentinel/ant2_test//S1A_EW_GRDM_1SDH_20191224T193600_20191224T193700_030494_037DD6_1F74.zip";
	//let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    
    let file_h = File::open(path);

    let zipval = get_data_from_zip_path(path, true);
    let archout = zipval.unwrap();
    let mut x:Option<Array2<f64>> = None;
    let mut y:Option<Array2<f64>> = None;
    let mut base:Option<Array2<f64>> = None;
    match archout {
	SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16) => {
            let mut x_:Option<Array2<f64>> = None;

            
	    let mut y_ = noisefield.data.view_mut();
	    let mut mv = prep_measurement(x16.view(), y_);
	    x_ = Some(mv.clone());

	    y = Some(noisefield.data);
            
            x = x_;
	    println!("stmv={}", mv[(0,100)]);
	    println!("stmv={}", mv[(100,0)]);
	    println!("mv={}", mv[(500,500)]);

	    let sb:Vec<&[SwathElem]> = swath_bounds.iter().map(|a| a.as_slice()).collect();
	    let k = arr1(&[1.0,1.0,1.0,1.0,1.0]);
	    apply_swath_scale(mv.view_mut(), y.unwrap().view(), k.view(), &sb);
	    base = Some(mv);
	    
	},
	_ => {}
    }
    
    match file_h {
        Ok(f) => {
            let mut ziparch = ZipArchive::new(f).unwrap();

            let mut sentid: Option<SentinelFormatId> = None;
            // TODO: find ids and stuff
            for i in 0..ziparch.len() {
                let filename = ziparch.by_index(i).unwrap();

                match get_id_prefix(filename.name()) {
                    Some(id_result) => {
                        sentid = Some(id_result);
                        break;
                    },
                    None => {}
                }
            }
	    match sentid {
                Some(id) => {
		    let xk = x.unwrap();
		    let original = xk.clone();
		    
		    let crosspol_anno = id.create_crosspol_annotation();
		    let crosspol_noise = id.create_crosspol_noise();
		    let mut noisefield:Option<NoiseField> = None;
		    for i in 0..ziparch.len() {
                        let mut file = ziparch.by_index(i).unwrap();
			if file.name() == crosspol_noise {
                            let mut buffer = String::new();
                            let xmldata = file.read_to_string(&mut buffer).unwrap();
                            noisefield = Some(NoiseField::compute_azimuth_field(&buffer, (xk.shape()[0], xk.shape()[1])));
			}
		    }
		    
		    let rewrap = Instant::now();
		    let xv = Arc::new(TwoDArray::from_ndarray(xk));
		    //let yv = TwoDArray::from_ndarray(y.unwrap());
		    let yv = Arc::new(TwoDArray::from_ndarray(noisefield.unwrap().data));
		    let base_v = Arc::new(TwoDArray::from_ndarray(base.unwrap()));

		    println!("rewraptime = {}", rewrap.elapsed().as_secs_f64());
		    for i in 0..ziparch.len() {
                        let mut file = ziparch.by_index(i).unwrap();

			if file.name() == crosspol_anno {
			    let mut buffer = String::new();
                            let xmldata = file.read_to_string(&mut buffer).unwrap();
			    let k = SwathElem::new(&buffer, &id);
			    let swath_bounds = k.0;

			    let rawpatt = RawPattern::new(&buffer, &id);
			    let lut = TimeRowLut::new(&buffer, &swath_bounds, &id);
			    let bt = BurstEntry::create_burst_coords(&buffer,
								     &lut,
								     &swath_bounds,
								     &id);
			    let hyper:Arc<HyperParams> = Arc::new(HyperParams::default());

			    let min_time = Instant::now();
			    /*
			    println!("first_burst = {} {} {} {}",
				     bt[0][0].fa,
				     bt[0][0].la,
				     bt[0][0].fr,
				     bt[0][0].lr);
			   
			    println!("lut {} {} {} {}", lut.lut[0][0].row, lut.lut[0][0].col, lut.lut[0][0].elevangle, lut.lut[0][0].aztime);
			    println!("lut {} {} {} {}", lut.lut[0][1].row, lut.lut[0][1].col, lut.lut[0][1].elevangle, lut.lut[0][1].aztime);*/

			    /*
			    println!("start={}", base_v[(0,0)]);
			    println!("start={}", base_v[(0,100)]);
			    println!("start={}", base_v[(100,0)]);
			    println!("{}", base_v[(500,500)]);*/
			    let mut mino_list:Vec<Vec<f64>> = compute_mino_list(
				base_v.clone(),
				//xv.clone(),
				&bt,
				hyper.clone(),
				&id);
			    
			    
			    println!("mino_time = {}", min_time.elapsed().as_secs_f64());
			    println!("mino_values={:?}", mino_list);
			    
			    panic!("");

			    let (mp_dict, split_indices) = get_interpolation_pattern(&buffer,
										     &bt,
										     &lut,
										     &rawpatt,
										     &id,
										     false
			    );



			    let now = Instant::now();
			    println!("STarting estimation");
			    
			    let params = select_and_estimate_segments(xv.clone(),
								      mp_dict,
								      &bt,
								      &split_indices,
								      mino_list,
								      hyper.clone(),
								      &id);
			    println!("esttime = {}", now.elapsed().as_secs_f64());
			    
			    for i in params.iter() {
				println!("Size = {}", i.len());
				for j in i.iter() {
				    println!("m={} b={}", j.m, j.b);
				}
			    }

			    let (r_mp_dict, r_split_indices) = get_interpolation_pattern(&buffer,
										     &bt,
										     &lut,
										     &rawpatt,
										     &id,
											 true
			    );

			    let applytime = Instant::now();
			    let sb = Arc::new(swath_bounds);
			    let sb1 = sb.clone();
			    LpApply::apply_lp_noisefield(xv.clone(),
							 yv.clone(),
							 r_mp_dict,
							 &bt,
							 &r_split_indices,
							 sb,
							 &params,
							 hyper.clone(),
							 &id);
			    println!("apply = {}", applytime.elapsed().as_secs_f64());

			    LpApply::apply_affine(xv.clone(),
						  TwoDArray::from_ndarray(original.clone()),
						  sb1,
						  hyper.clone(),
						  &id);
						  

			    

                        }
		    }
		},
		None=>{}
	    }
	},
	Err(E) => {}
    }

}
