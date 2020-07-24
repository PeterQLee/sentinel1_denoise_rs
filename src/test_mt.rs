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
    let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    
    let file_h = File::open(path);

    let zipval = get_data_from_zip_path(path, true);
    let archout = zipval.unwrap();
    let mut x:Option<Array2<f64>> = None;
    let mut y:Option<Array2<f64>> = None;
    match archout {
	SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16) => {
            let mut x_:Option<Array2<f64>> = None;
            {
		let mut y_ = noisefield.data.view_mut();
		x_ = Some(prep_measurement(x16.view(), y_));
		y = Some(noisefield.data);
            }
            x = x_;
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
		    let yv = TwoDArray::from_ndarray(noisefield.unwrap().data);

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
			    
			    let (mp_dict, split_indices) = get_interpolation_pattern(&buffer,
										     &bt,
										     &lut,
										     &rawpatt,
										     &id,
										     false
			    );

			    let mut o_list:Vec<Vec<f64>> = Vec::new();
			    for s in 0..mp_dict.len() {
				o_list.push(Vec::new());
				for t in 0..mp_dict[s].len() {
				    o_list[s].push(0.0);
				}
			    }
			    
			    let hyper:Arc<HyperParams> = Arc::new(HyperParams{
				box_l:51,
				add_pad:0
			    });

			    let now = Instant::now();
			    println!("STarting estimation");
			    
			    let params = select_and_estimate_segments(xv.clone(),
								      mp_dict,
								      &bt,
								      &split_indices,
								      o_list,
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
			    LpApply::apply_lp_noisefield(xv.clone(),
							 &yv,
							 r_mp_dict,
							 &bt,
							 &r_split_indices,
							 &swath_bounds,
							 &params,
							 hyper.clone(),
							 &id);
			    println!("apply = {}", applytime.elapsed().as_secs_f64());

                        }
		    }
		},
		None=>{}
	    }
	},
	Err(E) => {}
    }

}
