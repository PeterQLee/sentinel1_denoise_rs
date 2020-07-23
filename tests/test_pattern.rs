extern crate denoise_engine;
use denoise_engine::parse::{NoiseField, SwathElem, TimeRowLut, BurstEntry, RawPattern};
use denoise_engine::prep_lp::*;
use denoise_engine::read_from_archive::{SentinelFormatId,get_data_from_zip_path,  get_id_prefix, SentinelArchiveOutput};
use denoise_engine::apply::{apply_swath_scale, prep_measurement, restore_scale};
use denoise_engine::estimate::*;
use std::fs::File;
use zip::read::{ZipArchive};
use std::io::prelude::*;
use std::io::Cursor;
use std::sync::Arc;


#[test]
/// Test timerowlut extraction against the python implementation.
fn test_pattern_funcs(){
    let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    
    let file_h = File::open(path);
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
		    let crosspol_anno = id.create_crosspol_annotation();
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

			    assert!((mp_dict[0][5](&[600])[0] - 1.78497e-05).abs() <1e-7);
			    // println!("v={}",mp_dict[1][1](&[3600])[0]);
			    // println!("v={}",mp_dict[2][8](&[6800])[0]);
			    // println!("v={}",mp_dict[3][9](&[8000])[0]);
			    // println!("v={}",mp_dict[4][2](&[9000])[0]);
			    
			    assert!((mp_dict[1][1](&[3600])[0] - 5.0644e-05).abs() < 1e-7);
			    assert!((mp_dict[2][8](&[6800])[0] - 5.0173e-05).abs() < 1e-7);
			    assert!((mp_dict[3][9](&[8000])[0] - 0.00016382).abs() < 1e-7);
			    assert!((mp_dict[4][2](&[9000])[0] - 0.00010149).abs() < 1e-7);

			    // match split_indices {
			    // 	MidPoint::Est(m) => println!("{:?}\n\n", m),
			    // 	MidPoint::Test(m) => {},
			    // }

			    let (r_mp_dict, r_split_indices) = get_interpolation_pattern(&buffer,
										     &bt,
										     &lut,
										     &rawpatt,
										     &id,
										     true
			    );
			    // match r_split_indices {
			    // 	MidPoint::Est(m) => {},
			    // 	MidPoint::Test(m) => {println!("{:?}", m)},
			    // }


                        }
		    }
		},
		None=>{}
	    }
	},
	Err(E) => {}
    }
}



use ndarray::{Array1, Array2, arr1};

#[test]
fn test_segment_gathering() {
    let path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
    
    let file_h = File::open(path);

    let zipval = get_data_from_zip_path(path, true);
    let archout = zipval.unwrap();
    let mut x:Option<Array2<f64>> = None;
    match archout {
	SentinelArchiveOutput::BothPolOutput(swath_bounds, w, mut noisefield, x16, co16) => {
            let mut x_:Option<Array2<f64>> = None;
            {
		let mut y = noisefield.data.view_mut();
		x_ = Some(prep_measurement(x16.view(), y));
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
		    let xv = Arc::new(TwoDArray::from_ndarray(xk));
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
			    
			    let hyper = HyperParams{
				box_l:51,
				add_pad:0
			    };

			    
			    let params = select_and_estimate_segments(xv.clone(),
								      mp_dict,
								      &bt,
								      &split_indices,
								      o_list,
								      Arc::new(hyper),
								      &id);
			    for i in params.iter() {
				for j in i.iter() {
				    println!("m={} b={}", j.m, j.b);
				}
			    }


                        }
		    }
		},
		None=>{}
	    }
	},
	Err(E) => {}
    }

}
