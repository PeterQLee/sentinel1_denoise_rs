extern crate denoise_engine;
use denoise_engine::parse::{NoiseField, SwathElem, TimeRowLut, BurstEntry, RawPattern};
use denoise_engine::read_from_archive::{SentinelFormatId,get_data_from_zip_path,  get_id_prefix};
use denoise_engine::apply::{apply_swath_scale, prep_measurement, restore_scale};
use denoise_engine::estimate::*;
use std::fs::File;
use zip::read::{ZipArchive};
use std::io::prelude::*;
use std::io::Cursor;

    
#[test]
/// Test timerowlut extraction against the python implementation.
fn test_timerow(){
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

			    let lut = TimeRowLut::new(&buffer, &swath_bounds, &id);

			    assert!(lut.lut[0].len() == 126 &&
				    lut.lut[1].len() == 84 &&
				    lut.lut[2].len() == 84 &&
				    lut.lut[3].len() == 84 &&
				    lut.lut[4].len() == 63);

			    assert!((lut.lut[0][0].aztime - 147372572.172092).abs() <1e-5);

			    assert!(lut.lut[1][7].row == 500);
			    assert!(lut.lut[1][7].col == 4680);
			    assert!((lut.lut[1][7].aztime - 147372575.174805).abs() <1e-5);
			    assert!(lut.lut[1][7].elevangle != 0.0);
				    
			    
                        }
		    }
		},
		None=>{}
	    }
	},
	Err(E) => {}
    }
}


#[test]
/// Test timerowlut extraction against the python implementation.
fn test_burstderive(){
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

			    let lut = TimeRowLut::new(&buffer, &swath_bounds, &id);

			    let bt = BurstEntry::create_burst_coords(&buffer,
								     &lut,
								     &swath_bounds,
								     &id);

			    assert!(bt[0].len() == 18);
			    assert!(bt[1].len() == 19);
			    assert!(bt[2].len() == 19);
			    assert!(bt[3].len() == 19);
			    assert!(bt[4].len() == 19);

			    println!("{:?}", bt[1][7]);
			    println!("{:?}", bt[1][0]);
			    assert!(bt[1][7].fa == 3564);
			    assert!(bt[1][7].la == 4070);
			    assert!(bt[1][7].fr == 3015);
			    assert!(bt[1][7].lr == 4917);
			    

			    
                        }
		    }
		},
		None=>{}
	    }
	},
	Err(E) => {}
    }
}


#[test]
/// Test timerowlut extraction against the python implementation.
fn test_rawpattern(){
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

			    let lut = RawPattern::new(&buffer, &id);

			 


                        }
		    }
		},
		None=>{}
	    }
	},
	Err(E) => {}
    }
}



#[test]
/// Test timerowlut extraction against the python implementation.
/// Test results: good enough
fn test_interpolation(){
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

			    let lut = TimeRowLut::new(&buffer, &swath_bounds, &id);

			    lut.interp_angle_to_col();

				    
			    
                        }
		    }
		},
		None=>{}
	    }
	},
	Err(E) => {}
    }
}
