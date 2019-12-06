/// Module for reading from sentinel-1 zip archive.
use crate::parse::{NoiseField, SwathElem};
use zip::read::{ZipArchive};
use std::io::prelude::*;
use std::io::Cursor;
use std::fs::File;
use regex::Regex;
use std::string::String;
use tiff::decoder::{Decoder, DecodingResult};
use ndarray::{Array,Array2};

struct SentinelFormatId {
    upper_dateid:String,
    sentid:String,
    polarization:String,
    lower_dateid:String
}

impl SentinelFormatId {
    fn create_crosspol_annotation(&self) -> String {
        let polid;
        match self.polarization.as_str() {
            "h" => {polid = "hv"},
            _ => {polid = "vh"}
        }
        return format!("annotation/s1{}-ew-grd-{}-{}-002.xml", self.sentid, polid, self.lower_dateid);
    }

    fn create_crosspol_noise(&self) -> String {
        let polid;
        match self.polarization.as_str() {
            "h" => {polid = "hv"},
            _ => {polid = "vh"}
        }
        return format!("annotation/calibration/noise-s1{}-ew-grd-{}-{}-002.xml", self.sentid, polid, self.lower_dateid);
    }

    fn create_crosspol_measurement(&self) -> String {
        let polid;
        match self.polarization.as_str() {
            "h" => {polid = "hv"},
            _ => {polid = "vh"}
        }
        return format!("measurement/s1{}-ew-grd-{}-{}-002.tiff", self.sentid, polid, self.lower_dateid);
    }
    
}

/// Get the dateid prefix for the relevant files.
///ex: S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.SAFE-report-20180902T190505.pdf
/// Returns (Uppercase string, lowercase dateid, sentinelid[a or b])
fn get_id_prefix(token:&str) -> Option<SentinelFormatId>{
    // Takes the pdf name as input
    let pdf_re = Regex::new(r"^S1([A|B])_EW_GRD._1SD([H|V])_(.*)\.SAFE.*\.pdf$").unwrap();
    if pdf_re.is_match(token) {
        let mut upper_dateid:String = String::new();
        let mut sentid:String = String::new();
        let mut polarization:String = String::new();
        for t in pdf_re.captures_iter(token) {
            sentid.insert_str(0,&t[0]);
            polarization.insert_str(0,&t[1]);
            upper_dateid.insert_str(0,&t[2]);

        }
        polarization.make_ascii_lowercase();
        sentid.make_ascii_lowercase();
        
        let mut lower_dateid:String = upper_dateid.clone();
        lower_dateid.make_ascii_lowercase();
        lower_dateid = lower_dateid.replace("_","-");

        Some(SentinelFormatId {
            upper_dateid:upper_dateid,
            sentid:sentid,
            polarization:polarization,
            lower_dateid:lower_dateid
        });
    }
    None
}

/// Gets the original 16-bit image and noise field from a given path.
pub fn get_data_from_zip_path(path:&str) -> Option<(Vec<Vec<SwathElem>>, Vec<usize>, NoiseField, Array2<u16>)> {
    let file_h = File::open(path);

    match file_h {
        Ok(f) => {
            let mut ziparch = ZipArchive::new(f).unwrap();

            let mut sentid: Option<SentinelFormatId> = None;
            let mut found_id_flag = false;
            // TODO: find ids and stuff
            for i in 0..ziparch.len() {
                let filename = ziparch.by_index(i).unwrap();
                match (get_id_prefix(filename.name())) {
                    Some(id_result) => {
                        found_id_flag = true;
                        sentid = Some(id_result);
                        break;
                    },
                    None => {}
                }
            }

            match sentid {
                Some(id) => {
                    // generate the matching xml files for the calibration files we want.
                    let crosspol_anno = id.create_crosspol_annotation();
                    let crosspol_noise = id.create_crosspol_noise();
                    let crosspol_measurement = id.create_crosspol_measurement();

                    let mut swath_bounds:Option<Vec<Vec<SwathElem>>> = None;
                    let mut w:Option<Vec<usize>> = None;
                    let mut noisefield:Option<NoiseField> = None;
                    let mut measurement_array:Option<Array2<u16>> = None;
                    
                    for i in 0..ziparch.len() {
                        let mut file = ziparch.by_index(i).unwrap();
                        if file.name() == crosspol_anno {
                            let mut buffer = String::new();
                            let xmldata = file.read_to_string(&mut buffer).unwrap();
                            let k = SwathElem::new(&buffer);
                            swath_bounds = Some(k.0);
                            w = Some(k.1);
                            
                        }

                        else if file.name() == crosspol_noise {
                            let mut buffer = String::new();
                            let xmldata = file.read_to_string(&mut buffer).unwrap();
                            noisefield = Some(NoiseField::new(&buffer, true));
                            
                        }

                        else if file.name() == crosspol_measurement {
                            let mut buffer = Vec::new();
                            let xmldata = file.read_to_end(&mut buffer);
                            let virt_file = Cursor::new(buffer);
                            let mut tiff_file = Decoder::new(virt_file).unwrap();

                            let tiff_dims = tiff_file.dimensions();
                            let (x,y):(u32,u32);
                            match tiff_dims {
                                Ok(t) => {x = t.0; y = t.1;}
                                Err(_e) => {
                                    println!("The tiff file is not encoded properly");
                                    return None;
                                }
                            }
                            
                            
                            let tiff_result = tiff_file.read_image();
                            match tiff_result {
                                Ok(dec_result) => {
                                    if let DecodingResult::U16(dvec) = dec_result {
                                        measurement_array = Some(Array::from_shape_vec((x as usize,y as usize), dvec).unwrap());
                                    }
                                    else {
                                        println!("The tiff file is not encoded properly (Bitdepth)."); return None;
                                    }
                                },
                                Err(_e) => {println!("The tiff file is not encoded properly"); return None;}
                            }
                            
                        }
                    }

                    if swath_bounds.is_some() && w.is_some() && noisefield.is_some() && measurement_array.is_some() {
                        return Some((swath_bounds.unwrap(), w.unwrap(), noisefield.unwrap(), measurement_array.unwrap()));
                    }

                    println!("One or more files not found");
                    return None;
                }
                
                None => {
                    println!("The given zip archive is not formatted as Sentinel-1 archive");
                    return None;
                }
                    
            }
            
        }
        Err(_e) => {
            println!("Cannot open zipfile {}", path);
            return None;
        }
    }

}
