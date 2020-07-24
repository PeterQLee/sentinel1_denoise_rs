/// Module for reading from sentinel-1 zip archive.
use crate::parse::{NoiseField, SwathElem};
use zip::read::{ZipArchive};
use std::io::prelude::*;
use std::io::Cursor;
use std::fs::File;
use regex::Regex;
use std::string::String;
use tiff::decoder::{Decoder, DecodingResult, Limits};
use ndarray::{Array,Array2};


#[derive(Debug)]
pub struct SentinelFormatId {
    upper_dateid:String,
    sentid:String,
    pub sentmode:String,
    lower_sentmode:String,
    pub quality:String,
    polarization:String,
    lower_dateid:String
}

pub enum SentinelArchiveOutput{
    CrossPolOutput (
        Vec<Vec<SwathElem>>,
        Vec<usize>,
        NoiseField,
        Array2<u16>
    ),
    BothPolOutput (
        Vec<Vec<SwathElem>>,
        Vec<usize>,
        NoiseField,
        Array2<u16>,
        Array2<u16>
    )
}

 

impl SentinelFormatId {
    pub fn create_crosspol_annotation(&self) -> String {
        let polid;
        let upemit;
        
        match self.polarization.as_str() {
            "h" => {polid = "hv"; upemit = "H"},
            _ => {polid = "vh"; upemit = "V"}
        }

        let up_sentid;
        match self.sentid.as_str() {
            "a" => up_sentid = "A",
            _ =>  up_sentid = "B"
        }
        
        return format!("S1{up_sentid}_{sentmode}_GRD{quality}_1SD{upemit}_{upper_dateid}.SAFE/annotation/s1{sentid}-{lower_sentmode}-grd-{polid}-{lower_dateid}-002.xml", up_sentid = up_sentid, sentmode = self.sentmode, lower_sentmode = self.lower_sentmode, quality = self.quality, upemit = upemit, upper_dateid = self.upper_dateid, sentid = self.sentid, polid = polid, lower_dateid = self.lower_dateid);
    }

    pub fn create_crosspol_noise(&self) -> String {
        let polid;
        let upemit;
        
        match self.polarization.as_str() {
            "h" => {polid = "hv"; upemit = "H"},
            _ => {polid = "vh"; upemit = "V"}
        }

        let up_sentid;
        match self.sentid.as_str() {
            "a" => up_sentid = "A",
            _ =>  up_sentid = "B"
        }
        return format!("S1{up_sentid}_{sentmode}_GRD{quality}_1SD{upemit}_{upper_dateid}.SAFE/annotation/calibration/noise-s1{sentid}-{lower_sentmode}-grd-{polid}-{lower_dateid}-002.xml", up_sentid = up_sentid, sentmode = self.sentmode, lower_sentmode = self.lower_sentmode, quality = self.quality, upemit = upemit, upper_dateid = self.upper_dateid, sentid = self.sentid, polid = polid, lower_dateid = self.lower_dateid);
    }

    fn create_crosspol_measurement(&self) -> String {
        let polid;
        let upemit;
        
        match self.polarization.as_str() {
            "h" => {polid = "hv"; upemit = "H"},
            _ => {polid = "vh"; upemit = "V"}
        }

        let up_sentid;
        match self.sentid.as_str() {
            "a" => up_sentid = "A",
            _ =>  up_sentid = "B"
        }
        return format!("S1{up_sentid}_{sentmode}_GRD{quality}_1SD{upemit}_{upper_dateid}.SAFE/measurement/s1{sentid}-{lower_sentmode}-grd-{polid}-{lower_dateid}-002.tiff", up_sentid = up_sentid, quality = self.quality, sentmode = self.sentmode, lower_sentmode = self.lower_sentmode, upemit = upemit, upper_dateid = self.upper_dateid, sentid = self.sentid, polid = polid, lower_dateid = self.lower_dateid);
    }


    fn create_copol_measurement(&self) -> String {
        let polid;
        let upemit;
        
        match self.polarization.as_str() {
            "h" => {polid = "hh"; upemit = "H"},
            _ => {polid = "vv"; upemit = "V"}
        }

        let up_sentid;
        match self.sentid.as_str() {
            "a" => up_sentid = "A",
            _ =>  up_sentid = "B"
        }
        return format!("S1{up_sentid}_{sentmode}_GRD{quality}_1SD{upemit}_{upper_dateid}.SAFE/measurement/s1{sentid}-{lower_sentmode}-grd-{polid}-{lower_dateid}-001.tiff", up_sentid = up_sentid, quality = self.quality, sentmode = self.sentmode, lower_sentmode = self.lower_sentmode, upemit = upemit, upper_dateid = self.upper_dateid, sentid = self.sentid, polid = polid, lower_dateid = self.lower_dateid);
    }

}

/// Get the dateid prefix for the relevant files.
///ex: S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.SAFE-report-20180902T190505.pdf
/// Returns (Uppercase string, lowercase dateid, sentinelid[a or b])
pub fn get_id_prefix(token:&str) -> Option<SentinelFormatId>{
    // Takes the pdf name as inputg
    let pdf_re = Regex::new(r"^S1[A|B]_[IE]W_GRD._1SD[H|V]_.*\.SAFE/S1([A|B])_([EI]W)_GRD(.)_1SD([H|V])_(.*)\.SAFE.*\.pdf$").unwrap();
    if pdf_re.is_match(token) {
        let mut upper_dateid:String = String::new();
	let mut sentmode:String = String::new();
        let mut sentid:String = String::new();
        let mut polarization:String = String::new();
        let mut quality:String = String::new();
        for t in pdf_re.captures_iter(token) {
            sentid.insert_str(0,&t[1]);
	    sentmode.insert_str(0, &t[2]);
            quality.insert_str(0, &t[3]);
            polarization.insert_str(0,&t[4]);
            upper_dateid.insert_str(0,&t[5]);

        }
        polarization.make_ascii_lowercase();
        sentid.make_ascii_lowercase();
        
        let mut lower_dateid:String = upper_dateid.clone();
        lower_dateid.make_ascii_lowercase();
        lower_dateid = lower_dateid.replace("_","-");

	let mut lower_sentmode:String = sentmode.clone();
	lower_sentmode.make_ascii_lowercase();
	
        let spl = lower_dateid.rfind('-').unwrap();
        let _ = lower_dateid.split_off(spl);

        return Some(SentinelFormatId {
            upper_dateid:upper_dateid,
            sentid:sentid,
	    sentmode:sentmode,
	    lower_sentmode:lower_sentmode,
            quality:quality,
            polarization:polarization,
            lower_dateid:lower_dateid
        });
    }
    else {
        return None;
    }
}

/// Gets the original 16-bit image and noise field from a given path.
pub fn get_data_from_zip_path(path:&str, bothpol_flag:bool) -> Option<SentinelArchiveOutput> {
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
                    // generate the matching xml files for the calibration files we want.

                    let crosspol_anno = id.create_crosspol_annotation();
                    let crosspol_noise = id.create_crosspol_noise();
                    let crosspol_measurement = id.create_crosspol_measurement();
                    let copol_measurement = id.create_copol_measurement();


                    let mut swath_bounds:Option<Vec<Vec<SwathElem>>> = None;
                    let mut w:Option<Vec<usize>> = None;
                    let mut noisefield:Option<NoiseField> = None;
                    let mut measurement_array:Option<Array2<u16>> = None;
                    let mut copol_array:Option<Array2<u16>> = None;

                    let mut shape_o:Option<(usize, usize)> = None;

                    // Measurement files first, so we can get the required shape info
                    for i in 0..ziparch.len() {
                        let mut file = ziparch.by_index(i).unwrap();
                        // Get crosspol array from tiff file.
                        if file.name() == crosspol_measurement {
                            let mut buffer = Vec::new();
                            let xmldata = file.read_to_end(&mut buffer);
                            let virt_file = Cursor::new(buffer);
                            let mut tiff_file = Decoder::new(virt_file).unwrap();
                            let mut limits = Limits::default();
                            limits.decoding_buffer_size = 1_usize<<40;
                            limits.ifd_value_size = 1_usize<<40;
                            tiff_file = tiff_file.with_limits(limits);


                            let tiff_dims = tiff_file.dimensions();
                            let (x,y):(u32,u32);
                            match tiff_dims {
                                Ok(t) => {x = t.0; y = t.1;}
                                Err(e) => {
                                    println!("The tiff file is not encoded properly (dimensions)");
                                    println!("{}", e);
                                    return None;
                                }
                            }
                            shape_o = Some((y as usize, x as usize));
                            
                            let tiff_result = tiff_file.read_image();
                            match tiff_result {
                                Ok(dec_result) => {
                                    if let DecodingResult::U16(dvec) = dec_result {
                                        measurement_array = Some(Array::from_shape_vec((y as usize,x as usize), dvec).unwrap());
                                    }
                                    else {
                                        println!("The tiff file is not encoded properly (Bitdepth)."); return None;
                                    }
                                },
                                Err(e) => {
                                    println!("The tiff file is not encoded properly (readimg)");
                                    println!("{}", e);
                                    return None;
                                }
                            }
                        }

                        // Get crosspol array from tiff file.
                        if bothpol_flag && file.name() == copol_measurement {
                            let mut buffer = Vec::new();
                            let xmldata = file.read_to_end(&mut buffer);
                            let virt_file = Cursor::new(buffer);
                            let mut tiff_file = Decoder::new(virt_file).unwrap();
                            let mut limits = Limits::default();
                            limits.decoding_buffer_size = 1_usize<<40;
                            limits.ifd_value_size = 1_usize<<40;
                            tiff_file = tiff_file.with_limits(limits);

                            let tiff_dims = tiff_file.dimensions();
                            let (x,y):(u32,u32);
                            match tiff_dims {
                                Ok(t) => {x = t.0; y = t.1;}
                                Err(e) => {
                                    println!("The tiff file is not encoded properly (dimensions)");
                                    println!("{}", e);
                                    return None;
                                }
                            }
                            
                            
                            let tiff_result = tiff_file.read_image();
                            match tiff_result {
                                Ok(dec_result) => {
                                    if let DecodingResult::U16(dvec) = dec_result {
                                        copol_array = Some(Array::from_shape_vec((y as usize,x as usize), dvec).unwrap());
                                    }
                                    else {
                                        println!("The tiff file is not encoded properly (Bitdepth)."); return None;
                                    }
                                },
                                Err(e) => {
                                    println!("The tiff file is not encoded properly");
                                    println!("{}", e);
                                    return None;
                                }
                            }
                        }
                    }

                    let shape = shape_o.unwrap();
                    
                    for i in 0..ziparch.len() {
                        let mut file = ziparch.by_index(i).unwrap();


                        // Get swath bounds and period from anno file
                        if file.name() == crosspol_anno {
                            let mut buffer = String::new();
                            let xmldata = file.read_to_string(&mut buffer).unwrap();
                            let k = SwathElem::new(&buffer, &id);
                            swath_bounds = Some(k.0);
                            w = Some(k.1);
                            
                        }

                        // Get noise from noise calibration file.
                        else if file.name() == crosspol_noise {
                            let mut buffer = String::new();
                            let xmldata = file.read_to_string(&mut buffer).unwrap();
                            noisefield = Some(NoiseField::new(&buffer, shape, true));
                            
                        }

                    }

                    // Return the crosspol result only.
                    if !bothpol_flag && swath_bounds.is_some() && w.is_some() && noisefield.is_some() && measurement_array.is_some() {
                        return Some(
                            SentinelArchiveOutput::CrossPolOutput (
                                swath_bounds.unwrap(),
                                w.unwrap(),
                                noisefield.unwrap(),
                                measurement_array.unwrap()
                            )
                        );
                        //return Some((swath_bounds.unwrap(), w.unwrap(), noisefield.unwrap(), measurement_array.unwrap()));
                    }


                    // Both crosspol and co polarization
                    else if bothpol_flag && swath_bounds.is_some() && w.is_some() && noisefield.is_some() && measurement_array.is_some() && copol_array.is_some() {
                        return Some(
                            SentinelArchiveOutput::BothPolOutput (
                                swath_bounds.unwrap(),
                                w.unwrap(),
                                noisefield.unwrap(),
                                measurement_array.unwrap(),
                                copol_array.unwrap()
                            )
                        );

                    
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


