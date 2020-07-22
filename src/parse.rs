
use quick_xml::Reader;
use quick_xml::events::Event;
use std::str;
use ndarray::{Array, Array1, Array2, ArrayBase, Axis, ArrayViewMut1, ArrayView1, ArrayView2, Slice};
use ndarray::Zip;
use ndarray::prelude::*;
//use ndarray::parallel::prelude::*;
//use ndarray_parallel::prelude::*;
use rayon::prelude::*;
//use itertools::Itertools;
use chrono::prelude::*;
use crate::read_from_archive::SentinelFormatId;

use std::collections::HashSet;
// Linear algebra
use lapack::*;
use blas::*;


pub struct NoiseField {
    pub data:Array2<f64>
}

struct NoiseRangeEntry {
    line: usize,
    pixels: Vec<usize>,
    values: Vec<f64>
}


pub struct SwathElem {
    pub fa:usize,
    pub la:usize,
    pub fr:usize,
    pub lr:usize
}

macro_rules! min_row {
    ($sw:expr) => {
	$sw.iter().fold(99999, |acc,x| {
	    if x.fa < acc {return x.fa;}
	    acc});
    };
}

macro_rules! max_row {
    ($sw:expr) => {
	$sw.iter().fold(0, |acc,x| {
	    if x.la > acc {return x.la;}
	    acc});
    };
}



struct NoiseAzimuthEntry {
    lines: Vec<usize>,
    firstpixel: usize,
    lastpixel: usize,
    firstline: usize,
    lastline: usize,
    values: Vec<f64>
}


/// Seek the reader to the given "path" 
fn seek_to_list<T, F>(path_list:&[Box<&[u8]>],
                   reader:& mut Reader<&[u8]>,
                   parse_func:F)
                   -> Vec<T>
where F:Fn(&mut Reader<&[u8]>) -> (T){
    
    let mut index:usize = 0;
    let mut result:Vec<T> = Vec::new();
    let mut buf = Vec::new();
    let mut started = false;

    loop {
        match reader.read_event(&mut buf) {
            Ok(Event::Start(ref e)) => {

                // increment the hierarchy
                if index < path_list.len() && *path_list[index] == e.name() {
                    index += 1;
                    started = true;
                }
                
                // XML is malformed
                if index > path_list.len() {
                    panic!("XML file is malformed");
                }
                
                // While we've reached the hierarchy, call the function
                if index == path_list.len() {
                    result.push(parse_func(reader));
                }
            },
            Ok(Event::End(ref e)) => {

                // Decrement the hierarchy
                if started {
                    if *path_list[index-1] == e.name() {
                        index = index - 1 ;
                    }
                    if index == 0 {
                        break;
                    }
                }
            },
            e => {},
        }
    }
    
    return result;
}


impl NoiseField {
    ///Creates a new noise field given the xml file
    ///
    /// Arguments
    ///
    /// * `filebuffer` - The buffer holding the contents of the noise cal xml file
    ///
    /// * `two_eight_flag` - Flag indicates if product is from IPF 2.8 or later
    pub fn new(filebuffer:&str, shape:(usize, usize), two_eight_flag:bool) -> NoiseField {
        let mut reader = Reader::from_str(filebuffer);

        //if two_eight_flag {
        
        let rgst:[Box<&[u8]>;2] = [Box::new(b"noise"),
                                   Box::new(b"noiseRangeVectorList"),
        ];
        
        let azst:[Box<&[u8]>;2] = [Box::new(b"noise"),
                                   Box::new(b"noiseAzimuthVectorList"),
        ];
        
 
        let range_result = seek_to_list(&rgst, &mut reader, NoiseField::range_parse_func);
        reader = Reader::from_str(filebuffer);
        let azimuth_result = seek_to_list(&azst, &mut reader, NoiseField::azimuth_parse_func);
        
        return NoiseField {
            data:NoiseField::compute_field(range_result, azimuth_result, shape)
        };

        
        //     }
       // }
        
      //  else {
       //     let rgst:[Box<&[u8]>;2] = [Box::new(b"noise"),
        //                               Box::new(b"noiseVectorList"),
       //     ];
       //     let range_result = seek_to_list(&rgst, &mut reader, range_parse_func);
       //     return NoiseField {
       //         data:NoiseField::compute_field(range_result, )
       //     };
            
    }



    
    fn old_compute_field(rg_result:Vec<NoiseRangeEntry>) {
        
    }
    
    fn compute_field(rg_result:Vec<NoiseRangeEntry>, az_result:Vec<NoiseAzimuthEntry>, shape:(usize, usize)) -> Array2<f64>{
        let rgarr = NoiseField::interpolate_row(rg_result, shape);
        let azarr = NoiseField::interpolate_col(az_result, shape);
        let mut arr = Array2::zeros(shape);
        Zip::from(&mut arr)
            .and(&rgarr)
            .and(&azarr)
            .apply(|a, &rg, &az| {
            //.par_apply(|a, &rg, &az| {
                *a = rg*az;
            });
        return arr;
    }
    

    fn interpolate_row(rg_result:Vec<NoiseRangeEntry>, shape:(usize, usize)) -> Array2<f64> {
        let mut arr:Array2<f64> = Array2::zeros(shape);

        
        let mut arrlist:Vec<Array1<f64>> = vec![Array1::zeros(shape.1);rg_result.len()]; //this is really bad


        arrlist.par_iter_mut()
            .zip(rg_result.par_iter())
            .for_each(|a|{
                let (val,rg) = a;
                let line = rg.line;
                let mut first_pixel = rg.pixels[0];
                let mut last_pixel:usize;
                let mut first_value = rg.values[0];
                let mut last_value:f64;
                let mut slope:f64;
                for (p,v) in rg.pixels[1..].iter().zip(rg.values[1..].iter()) {
                    last_pixel = *p;
                    last_value = *v;
                    slope = (last_value - first_value) / ((last_pixel - first_pixel) as f64);
                    for i in first_pixel..last_pixel+1 {
                        val[i] = slope*((i-first_pixel) as f64) + first_value;
                    }
                    first_pixel = last_pixel;
                    first_value = last_value;
                        
                }
            });





        //copy the data in the right rows
        for (e,rg) in rg_result.iter().enumerate() {
            arr.index_axis_mut(Axis(0),rg.line).assign(&arrlist[e]);
        }



        //interpolate between the rows
        let mut rg_iter = rg_result.iter();
        let mut start:usize = rg_iter.next().unwrap().line;
        let mut end:usize;
        for rg in rg_iter{
            end = rg.line;
            let endpoints = arr.index_axis(Axis(0),end).into_owned();
            let startpoints = arr.index_axis(Axis(0),start).into_owned();
            Zip::indexed(&mut arr.slice_axis_mut(Axis(0), Slice::from(start..end+1)))
            //.par_apply(|ind, a| {
		.apply(|ind, a| {
                    let (x,y) = ind;
                    let slope = ((&endpoints)[y] - (&startpoints)[y])/((end-start) as f64);
                    *a = slope*(x as f64) + (&startpoints)[y];
                    
                });
            start = end;
        }
        


        
        return arr;
    }

    fn interpolate_col(az_result:Vec<NoiseAzimuthEntry>, shape:(usize, usize)) -> Array2<f64>{
        let mut arr = Array2::zeros(shape);
        let mut arrlist:Vec<Array2<f64>> = az_result
            .iter()
            .map(|a| {
                Array2::zeros((a.lastline-a.firstline+1, a.lastpixel-a.firstpixel+1))
            }).collect();


        arrlist.par_iter_mut()
            .zip(az_result.par_iter())
            .for_each(|a| {
                let (val, az) = a;
                let lowerline = az.firstline;
                let upperline = az.lastline;
                let lowerpixel = az.firstpixel;
                let upperpixel = az.lastpixel;
                
                let mut line_iter = az.lines.iter();
                let mut first_line:usize = *(line_iter.next().unwrap());
                let mut val_iter = az.values.iter();
                let mut first_value:f64 = *(val_iter.next().unwrap());

                
                let mut last_line:usize;
                let mut last_value:f64;
                let mut slope:f64 = 0.0;

                // extrapolate
                if lowerline != first_line {
                    slope = (az.values[1]-first_value)/(az.lines[1] as f64 -first_line as f64);
                    for i in lowerline..first_line+1 {
                        for j in 0..upperpixel-lowerpixel+1 {
                            val[(i-lowerline,j)] = slope * (i as f64 -first_line as f64) + first_value;
                        }
                    }
                }

                // interpolate
                for (l,v) in line_iter.zip(val_iter) {
                    last_line = *l;
                    last_value = *v;

                    slope = (last_value- first_value)/((last_line - first_line) as f64);
                    for i in first_line..last_line+1 {
                        for j in 0..upperpixel-lowerpixel+1 {
                            val[(i-lowerline,j)] = slope*((i-first_line) as f64) + first_value;
                        }
                    }
                    first_line = last_line;
                    first_value = last_value;
                }

                // extrapolate
                if upperline != first_line {
                    //slope remains the same
                    //for i in first_line..upperline+1 {
                    for i in first_line..upperline+1 {
                        for j in 0..upperpixel-lowerpixel+1 {
                            val[(i-lowerline,j)] = slope * (i as f64 -first_line as f64) + first_value;
                        }
                    }
                }
                
            });
        // fill in the gaps (bottleneck)

        for (colarr, az) in arrlist.iter().zip(az_result.iter()){
            let mut slice = arr.slice_mut(s![az.firstline..az.lastline+1, az.firstpixel..az.lastpixel+1]);
            //let bcast = colarr.broadcast((az.lastpixel+1-az.firstpixel,shape.0)).unwrap();
            //slice.assign(&bcast.reversed_axes());
            slice.assign(&colarr);

        }

        
        
        return arr;

    }

    /// Parses a single range noise vector
    fn range_parse_func(reader:&mut Reader<&[u8]>) -> NoiseRangeEntry {
        let mut line = 0;
        let mut buf = Vec::new();
        enum ExtractState {Line, Pixel, Values, None};
        let mut state:ExtractState = ExtractState::None;
        let mut pixels:Vec<usize> = Vec::new();
        let mut values:Vec<f64> = Vec::new();
        loop {
            match reader.read_event(&mut buf) {
                // Determine the current field state
                Ok(Event::Start(ref e)) => {
                    match e.name() {
                        b"line" => {state = ExtractState::Line},
                        b"pixel" => {state = ExtractState::Pixel},
                        b"noiseRangeLut" => {state = ExtractState::Values},
                        _ => {}
                    }
                    
                },
                // Read and parse the text
                Ok(Event::Text(ref e)) =>{
                    match state {
                        ExtractState::Line => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => line = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        },
                        ExtractState::Pixel => {
                            match str::from_utf8(&e.unescaped().unwrap()) {
                                Ok(v) => {
                                    pixels = v.split(|a| a == ' ')
                                        .map(|a| a.parse::<usize>().unwrap())
                                        .collect();
                                }
                                Err(_e) => panic!("Malformed xml file")
                            };

                        },
                        ExtractState::Values => {
                            match str::from_utf8(&e.unescaped().unwrap()) {
                                Ok(v) => {
                                    values = v.split(|a| a == ' ')
                                        .map(|a| a.parse::<f64>().unwrap())
                                        .collect();
                                }
                                Err(_e) => panic!("Malformed xml file")
                            };

                        },
                        ExtractState::None => {
                        }
                    }
                }
                // End the field state
                Ok(Event::End(ref e)) => {
                    match e.name() {

                        b"noiseRangeVector" => {break;},
                        b"line"=> {state = ExtractState::None},
                        b"pixel"=> {state = ExtractState::None},
                        b"noiseRangeLut"=> {state = ExtractState::None},
                        _ => {},
                    }
                }
                _ => ()
            }
        }

        NoiseRangeEntry {
            line: line,
            pixels: pixels,
            values: values
        }
    }

    /// Parses a single range noise vector 
    fn azimuth_parse_func(reader:&mut Reader<&[u8]>) -> NoiseAzimuthEntry {

        let mut buf = Vec::new();
        enum ExtractState {Line, FirstPixel, LastPixel, FirstLine, LastLine, Values, None};
        let mut state:ExtractState = ExtractState::None;
        
        let mut firstpixel:usize = 0;
        let mut lastpixel:usize = 0;
        let mut firstline:usize = 0;
        let mut lastline:usize = 0;
        let mut lines:Vec<usize> = Vec::new();
        let mut values:Vec<f64> = Vec::new();
        loop {
            match reader.read_event(&mut buf) {
                // Determine the current field state
                Ok(Event::Start(ref e)) => {
                    match e.name() {
                        b"line" => {state = ExtractState::Line},
                        b"firstRangeSample" => {state = ExtractState::FirstPixel},
                        b"lastRangeSample" => {state = ExtractState::LastPixel},
                        b"firstAzimuthLine" => {state = ExtractState::FirstLine},
                        b"lastAzimuthLine" => {state = ExtractState::LastLine},
                        b"noiseAzimuthLut" => {state = ExtractState::Values},
                        _ => {}
                    }
                    
                },
                // Read and parse the text
                Ok(Event::Text(ref e)) =>{
                    match state {
                        ExtractState::FirstPixel => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => firstpixel = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        },
                        ExtractState::LastPixel => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => lastpixel = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        },

                        ExtractState::Line => {
                            match str::from_utf8(&e.unescaped().unwrap()) {
                                Ok(v) => {
                                    lines = v.split(|a| a == ' ')
                                        .map(|a| a.parse::<usize>().unwrap())
                                        .collect();
                                }
                                Err(_e) => panic!("Malformed xml file")
                            };

                        },
                        ExtractState::FirstLine => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => firstline = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        },
                        
                        ExtractState::LastLine => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => lastline = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        },

                        
                        ExtractState::Values => {
                            match str::from_utf8(&e.unescaped().unwrap()) {
                                Ok(v) => {
                                    values = v.split(|a| a == ' ')
                                        .map(|a| a.parse::<f64>().unwrap())
                                        .collect();
                                }
                                Err(_e) => panic!("Malformed xml file")
                            };

                        },
                        ExtractState::None => {
                        }
                    }
                }
                // End the field state
                Ok(Event::End(ref e)) => {
                    match e.name() {
                        b"noiseAzimuthVector" => {break;},
                        b"line"=> {state = ExtractState::None},
                        b"firstRangeSample"=> {state = ExtractState::None},
                        b"lastRangeSample"=> {state = ExtractState::None},
                        b"firstAzimuthLine" => {state = ExtractState::None},
                        b"lastAzimuthLine" => {state = ExtractState::None},
                        b"noiseAzimuthLut"=> {state = ExtractState::None},
                        _ => {},
                    }
                }
                _ => ()
            }
        }

        NoiseAzimuthEntry {
            lines: lines,
            firstpixel: firstpixel,
            lastpixel: lastpixel,
            firstline:firstline,
            lastline:lastline,
            values: values
        }
    }


}


impl SwathElem {
    /// Creates an nested vector of swath elems and the period for each subswath
    pub fn new(filebuffer:&str, sentformat:&SentinelFormatId) -> (Vec<Vec<SwathElem>>, Vec<usize>) {
        let mut reader = Reader::from_str(filebuffer);

	let num_subswaths:usize = match sentformat.sentmode.as_str() {
	    "EW" => {5},
	    "IW" => {3},
	    _ => panic!("Incorrect sentmode")
	};

        let swath_keys:[Box<&[u8]>; 3] = [Box::new(b"product"), 
                                          Box::new(b"swathMerging"),
                                          Box::new(b"swathMergeList")];

        let burst_keys:[Box<&[u8]>; 3] = [Box::new(b"product"),
                                             Box::new(b"antennaPattern"),
                                             Box::new(b"antennaPatternList")];
        

        
        // get subswaths
        let swath_bounds = seek_to_list(&swath_keys, &mut reader, SwathElem::parse_swath_elem);
        reader = Reader::from_str(filebuffer);
        let mut w = vec![0_usize;num_subswaths];
        
        // get number of bursts
        let number_count = seek_to_list(&burst_keys, &mut reader, SwathElem::parse_burst);


        let number_burst:Vec<usize> = (1..num_subswaths+1).map(
            |a| number_count.iter().fold(
                0, |acc, x| if *x == a {acc+1} else {acc}) - 1).collect(); //number of antenna elements - 1



        // compute period.

        for a in 0..num_subswaths {
            let min_az = swath_bounds[a].iter().fold(99999, |acc, x| if x.fa <= acc {x.fa} else {acc});// min
            let max_az = swath_bounds[a].iter().fold(0, |acc, x| if x.la >= acc {x.la} else {acc});// max

            w[a] = ((max_az - min_az) / number_burst[a])/2;
        }
        

        return (swath_bounds, w);
    }

    /// At a position, extract a single swathelem
    fn parse_swath_elem(reader:&mut Reader<&[u8]>) -> Vec<SwathElem> {
        let mut swathboundlist:Vec<SwathElem> = Vec::new();
        let mut buf = Vec::new();
        enum SwathState {Swath, FirstAzimuth, FirstRange, LastAzimuth, LastRange, None};
        //let mut current_ss:usize = 1000; // current subswath that is being visited
        let mut state = SwathState::None;

        let mut fa:usize=99999;
        let mut fr:usize=99999;
        let mut la:usize=99999;
        let mut lr:usize=99999;

        loop {
            match reader.read_event(&mut buf) {
                // Set state to the correct mode
                Ok(Event::Start(ref e)) => {
                    match e.name() {
                        b"swath" => {state = SwathState::Swath}, // change swaths
                        b"firstAzimuthLine" => {state = SwathState::FirstAzimuth},
                        b"firstRangeSample" => {state = SwathState::FirstRange},
                        b"lastAzimuthLine" => {state = SwathState::LastAzimuth},
                        b"lastRangeSample" => {state = SwathState::LastRange},
                        _ => {}
                    }
                },

                // Parse numeric values
                Ok(Event::Text(ref e)) => {
                    match state {
                        SwathState::Swath => {
                            /*
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap();
                            if val.len() == 3 {
                                let v = (val.bytes().nth(2).unwrap() - '0' as u8) as usize;
                                current_ss = v;
                            }
                            else{
                                panic!("Malformed xml file");
                            }*/
                        },
                        SwathState::FirstAzimuth => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => fa = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        },
                        SwathState::FirstRange => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => fr = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        },
                        SwathState::LastAzimuth => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => la = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        },
                        SwathState::LastRange => {
                            let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => lr = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        },

                        SwathState::None => {}
                    }
                },

                // Reset the state
                Ok(Event::End(ref e)) => {
                    match e.name() {
                        b"swathMerge" => {break;} // End of data
                        b"firstAzimuthLine" => {state = SwathState::None},
                        b"firstRangeSample" => {state = SwathState::None},
                        b"lastAzimuthLine" => {state = SwathState::None},
                        b"lastRangeSample" => {state = SwathState::None},
                        b"swathBounds" => {
                            // push completed subswath
                            swathboundlist.push(SwathElem{fa:fa, fr:fr, la:la, lr:lr}); 
                            state = SwathState::None;
                        },
                        b"swath" => {state = SwathState::None}, // change swaths
                        _ => {},
                    }
                },

                _ => {}
            }
        }
        return swathboundlist;
    }

    
    fn parse_burst(reader:&mut Reader<&[u8]>) -> usize { //return subswath 
        
        let mut subswath:usize = 0;
        let mut buf = Vec::new();
        enum BurstState {Swath, None};
        let mut state = BurstState::None;
        loop {
            match reader.read_event(&mut buf) {
                Ok(Event::Start(ref e)) => {
                    match e.name() {
                        b"swath" => {state = BurstState::Swath},
                        _ => {},
                    }
                }

                Ok(Event::Text(ref e)) => {
                    match state {
                        BurstState::Swath => {
                            match  str::from_utf8(&e.unescaped().unwrap()) {
                                Ok(val) =>
                                    match val.bytes().nth(2) {
                                        Some(u) => {
                                            let v = (u - ('0' as u8)) as usize;
                                            subswath = v;
                                        }
                                        None => panic!("Malformed xml file")
                                    }
                                Err(_e) => panic!("Malformed xml file")
                            }
                        }
                        BurstState::None => {}
                    }
                }
                    

                Ok(Event::End(ref e)) => {
                    match e.name() {
                        b"swath" => {state = BurstState::None},
                        b"antennaPattern" => {break;},
                        _ => {},
                    }
                },
                
                _ => {}
            }
        }

        subswath
    }
}




#[derive(Debug,Clone)]
 pub struct TimeRow {
     pub row:usize,
     pub col:usize,
     pub elevangle:f64,
     pub aztime:f64
 }



pub struct TimeRowLut {
    pub lut:Vec<Vec<TimeRow>>
}

impl TimeRowLut {
   
    
    /// Time delta of the given azimuth time
  //  fn strtime_to_f64(aztime:&str) -> Result<f64, >{
    fn strtime_to_f64(aztime:&str) -> Result<f64, chrono::format::ParseError>{
	let base =  Utc.ymd(2014, 1, 1).and_hms_milli(0, 0, 0, 0);
	let dt = Utc.datetime_from_str(aztime, "%Y-%m-%dT%H:%M:%S%.f")?;
	let g = dt.signed_duration_since(base);
	match g.to_std() {
	    Ok(s) => Ok(s.as_secs_f64()),
	    Err(e) => panic!("XML file malformed (impossible date time): {}",e)
	}

    }
    ///Return a full lookup table
    pub fn new(filebuffer:&str, swath_bounds:&Vec<Vec<SwathElem>>, sentformat:&SentinelFormatId) -> TimeRowLut{
	let mut reader = Reader::from_str(filebuffer);
	let keys:[Box<&[u8]>;2] = [Box::new(b"geolocationGrid"),
				   Box::new(b"geolocationGridPointList")

	];

	// Parse out TimeRow info from geo-coordinates
	let time_list = seek_to_list(&keys, &mut reader, TimeRowLut::parse_geo);

	let num_subswaths = match sentformat.sentmode.as_str() {
	    "EW" => {5},
	    "IW" => {3},
	    _ => panic!("Incorrect sentmode")
	};

	let mut lut:Vec<Vec<TimeRow>> = vec![Vec::new();num_subswaths];
	let _g = time_list.iter().for_each(|a| (TimeRowLut::find_best_subswath(&a, swath_bounds, &mut lut)));

	TimeRowLut{lut:lut}

    }


    ///Returns the orignal lookuptable from the geocoordinate points.
    fn parse_geo(reader:&mut Reader<&[u8]>) -> TimeRow {
	enum ExtractState { Line, Pixel, ElevationAngle, AzimuthTime, None};
	let mut state = ExtractState::None;

	let mut line:usize = 99999;
	let mut pixel:usize = 99999;
	let mut az_time:f64 = 0.0;
	let mut angle:f64 = 0.0;

	let mut buf = Vec::new();
	
	
	loop {
	    match reader.read_event(&mut buf) {
                Ok(Event::Start(ref e)) => {
                    match e.name() {
                        b"line" => {state = ExtractState::Line},
			b"pixel" => {state = ExtractState::Pixel},
			b"elevationAngle" => {state = ExtractState::ElevationAngle},
			b"azimuthTime" => {state = ExtractState::AzimuthTime},
                        _ => {},
                    }
                }

                Ok(Event::Text(ref e)) => {
                    match state {
                        ExtractState::Line => {
			    let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => line = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        }
			ExtractState::Pixel => {
			    let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<usize>();
                            match val {
                                Ok(v) => pixel = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
                        }
			ExtractState::AzimuthTime => {
			    
			    match TimeRowLut::strtime_to_f64(str::from_utf8(&e.unescaped().unwrap()).unwrap()) {
				Ok(v) => az_time = v,
				Err(_e) => panic!("Malformed xml file (Azimuth time)")
			    }
                        }
			ExtractState::ElevationAngle => {
			    let val = str::from_utf8(&e.unescaped().unwrap()).unwrap().parse::<f64>();
                            match val {
                                Ok(v) => angle = v,
                                Err(_e) => panic!("Malformed xml file")
                            };
			}
                        ExtractState::None => {}
                    }
                }
                    
//		b"geolocationGridPointList" => {break;},
                Ok(Event::End(ref e)) => {
                    match e.name() {
			b"line" => {state = ExtractState::None},
			b"pixel" => {state = ExtractState::None},
			b"azimuthTime" => {state = ExtractState::None},
			b"elevationAngle" => {state = ExtractState::None},
                        b"geolocationGridPoint" => {break},
                        _ => {},
                    }
                },
                
                _ => {}

	    }
	}
	
	TimeRow {
	    row:line,
	    col:pixel,
	    elevangle:angle,
	    aztime:az_time,
	}
    }

    /// Reorders timerow to most appopriate subswath
    fn find_best_subswath(a:&TimeRow, swath_bounds:&Vec<Vec<SwathElem>>, lut:&mut Vec<Vec<TimeRow>>) {
	for (e, swath) in swath_bounds.iter().enumerate() {
	    for s_elem in swath.iter() {
		if a.row >= s_elem.fa && a.row <= s_elem.la && a.col >= s_elem.fr && a.col <= s_elem.lr {
		    lut[e].push(a.clone());
		    return;
		}
	    }
	}

	for (e, swath) in swath_bounds.iter().enumerate() {
	// we did not find the the subswath in bounds. Now search from the most appropriate row.
	    let s_elem = &swath[0];
	    if a.col >= s_elem.fr && a.col <= s_elem.lr {
		lut[e].push(a.clone());
		return;
	    }
	}
	
	panic!("Could not find suitable subswath for geo-coordinates");
    }
    
    fn polyfit_4(x:&[f64], y:&[f64]) -> (f64, f64, f64, f64, f64) {
	let n = x.len(); // number of equations
	let k = 4; // degree of polynomial
	let mut A = vec![0.0;(k+1)*(k+1)];
	let mut b = vec![0.0;(k+1)];

	for i in 0..k+1 {
	    b[i] = x.iter().zip(y.iter()).fold(0.0, |acc, s| acc+(s.0.powi(i as i32))*s.1);
	    for j in 0..k+1 {
		// Note the reverse row/col order to satisfy the fortran routines.
		A[j*(k+1) + i] = x.iter().fold(0.0, |acc, s| acc + s.powi((i+j) as i32));
	    }
	}

	// Solve linear system.
	let mut INFO:i32 = 0;
	let mut IPIV:Vec<i32> = vec![0;k+1];
	// dgesv (square system)
	unsafe {
	    dgesv((k+1) as i32,//num eqs
		  1 as i32,//num eqs
		  &mut A,
		  (k+1) as i32,//leading dim of A
		  &mut IPIV,//pivot matrix
		  &mut b,/////// right hand side
		  (k+1) as i32,//LDB
		  &mut INFO);
	}
	assert!(INFO == 0);
	(b[0], b[1], b[2], b[3], b[4])
    }

    /// Return a closure that is applied to interpolate the values.
    pub fn interp_angle_to_col(&self)  -> Box<Fn(&[f64], usize) -> Vec<f64>>
    {
	// determine unique row entries
	let mut unique_rows = HashSet::new();
	for swath in self.lut.iter() {
	    for ent in swath.iter() {
		unique_rows.insert(ent.row);
	    }
	}
	let mut s_rows:Vec<_> =  unique_rows.into_iter().collect();
	s_rows.sort();
	let mut coef = vec![(0.0,0.0,0.0,0.0,0.0); s_rows.len()];
	    
	for (i,r) in s_rows.iter().enumerate() {
	    let mut ang_arr = Vec::new();
	    let mut col_arr = Vec::new();
	    
	    for swath in self.lut.iter() {
		for ent in swath.iter() {
		    if ent.row == *r {
			ang_arr.push(ent.elevangle);
			col_arr.push(ent.col as f64);
		    }
		}
	    }
	    // fit parameters for a 4th degree polynomial
	   coef[i] = TimeRowLut::polyfit_4(&ang_arr, &col_arr);
	}

	
	fn get_two_best(r_list:&[usize], r:i64) -> (usize, usize) {
	    let mut first_best:i64 = 9999999;
	    let mut first_ind = 0;
	    let mut second_best:i64 = 9999999;
	    let mut second_ind = 0;
	    for (e,a) in r_list.iter().enumerate() {
		let val = ((*a as i64) - r).abs();
		if  val < first_best {
		    second_best = first_best;
		    second_ind = first_ind;
		    first_best = val;
		    first_ind = e;
		}
		else if val < second_best {
		    second_best = val;
		    second_ind = e;
		}
	    }

	    return (first_ind, second_ind);
	    
	    
	}
	let s_rows_ex = s_rows.clone();
	let coef_ex = coef.clone();

	Box::new(move |angle:&[f64], row:usize| {
	     let (a_ind, b_ind) = get_two_best(&s_rows_ex, row as i64);
	     let x_l = angle;
	     let mut cf = coef_ex[a_ind];
	     let v1:Vec<_> = x_l.iter().map(|x| (cf.4)*x.powi(4) + (cf.3)*x.powi(3) + (cf.2)*x.powi(2) + (cf.1)*x.powi(1) + cf.0).collect();
	     cf = coef_ex[b_ind];
	     let v2:Vec<_> = x_l.iter().map(|x| (cf.4)*x.powi(4) + (cf.3)*x.powi(3) + (cf.2)*x.powi(2) + (cf.1)*x.powi(1) + cf.0).collect();

	     let slope:Vec<_> = v2.iter().zip(v1.iter()).map(|z|
						      (z.0 - z.1)/(s_rows_ex[b_ind] as f64 - s_rows_ex[a_ind] as f64)).collect();

	     v1.iter().zip(slope.iter()).map(|z| z.0 + z.1*(row as f64 - s_rows_ex[a_ind] as f64)).collect()

	 })
	
    }
}




#[derive(Clone, Debug)]
pub struct BurstEntry {
    pub fa:usize,
    pub la:usize,
    pub fr:usize,
    pub lr:usize
}

struct RawBurst {
    aztime:f64
}


/// Tracks the givne burst entries.
impl BurstEntry {
    pub fn create_burst_coords (filebuffer:&str, time_row_lut:&TimeRowLut, swath_bounds:&Vec<Vec<SwathElem>>, sentformat:&SentinelFormatId) -> Vec<Vec<BurstEntry>> {
	let mut reader = Reader::from_str(filebuffer);
	let burst_keys:[Box<&[u8]>; 3] = [Box::new(b"product"),
                                          Box::new(b"antennaPattern"),
                                          Box::new(b"antennaPatternList")];

	// Raw burst entries that are still mapped to aztimes
        let raw_burst_entries = seek_to_list(&burst_keys, &mut reader, BurstEntry::parse_burst);

	let num_subswaths:usize = match sentformat.sentmode.as_str() {
	    "EW" => {5},
	    "IW" => {3},
	    _ => panic!("Incorrect sentmode")
	};

	let mut burst_coords:Vec<Vec<BurstEntry>> = vec![Vec::new();num_subswaths];

	for cur_swath in 0..num_subswaths {
	    let mut raw_it = raw_burst_entries.iter().filter(|s| s.0 == cur_swath).peekable();
	    // TODO: ensure that last element is not taken.
	    loop {
		let m = raw_it.next().unwrap();
		let item = &m.1;
		let aztime = item.aztime;
		
		let next_item = raw_it.peek();
		let next_aztime = match next_item {
		    Some(nt) => {nt.1.aztime},
		    None => {break}
		};
		
		let (startrow, endrow) = BurstEntry::lookup_row_from_aztime(aztime, next_aztime, cur_swath, time_row_lut);

		// Identify columns.
		let mut startcol_:Option<usize> = None;
		let mut endcol_:Option<usize> = None;

		for sw_elem in swath_bounds[cur_swath].iter() {
		    if startrow >= sw_elem.fa && startrow <= sw_elem.la {
			startcol_ = Some(sw_elem.fr);
			endcol_ = Some(sw_elem.lr);
		    }
		}

		let minrow = min_row!((swath_bounds[cur_swath]));
		let maxrow = max_row!((swath_bounds[cur_swath]));
		// Couldn't find one in range.
		if startcol_.is_none() {
		    if startrow < minrow {
			startcol_ = Some(swath_bounds[cur_swath][0].fr);
			endcol_ = Some(swath_bounds[cur_swath][0].lr);
		    }
		    else if startrow > maxrow{
			startcol_ = Some(swath_bounds[cur_swath].last().unwrap().fr);
			endcol_ = Some(swath_bounds[cur_swath].last().unwrap().lr);

		    }
		    else {
			panic!("Failure in matching subswaths in burst resolution");
		    }
			
		}

		
		burst_coords[cur_swath].push(
		    BurstEntry {
			fa:startrow,
			la:endrow,
			fr:startcol_.unwrap(),
			lr:endcol_.unwrap()});
	    }
	}
	burst_coords
    }

    /// TODO : lookup the row based on aztime.
    fn lookup_row_from_aztime(aztime:f64, nextaztime:f64, swath:usize, time_row_lut:&TimeRowLut) -> (usize, usize) {

	// Obtain the change in time between each of the lut entries
	let delta_vals:Vec<f64> = time_row_lut.lut[swath].iter().map(|s| aztime - s.aztime ).collect();

	let mut best_ind = 0;
	let mut best_val:f64 = 99999.0;
	let mut found:bool = false;
	for (e,d) in delta_vals.iter().enumerate() {
	    if d.abs() < best_val.abs() {
		best_ind = e;
		best_val = *d;
		found = true;
	    }
	}
	if !found {panic!("Could not resolve lookup table");}
	
	let cur_row = time_row_lut.lut[swath][best_ind].row;

	// find the next closest row with matching column
	let mut next_row = 9999999;
	let mut next_ind = 9999999;
	found = false;
	    
	for (e, xd) in time_row_lut.lut[swath].iter().enumerate() {
            if e == best_ind || xd.row == cur_row {continue}
            if xd.row > cur_row {
		if xd.row < next_row {
                    next_row = xd.row;
		    next_ind = e;
		    found = true;
		}
	    }
	}
	if !found {panic!("Could not resolve lookup table");}

	// find the lowest corresponding time column that matches this row
	// this ensures that our interpolation will be more accurate
	let mut cur_ind = 0;
	let mut cur_col = 10000000;
	let mut next_ind = 0;
	let mut next_col = 10000000;
	found = false;
	let mut found0 = false;
	for (e, sb) in time_row_lut.lut[swath].iter().enumerate() {
            if sb.row == cur_row {
		if sb.col < cur_col{
		    cur_ind = e;
                    cur_col = sb.col;
		    found = true;
		}
	    }
            if sb.row == next_row {
		if sb.col < next_col {
                    next_ind = e;
                    next_col = sb.col;
		    found0 = true;
		}
	    }
	}
	if !found || !found0 {panic!("Could not resolve lookup table");}

	let slope:f64 = (time_row_lut.lut[swath][next_ind].row - time_row_lut.lut[swath][cur_ind].row) as f64/
            (time_row_lut.lut[swath][next_ind].aztime - time_row_lut.lut[swath][cur_ind].aztime);
	
        

	let startrow = (time_row_lut.lut[swath][cur_ind].row as f64 + slope*(aztime - time_row_lut.lut[swath][cur_ind].aztime as f64).max(0.0)) as usize;
	
	let endrow = (time_row_lut.lut[swath][cur_ind].row as f64 + slope*(nextaztime - time_row_lut.lut[swath][cur_ind].aztime) as f64) as usize;
	//Interpolate between rows.

	(startrow, endrow)
    }


    
    /// Parse burst information
    fn parse_burst(reader:&mut Reader<&[u8]>) -> (usize, RawBurst) {
        
        let mut subswath:usize = 0;
	let mut az_time:f64 = 0.0;
        let mut buf = Vec::new();
        enum BurstState {Swath, At, None};
        let mut state = BurstState::None;
        loop {
            match reader.read_event(&mut buf) {
                Ok(Event::Start(ref e)) => {
                    match e.name() {
                        b"swath" => {state = BurstState::Swath},
			b"azimuthTime" => {state = BurstState::At},
                        _ => {},
                    }
                }

                Ok(Event::Text(ref e)) => {
                    match state {
                        BurstState::Swath => {
                            match  str::from_utf8(&e.unescaped().unwrap()) {
                                Ok(val) =>
                                    match val.bytes().nth(2) {
                                        Some(u) => {
                                            let v = (u - ('0' as u8)) as usize;
                                            subswath = v - 1;
                                        }
                                        None => panic!("Malformed xml file")
                                    }
                                Err(_e) => panic!("Malformed xml file")
                            }
                        },
			
			BurstState::At => {
			    match TimeRowLut::strtime_to_f64(str::from_utf8(&e.unescaped().unwrap()).unwrap()) {
				Ok(v) => az_time = v,
				Err(_e) => panic!("Malformed xml file (Azimuth time)")
			    }
			},
			
                        BurstState::None => {}
                    }
                }
                    

                Ok(Event::End(ref e)) => {
                    match e.name() {
                        b"swath" => {state = BurstState::None},
			b"azimuthTime" => {state = BurstState::None},
                        b"antennaPattern" => {break;},
                        _ => {},
                    }
                },
                
                _ => {}
            }
        }
	(subswath, RawBurst{aztime:az_time})
        
    }

}


pub struct RawPattern {
    pub angle:Vec<Vec<Vec<f64>>>,
    pub pattern:Vec<Vec<Vec<f64>>>
}

/// Raw Pattern directly from the xml files.
impl RawPattern {
    pub fn new(filebuffer:&str, sentformat:&SentinelFormatId) -> RawPattern{
	let mut reader = Reader::from_str(filebuffer);
	let burst_keys:[Box<&[u8]>; 3] = [Box::new(b"product"),
                                          Box::new(b"antennaPattern"),
                                          Box::new(b"antennaPatternList")];

	let num_subswaths:usize = match sentformat.sentmode.as_str() {
	    "EW" => {5},
	    "IW" => {3},
	    _ => panic!("Incorrect sentmode")
	};
	
	let raw_pattern_entries = seek_to_list(&burst_keys, &mut reader, RawPattern::parse_burst);
	
	let mut angle:Vec<Vec<Vec<f64>>> = vec![Vec::new(); num_subswaths];
	let mut pattern:Vec<Vec<Vec<f64>>> = vec![Vec::new(); num_subswaths];

	for (ss, ang, patt) in raw_pattern_entries {
	    angle[ss].push(ang.clone());
	    pattern[ss].push(patt.clone());
	}

	RawPattern {
	    angle:angle,
	    pattern:pattern
	}
	
    }

    /// Parse burst information
    fn parse_burst(reader:&mut Reader<&[u8]>) -> (usize, Vec<f64>, Vec<f64>) {
	let ANTNORM:f64 = 43.3_f64.exp();
        let mut subswath:usize = 0;
	let mut angle:Vec<f64> = Vec::new();
	let mut pattern:Vec<f64> = Vec::new();
        let mut buf = Vec::new();
        enum BurstState {Swath, Pattern, Angle, None};
        let mut state = BurstState::None;
        loop {
            match reader.read_event(&mut buf) {
                Ok(Event::Start(ref e)) => {
                    match e.name() {
                        b"swath" => {state = BurstState::Swath},
			b"elevationPattern" => {state = BurstState::Pattern},
			b"elevationAngle" => {state = BurstState::Angle},
                        _ => {},
                    }
                }

                Ok(Event::Text(ref e)) => {
                    match state {
                        BurstState::Swath => {
                            match  str::from_utf8(&e.unescaped().unwrap()) {
                                Ok(val) =>
                                    match val.bytes().nth(2) {
                                        Some(u) => {
                                            let v = (u - ('0' as u8)) as usize;
                                            subswath = v - 1;
                                        }
                                        None => panic!("Malformed xml file")
                                    }
                                Err(_e) => panic!("Malformed xml file")
                            }
                        },
			/* Derive the amplitude of the wave by taking absolute
			of intensity and phase components */
			BurstState::Pattern => {
			    match  str::from_utf8(&e.unescaped().unwrap()) {
                                Ok(val) => {
				    let mut v_vals = val.split_whitespace();
				    loop {
					let real_ = v_vals.next();
					let imag_ = v_vals.next();
					match real_ {
					    Some(real_s) => {
						let imag_s = imag_.unwrap();
						let real = real_s.parse::<f64>().unwrap()/ANTNORM;
						let imag = imag_s.parse::<f64>().unwrap()/ANTNORM;
						
						pattern.push((real*real + imag*imag).sqrt());
					    }
					    None => {
						break;
					    }
					}
					
				    }
				}
                                Err(_e) => panic!("Malformed xml file")
                            }
			}
			
			BurstState::Angle => {
			    match  str::from_utf8(&e.unescaped().unwrap()) {
				Ok(val) => {
				    let mut v_vals = val.split_whitespace();
				    loop {
					let real_ = v_vals.next();
					match real_ {
					    Some(real_s) => {
						let real = real_s.parse::<f64>().unwrap();
						angle.push(real);
					    }
					    None => {
						break;
					    }
					}
					
				    }
				}
                                Err(_e) => panic!("Malformed xml file")
                            }
			}

                        BurstState::None => {}
                    }
                }
                    

                Ok(Event::End(ref e)) => {
                    match e.name() {
                        b"swath" => {state = BurstState::None},
			b"elevationPattern" => {state = BurstState::None},
			b"elevationAngle" => {state = BurstState::None},
                        b"antennaPattern" => {break;},
                        _ => {},
                    }
                },
                
                _ => {}
            }
        }
	(subswath, angle, pattern)
        
    }

}

