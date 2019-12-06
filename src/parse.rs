
use quick_xml::Reader;
use quick_xml::events::Event;
use std::str;
use ndarray::{Array, Array1, Array2, ArrayBase, Axis, ArrayViewMut1, ArrayView1, ArrayView2, Slice};
use ndarray::Zip;
use ndarray::prelude::*;
//use ndarray::parallel::prelude::*;
use ndarray_parallel::prelude::*;
use rayon::prelude::*;
//use itertools::Itertools;

use time::PreciseTime; //DEBUG


const NUM_SUBSWATHS:usize = 5;


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
    pub fn new(filebuffer:&str, two_eight_flag:bool) -> NoiseField {
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
                data:NoiseField::compute_field(range_result, azimuth_result, (9992,10400))
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
        //.par_apply(|a, &rg, &az| {
            .par_apply(|a, &rg, &az| {
                *a = rg*az;
            });
        return arr;
    }
    

    fn interpolate_row(rg_result:Vec<NoiseRangeEntry>, shape:(usize, usize)) -> Array2<f64> {
        let mut arr:Array2<f64> = Array2::zeros(shape);

        
        let mut arrlist:Vec<Array1<f64>> = vec![Array1::zeros(shape.1);rg_result.len()]; //this is really bad

        let mut t0 = PreciseTime::now();
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
        let mut t1 = PreciseTime::now();


        println!("First loop {}", t0.to(t1));

        t0 = PreciseTime::now();
        //copy the data in the right rows
        for (e,rg) in rg_result.iter().enumerate() {
            arr.index_axis_mut(Axis(0),rg.line).assign(&arrlist[e]);
        }
        t1 = PreciseTime::now();
        println!("Second loop {}", t0.to(t1));

        t0 = PreciseTime::now();
        //interpolate between the rows
        let mut rg_iter = rg_result.iter();
        let mut start:usize = rg_iter.next().unwrap().line;
        let mut end:usize;
        for rg in rg_iter{
            end = rg.line;
            let endpoints = arr.index_axis(Axis(0),end).into_owned();
            let startpoints = arr.index_axis(Axis(0),start).into_owned();
            Zip::indexed(&mut arr.slice_axis_mut(Axis(0), Slice::from(start..end+1)))
                .par_apply(|ind, a| {
                    let (x,y) = ind;
                    let slope = ((&endpoints)[y] - (&startpoints)[y])/((end-start) as f64);
                    *a = slope*(x as f64) + (&startpoints)[y];
                    
                });
            start = end;
        }
        
        t1 = PreciseTime::now();
        println!("Third loop {}", t0.to(t1));
        println!("{:?}", arr[(0,0)]);
        
        return arr;
    }

    fn interpolate_col(az_result:Vec<NoiseAzimuthEntry>, shape:(usize, usize)) -> Array2<f64>{
        let mut arr = Array2::zeros(shape);
        let mut arrlist:Vec<Array2<f64>> = az_result
            .iter()
            .map(|a| {
                Array2::zeros((a.lastline-a.firstline+1, a.lastpixel-a.firstpixel+1))
            }).collect();
        let mut t0 = PreciseTime::now();

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
                    for i in lowerline..first_line {
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
        let mut t1 = PreciseTime::now();
        println!("Az First loop {}", t0.to(t1));
        // fill in the gaps (bottleneck)
        t0 = PreciseTime::now();
        for (colarr, az) in arrlist.iter().zip(az_result.iter()){
            let mut slice = arr.slice_mut(s![az.firstline..az.lastline+1, az.firstpixel..az.lastpixel+1]);
            //let bcast = colarr.broadcast((az.lastpixel+1-az.firstpixel,shape.0)).unwrap();
            //slice.assign(&bcast.reversed_axes());
            slice.assign(&colarr);

        }
        t1 = PreciseTime::now();
         println!("Az Second loop {}", t0.to(t1));
        
        
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
        //println!("{}\n{:?}\n{:?}\n\n", line, pixels, values);
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
        //println!("{:?}\n{}\n{}\n{}\n{}\n{:?}\n\n", lines, firstpixel, lastpixel, firstline, lastline,  values);
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
    pub fn new(filebuffer:&str) -> (Vec<Vec<SwathElem>>, Vec<usize>) {
        let mut reader = Reader::from_str(filebuffer);

        let swath_keys:[Box<&[u8]>; 3] = [Box::new(b"product"), 
                                          Box::new(b"swathMerging"),
                                          Box::new(b"swathMergeList")];

        let burst_keys:[Box<&[u8]>; 3] = [Box::new(b"product"),
                                             Box::new(b"antennaPattern"),
                                             Box::new(b"antennaPatternList")];
        

        
        // get subswaths
        let swath_bounds = seek_to_list(&swath_keys, &mut reader, SwathElem::parse_swath_elem);
        reader = Reader::from_str(filebuffer);
        let mut w = vec![0_usize;NUM_SUBSWATHS];
        
        // get number of bursts
        let number_count = seek_to_list(&burst_keys, &mut reader, SwathElem::parse_burst);


        let number_burst:Vec<usize> = (1..NUM_SUBSWATHS+1).map(
            |a| number_count.iter().fold(
                0, |acc, x| if *x == a {acc+1} else {acc}) - 1).collect(); //number of antenna elements - 1



        // compute period.

        for a in 0..NUM_SUBSWATHS {
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
