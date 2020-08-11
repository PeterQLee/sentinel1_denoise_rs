extern crate s1_noisefloor_engine;
use s1_noisefloor_engine::parse::{LinearConfig, HyperParams};
#[macro_use]
extern crate clap;
use clap::{Arg, App};

use s1_noisefloor_engine::interface;

arg_enum! {
    #[derive(Debug)]
    pub enum OpModes {
	LinearEst,
	LinearApply,
	Raw,
	LPEst
    }
}

fn check_path(path:&str) -> bool {
    true
}

fn main() {
    let matches = App::new("Sentinel-1 noisefloor removal")
	.version("1.0")
	.author("Peter Q. Lee <pqjlee@uwaterloo.ca>")
	.about("Sentinel-1 noise floor removal tool. 
Noise floor removal engine for Sentinel-1
Two types of noise removal methods are available.
Note that the engine currently only removes noise floor from cross-pol images.


 1. A linear noise floor removal method that rescales the ESA provided noise floor
    that is provided in each Sentinel-1 product. This is the application of the method in
    P. Q. Lee, L. Xu, D. A. Clausi. 2020. Sentinel-1 additive noise removal from cross-polarization extra-wide TOPSAR with dynamic least-squares. Remote Sensing of Environment. 248. https://doi.org/10.1016/j.rse.2020.111982
    Currently only available in GRD EW mode images.
    OpMode LinearEst for applying the full algorithm.
    OpMode LinearApply can apply the algorithm if you apply the parameters (k)

 2. A non-linear noise floor removal method that computes the noise floor as a power function of
    the antenna pattern. Parameters are estimated with linear programming. This is the application
    of the method in [].
    Can be applied to both EW and IW GRD mode images.
    OpMode LpEst, applies estimation and extraction.

Consists of four different methods (LinearEst, LinearApply, Raw, LPEst).
Allows usage of postprocessing routines for multilooking, converting to amplitude, etc.
Outputs results into an hdf5 file, with groups
x -> processed crosspol
co-> unprocessed copol
k -> estimated values of linear model, 
m -> estimated slope parameters for LP method
b -> estimated intercept parameters for LP method
")
	.arg(Arg::with_name("OpMode")
	     .help("Operational mode.


LinearEst:
  Applies the lstsquares estimation method to retrieve scaling parameters, k,
  rescales the noise floor, y, and subtracts it from the image, x.
  Returns the values back in square intensity units.

LinearApply:
   Applies the linear scaling method using custom user provided scales, k.
   Returns the values back in square intensity units.

Raw:
     Returns the original cross pol, co pol, and noise field from the archive.
     Note that arrays are in linear units.

LpEst:
     Applies the linear programming method to restimate a noise floor based
     on the characteristics of the original image and the 
     Returns the values back in square intensity units.

")
	     .required(true)
	     .possible_values(&OpModes::variants()))
	.arg(Arg::with_name("InputSAR")
	     .help("Input SAR archive or directory")
	     //.value_name("FILE")
	     .required(true))
	.arg(Arg::with_name("OutputHDF")
	     .help("Output hdf5 file")
	     .required(true))
    	.arg(Arg::with_name("config")
	     .short("c")
	     .value_name("Config")
	     .help("Configuration file for default arguments to the algorithm and solvers")
	     .required(false))
	.arg(Arg::with_name("k")
	     .short("k")
	     .help("Comma seperated values for linear scalars to apply")
	     .required(false)).get_matches();
	     

    let config:Option<&str> = matches.value_of("config");
    let opmode:OpModes = value_t!(matches.value_of("OpMode"), OpModes).unwrap_or_else(|e| e.exit());
    let out_path:&str = matches.value_of("OutputHDF").unwrap();
    let sarpath:&str = matches.value_of("InputSAR").unwrap();

    if check_path(out_path){
	match opmode {
	    LinearEst => {linear_get_dualpol_data(sarpath, out_path, config);},
	    LinearApply => {},
	    Raw => {linear_get_raw_data(sarpath, out_path);},
	    LPEst => {}//{lp_get_dualpol_data();}
	}
    }
    else {
	println!("Invalid output path={}. Parent directory does not exist.", out_path);
    }
	     
}


fn linear_get_dualpol_data(archpath:&str, outpath:&str, config:Option<&str>)  {
    let lin_param:LinearConfig = match config {
	Some(s) =>  match LinearConfig::parse_config(s) => {
	    Ok(d) => d,
	    Err(e) => {
		println!("Could not parse config {}",e);
		return exceptions::ValueError.into();
	    }
	},
	None => {LinearConfig::default()}
    };
    
    match interface::linear_get_dualpol_data(archpath, &lin_param) {
	Ok((x, co16, k)) => {
	    //write output directories.
	}
	 Err(e) => {
	     eprintln!("An error occurred. No output {}",e);
	     std::process::exit(1);
	 }
    }
    

}

fn linear_get_customscale_data(zippath:&str,  outpath:&str, py_k:Vec<f64>) {
}

fn linear_get_raw_data(zippath:&str, outpath:&str,) {
    
}

fn lp_get_dualpol_data(zippath:&str,  outpath:&str, lstsq_rescale:bool){
    
}
