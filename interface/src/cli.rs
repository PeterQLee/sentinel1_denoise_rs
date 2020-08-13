extern crate s1_noisefloor_engine;

use s1_noisefloor_engine::parse::{LinearConfig, HyperParams};
#[macro_use]
extern crate clap;
use clap::{Arg, App};
use ndarray;
use std::sync::Arc;
use hdf5;

use s1_noisefloor_engine::interface;

arg_enum! {
    #[derive(Debug)]
    pub enum OpModes {
        LinearEst,
        LinearApply,
        Raw,
        LPEst,
        LPApply
    }
}

fn check_path(path:&str) -> hdf5::Result<hdf5::File> {
    hdf5::File::create(path)
}

fn main() -> hdf5::Result<()> {
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
    OpMode LPEst, applies estimation and extraction.
    OpMode LPApply applies extraction, given parameters

Consists of four different methods (LinearEst, LinearApply, Raw, LPEst, LPApply).
Allows usage of postprocessing routines for multilooking, converting to amplitude, etc.
Outputs results into an hdf5 file, with groups
x -> processed crosspol
co-> unprocessed copol
k -> estimated values of linear model, 
m -> estimated slope parameters for LP method
b -> estimated intercept parameters for LP method
")
	.arg(Arg::with_name("opmode")
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

  LpApply:
     Applies the power function noise floor obtained from linear programming,
     with parameters given by the user.
     Returns the values back in square intensity units.

")
	     .required(true)
	     .possible_values(&OpModes::variants()))
	.arg(Arg::with_name("inputsar")
	     .help("Input SAR archive or directory")
	     .value_name("FILE")
	     .required(true))
	.arg(Arg::with_name("outputhdf")
	     .help("Output hdf5 file")
	     .required(true))
    	.arg(Arg::with_name("config")
	     .short("c")
	     .takes_value(true)
	     .help("Configuration .ini file for default arguments to the algorithm and solvers")
	     .required(false))
	.arg(Arg::with_name("paramhdf")
	     .short("p")
	     .help("HDF5 archive holding the parameters you want to apply for overriding the method\nOnly required for LinearApply")
	     .takes_value(true)
	     .required(false))
	.arg(Arg::with_name("lstsq_rescale")
	     .short("r")
	     .help("Flag for whether to apply least squares rescaling for LP method. Defaults to true")
	     .takes_value(true)
	     .required(false))
	.get_matches();

    

    let config:Option<&str> = matches.value_of("config");
    let paramhdf_s:Option<&str> = matches.value_of("paramhdf");
    let opmode:OpModes = value_t!(matches.value_of("opmode"), OpModes).unwrap_or_else(|e| e.exit());
    let out_path:&str = matches.value_of("outputhdf").unwrap();
    let sarpath:&str = matches.value_of("inputsar").unwrap();


    let lstsq_rescale:bool = value_t!(matches.value_of("lstsq_rescale"), bool).unwrap_or_else(|_e| {true});
    // Debug statement
    if !lstsq_rescale {println!("not estimating least squares");}

    let paramhdf:hdf5::Result<hdf5::File> = match paramhdf_s{
	Some(s) => hdf5::File::open(s),
	None => Err(hdf5::Error::Internal(format!("Could not open parameter file. You must specify an HDF5 files using the -p option.")))
    };

    match check_path(out_path){
	Ok(outf) =>  match opmode {
	    OpModes::LinearEst => {linear_get_dualpol_data(sarpath, outf, config)?;},
	    OpModes::LinearApply => {linear_get_customscale_data(sarpath, outf, paramhdf)?;},
	    OpModes::Raw => {linear_get_raw_data(sarpath, outf)?;},
	    OpModes::LPEst => {lp_get_dualpol_data(sarpath, outf, lstsq_rescale, config)?;},
	    OpModes::LPApply => {lp_get_customscale_data(sarpath, outf, paramhdf, config)?;}
	},
	Err(e) => {
	    eprintln!("Cannot open HDF5 file: {:?}.", e);
	    std::process::exit(1);
	}
    }
    Ok(())
}

	     



fn write_arrayf64<D : ndarray::Dimension>(outf:&hdf5::File, fieldname:&str, arr:ndarray::ArrayView<f64, D>) -> hdf5::Result<()>

{
    let s = arr.shape().to_vec();
    let group = outf.new_dataset::<f64>().create(fieldname, &s)?;
    group.write(arr)?;
    Ok(())
}

fn write_arrayu16<D : ndarray::Dimension>(outf:&hdf5::File, fieldname:&str, arr:ndarray::ArrayView<u16, D>) -> hdf5::Result<()>

{
    let s = arr.shape().to_vec();
    let group = outf.new_dataset::<u16>().create(fieldname, &s)?;
    group.write(arr)?;
    Ok(())
}

macro_rules! check_link {
    ($outf:expr, $g:expr) => {
	if $outf.link_exists($g) {
	    return Err(hdf5::Error::Internal(format!("Cannot write to file. Dataset {} already exists", $g)));
	}
    }	
}

fn linear_get_dualpol_data(archpath:&str, outf:hdf5::File, config:Option<&str>) -> hdf5::Result<()> {
    let lin_param:LinearConfig = match config {
	Some(s) =>  match LinearConfig::parse_config(s) {
	    Ok(d) => d,
	    Err(e) => {
		println!("Could not parse config: {}",e);
		std::process::exit(1);
	    }
	},
	None => {println!("Could not parse config path or was not provided. Using default.");
		 LinearConfig::default()}
    };
    // ensure that none of these groups exist.
    check_link!(outf,"crosspol");
    check_link!(outf,"copol");
    check_link!(outf,"k");
    
    match interface::linear_get_dualpol_data(archpath, &lin_param) {
	Ok((x, co16, k)) => {
	    //write output directories.
	    write_arrayf64(&outf, "crosspol", x.view())?;
	    write_arrayu16(&outf, "copol", co16.view())?;
	    write_arrayf64(&outf, "k", k.view())?;
	},
	 Err(e) => {
	     eprintln!("An error occurred. No output written: {}",e);
	     std::process::exit(1);
	 }
    }
    Ok(())
}

fn read_k(datafile:&hdf5::File) -> hdf5::Result<ndarray::Array1<f64>> {
    let kgroup = datafile.dataset("k")?;
    kgroup.read_1d::<f64>()
}

fn linear_get_customscale_data(archpath:&str, outf:hdf5::File, datafile:hdf5::Result<hdf5::File>) -> hdf5::Result<()> {
    // ensure that none of these groups exist.
    check_link!(outf,"crosspol");
    check_link!(outf,"copol");
    
    // Find the scale data in the datafile
    let df = datafile?;
    let k = read_k(&df)?;
    
    match interface::linear_get_customscale_data(archpath, k.view()) {
	Ok((x, co16)) => {
	    write_arrayf64(&outf, "crosspol", x.view())?;
	    write_arrayu16(&outf, "copol", co16.view())?;
	},
	Err(e) => {
	     eprintln!("An error occurred. No output written: {}",e);
	     std::process::exit(1);
	 }
    }
    Ok(())
}

fn linear_get_raw_data(archpath:&str, outf:hdf5::File) -> hdf5::Result<()> {
    check_link!(outf,"crosspol");
    check_link!(outf,"copol");
    check_link!(outf,"y");
    match interface::linear_get_raw_data(archpath) {
	Ok((x, co16, y)) => {
	    write_arrayu16(&outf, "crosspol", x.view())?;
	    write_arrayu16(&outf, "copol", co16.view())?;
	    write_arrayf64(&outf, "y", y.view())?;
	},
	Err(e) => {
	     eprintln!("An error occurred. No output written: {}",e);
	     std::process::exit(1);
	 }
    }
    Ok(())
}

fn lp_get_dualpol_data(archpath:&str, outf:hdf5::File, lstsq_rescale:bool, config:Option<&str>) -> hdf5::Result<()> {


    let (lin_param, lp_param) = match config{
	Some(s) => {
		//let s = pth.to_string_lossy();
	    (match LinearConfig::parse_config(&s) {
		Ok(d) => d,
		Err(e) => {
		    eprintln!("Error parsing config: {}",e);
		    std::process::exit(1);
		}
	    }, match HyperParams::parse_config(&s) {
		Ok(d) => d,
		Err(e) => {
		    eprintln!("Error parsing config: {}",e);
			std::process::exit(1);
		}
	    })
	}
	None => {
	    println!("Could not parse config path or was not provided. Using default.");
	    (LinearConfig::default(), HyperParams::default())
	}
    };
    check_link!(outf,"crosspol");
    check_link!(outf,"copol");
    check_link!(outf,"m");
    check_link!(outf,"b");
    check_link!(outf,"k");
    check_link!(outf,"subswaths");
    
    match interface::lp_get_dualpol_data(archpath, lstsq_rescale, &lin_param, lp_param) {
	Ok((xv, co16, params, k)) => {
	    let xout = Arc::try_unwrap(xv).expect("Could not unwrap");
	    let xview = xout.to_ndarray();
	    write_arrayf64(&outf, "crosspol", xview)?;
	    write_arrayu16(&outf, "copol", co16.view())?;
	    write_arrayf64(&outf, "k", k.view())?;
	    let m:Vec<f64> = params.iter().map(|i| i.iter()
				      .map(|j| j.m))
		.flatten().collect();
	    let b:Vec<f64> = params.iter().map(|i| i.iter()
				      .map(|j| j.b))
		.flatten().collect();
	    let num_ss:u32 = params.len() as u32;

	    write_arrayf64(&outf, "m", ndarray::ArrayView::from(&m))?;
	    write_arrayf64(&outf, "b", ndarray::ArrayView::from(&b))?;
	    let group = outf.new_dataset::<u32>().create("subswaths", &[1])?;
	    group.write(&[num_ss])?;
	}
	Err(e) => {
	    eprintln!("An error occurred. No output written: {}",e);
	    std::process::exit(1);
	}
    }
    Ok(())
}
fn read_mb(datafile:&hdf5::File) -> hdf5::Result<(Vec<f64>,Vec<f64>,Vec<u32>)> {
    let mgroup = datafile.dataset("m")?;
    let bgroup = datafile.dataset("b")?;
    let sgroup = datafile.dataset("subswaths")?;
    Ok((mgroup.read_1d::<f64>()?.to_vec(),
	bgroup.read_1d::<f64>()?.to_vec(),
	sgroup.read_1d::<u32>()?.to_vec()))
}
fn lp_get_customscale_data(archpath:&str, outf:hdf5::File, datafile:hdf5::Result<hdf5::File>, config:Option<&str>) -> hdf5::Result<()> {

    let lp_param:HyperParams = match config{
	Some(s) => {
	    match HyperParams::parse_config(&s) {
		Ok(d) => d,
		Err(e) => {
		    eprintln!("Error parsing config: {}",e);
		    std::process::exit(1);
		}
	    }
	}
	None => {
	    println!("Could not parse config path or was not provided. Using default.");
	    HyperParams::default()
	}
    };
    check_link!(outf,"crosspol");
    check_link!(outf,"copol");
    
    // Find the scale data in the datafile
    let df = datafile?;
    let (m,b, num_subswaths) = read_mb(&df)?;

    match interface::lp_get_customscale_data(archpath, &m, &b, num_subswaths[0] as usize, lp_param) {
	Ok((xv, co16)) => {
	    let xout = Arc::try_unwrap(xv).expect("Could not unwrap");
	    let xview = xout.to_ndarray();
	    write_arrayf64(&outf, "crosspol", xview)?;
	    write_arrayu16(&outf, "copol", co16.view())?;
	},
	
	Err(e) => {
	    eprintln!("An error occurred. No output written: {}",e);
	    std::process::exit(1);
	}
    }
	
    Ok(())
}
