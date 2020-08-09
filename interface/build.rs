extern crate pkg_config;

use std::env;
use std::path::Path;
use std::process::Command;
fn main() {
    //pkg_config::probe_library("scsdir").unwrap();
    
    let out_dir = env::var("OUT_DIR").unwrap();

    //println!("cargo:rustc-link-lib=scsdir");
    // Command::new("make")
    // 	.current_dir(&Path::new("src/c_solve/"))
    // 	.status().unwrap();

    println!("cargo:rustc-link-search=native=/home/peter/SAR/w_sentinel1_denoise_rs/engine/src/c_solve/");
    println!("cargo:rustc-link-lib=static=lp_solve");
    println!("cargo:rustc-link-search=native=/usr/local/lib/");
    println!("cargo:rustc-link-lib=static=scsdir");
    // println!("cargo:rerun-if-changed=src/c_solve/scs_solve.c");
    // println!("cargo:rerun-if-changed=src/c_solve/scs_solve.h");
    
}
