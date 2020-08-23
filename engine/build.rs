/*
MIT License

Copyright (c) 2020 Peter Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
extern crate pkg_config;

use std::env;
use std::path::Path;
use std::process::Command;
fn main() {
    match pkg_config::probe_library("scs"){
	Ok(_) => {}
	Err(_e) => {
	    match env::var("SCS_LIB") {
		Ok(s) => {
		    println!("cargo:rustc-link-search=native={}", s);
		},
		Err(_e) => {
		    eprintln!("---------\n---------\n---------\nCouldn't find the default location for scs libraries\nGoing to try the usual spot (i.e. /usr/local/lib).\nIf you have it installed somewhere else, specify this the system environment vars with SCS_INC=path for the directory of the scsdir library\nAlong with SCS_DIR=path for the location of the include files of scs.\n---------\n---------");
		    println!("cargo:rustc-link-search=native=/usr/local/lib/");
			
		}
	    }
	    println!("cargo:rustc-link-lib=static=scsdir");
	}
    }
    
    //let out_dir = env::var("OUT_DIR").unwrap();
    let base = env::var("CARGO_MANIFEST_DIR").unwrap();
    Command::new("make")
	.current_dir(&Path::new(format!("{}/src/c_solve/", base).as_str()))
	.status().unwrap();

    println!("cargo:rustc-link-search=native={}/src/c_solve/", base);
    println!("cargo:rustc-link-lib=static=lp_solve");

    println!("cargo:rerun-if-changed=src/c_solve/scs_solve.c");
    println!("cargo:rerun-if-changed=src/c_solve/scs_solve.h");
    
}
