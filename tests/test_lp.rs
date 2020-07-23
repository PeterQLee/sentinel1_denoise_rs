
extern crate denoise_engine;
use denoise_engine::prep_lp::*;
use denoise_engine::est_lp::*;


#[test]
fn test_simple_lp () {
    let lant:&[f64] = &[-9.0,-8.0,-8.0];
    let lreal:&[f64] = &[8.0,9.0,7.0];
    let params = solve_lp(&lreal[..], &lant[..]);
    println!("{} {}", params.m, params.b);
    
}
