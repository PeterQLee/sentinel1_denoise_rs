use quick_xml::Reader;
use quick_xml::events::Event;
use std::str;
use ndarray::{Array, Array1, Array2, ArrayBase, Axis, ArrayViewMut1, ArrayView1, ArrayView2, Slice};
use ndarray::Zip;
use ndarray_parallel::prelude::*;
//use rayon::prelude::*;


struct swath_elem {
    fa:usize,
    la:usize,
    fr:usize,
    lr:usize
}

const EXTENT:usize

/// Estimates k values for image x and noise field y
///
/// Arguments
///
/// - `x` input hv-SAR image
/// - `y` ESA calibrated noise for hv image
/// - `w` half-period in terms of row spacing
/// - `swath_bounds` boundaries of subswaths in col indices
pub fn estimate_k_values(x:Array2<f64>,
                     y:Array2<f64>,
                     w:usize,
                     swath_bounds:&[&[swath_elem]],//list
                     mu:f64,
                     gamma:f64,
                     lambda_:Array1<f64>,
                     lambda2_:f64) {

    // get size of equations
    let n_row_eqs = _num_row_equations(x.shape, w, swath_bounds);
    let n_col_eqs = _num_col_equations(swath_bounds);
    let n_reg = _num_regularization();
    let n_thermal = _num_thermal(swath_bounds);

    let n_eqs = n_row_eqs + n_col_eqs + n_reg + n_thermal;

    // allocate matrices/vectors
    let C = Array2::zeros((n_eqs, 5));
    let m = Array1::zeros((n_eqs));

    // fill matrices and vectors
    
    
    
}
                     

fn _num_row_equations(shape:(usize, usize), w:usize, swath_bounds:&[&[swath_elem]]) -> usize{
    let mut n = 0;
    for a in 0..swath_bounds.len() {
        let half_period = w[a];
        let rowlim = swath_bounds.filter(|x| x > 0);
        
    }
}
fn _num_col_equations(swath_bounds:&[&[swath_elem]]) -> usize{
}

fn _num_regularization() -> usize{
    return 5;
}
fn _num_thermal(swath_bounds:&[&[swath_elem]]) -> usize{
    let mut tot = 0;

    for elem in swath_bounds[0].iter() { //first subswath
        if elem.la - elem.fa >= 40 { tot+=4*4; }
    }

    for a in 1..5{
        for elem in swath_bounds[a].iter() {
            if elem.la - elem.fa >= 40 {tot+=2*4}
        }
    }
    
    return tot;
}

fn _fill_row_equations() {
}
fn _fill_col_equations() {
}
fn _fill_regularizationo() {
}
fn _fill_thermal_equations() {
}

