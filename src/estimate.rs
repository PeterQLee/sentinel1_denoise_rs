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

const EXTENT:usize = 35;
const NUM_SUBSWATH:usize = 5;
    

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
        let rowlim = swath_bounds[a].fold(0, |ac, y| if y > ac { y } else {ac}); // find the max row

        for i in 0..swath_bounds[a].len() {
            let curswth = swath_bounds[a][i];
            let mut hf_add = 0;
            if curswth.la + half_period >= rowlim {
                hf_add = (la + half_period) - rowlim;
            }
            let increment = curswth.la - curswth.fa - hf_add;
            if increment < 0 {continue;}
            n += increment;
        }
    }
    
    return n;
}
fn _num_col_equations(swath_bounds:&[&[swath_elem]]) -> usize{
    let mut tot = 0;

    for a in 0..NUM_SUBSWATH-1 {
        for swth in swath_bounds[a].iter() {
            if swth.la - swth.fa >= 40;
            tot+=4;
        }
    }

    return tot;
}

fn _num_regularization() -> usize{
    return NUM_SUBSWATH;
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

fn _fill_row_equations(m:Array1<f64>,
                       C:Array2<f64>,
                       x:Array2<f64>,
                       y:Array2<f64>,
                       w:&[usize],
                       swath_bounds:&[&[swath_elem]],
                       n_row_eqs:usize
) {
    let x_row, _ = _gather_row(x, w, swath_bounds, n_row_eqs);
    let y_row, eq_per_swath = _gather_row(y, w, swath_bounds, n_row_eqs);

    let mut n = 0;
    m[:n_row_eqs] = x_row[:,0] - x_row[:,1];


    for a in 0..NUM_SUBSWATH {
        let i = eq_per_swath[a];

        Zip::from(m.mut_slice(s![n..n+i]))
            .and(C.mut_slice(s![n..n+i, a]))
            .and(x_row.slice(s![n..n+1,0]))
            .and(x_row.slice(s![n..n+1,1]))
            .and(y_row.slice(s![n..n+1,0]))
            .and(y_row.slice(s![n..n+1,1]))
            .apply( |m, c, x0, x1, y0, y1| {
                let m_ = x0 - x1;
                let c_ = y0 - y1;
                let kratio = c_*m_/(c_*c_);
                if kratio >=0 && kratio <=2.5 {
                    *m = m_;
                    *c = c_;
                }
                    
            });
        
        //C[n:n+i, a] = y_row[n:n+i,0] - y_row[n:n+i,1];
        //kratio[n:n+i] = C[n:n+i,a]*m[n:n+i]/np.square(C[n:n+i,a]);
        n+=i;
    }
    
    // Mask according to ratio
    //m[:n_row_eqs][(kratio < 0) | (kratio > 2.5)] = 0;
    //C[:n_row_eqs][(kratio < 0) | (kratio > 2.5)] = 0;

}

fn _fill_col_equations() {
    let x_col = _gather_interswathcol(x, swath_bounds, EXTANT);
    let y_col = _gather_interswathcol(y, swath_bounds, EXTANT);
    let N = _num_col_equations(EXTANT, swath_bounds);
    let mut n=0;

    let kratio = Array1::zeros(N);

    for a in 0..len(swath_bounds)-1 {
        let mut n_add = 0;
        for i in 0..swath_bounds[a].len() {
            let swth = swath_bounds[a][i];
            if swth.la-swth.fa >= 40 {
                n_add += 4;
            }
        }
        //m[n_row_eqs + n:n_row_eqs + n + n_add] = mu*(x_col[n:n+n_add,0] - x_col[n:n+n_add,1]);
        //C[n_row_eqs + n:n_row_eqs + n + n_add, a] = mu*y_col[n:n+n_add,0];
        //C[n_row_eqs + n:n_row_eqs + n + n_add, a+1] = -mu*y_col[n:n+n_add,1];

        Zip::from(m.mut_slice(s![n_row_eqs + n .. n_row_eqs+ n + n_add]))
            .and(C.mut_slice(s![n_row_eqs + n .. n_row_eqs + n + n_add, a])) //TODO: PROBLEM multiple mut borrow PROBLEM
            .and(C.mut_slice(s![n_row_eqs + n .. n_row_eqs + n + n_add, a+1]))
            .and(x_col.slice(s![n..n+n_add,0]))
            .and(x_col.slice(s![n..n+n_add,1]))
            .and(y_col.slice(s![n..n+n_add,0]))
            .and(y_col.slice(s![n..n+n_add,1]))
            .apply(|m_, ca, ca1, x0, x1, y0, y1| {
                *m_ = mu*(x0 - x1);
                *ca = mu*y0;
                *ca1 = -mu*y1
            });
        n += n_add;
    }
}
fn _fill_regularizationo() {
}
fn _fill_thermal_equations() {
}



fn _gather_row(x:Array2<f64>, w:&[usize], swath_bounds:&[&[swath_elem]], n_row_eqs:usize) -> (Array2<f64>, Vec<usize>) {
    //TODO: convert using zip/iterator functions.
    let mut n = 0;
    let mut x_row_eqs = Array2::zeros((n_row_eqs, 2));
    let eq_per_swath = Vec::with_capacity(NUM_SUBSWATHS)
    for a in 0..NUM_SUBSWATHS {
        let rowlim = swath_bounds[a].fold(0, |ac, y| if y > ac { y } else {ac}); // find the max row
        let half_period = w[a];
        let o = n;
        for i in 0..swath_bounds[a].len() {
            let swth = swath_bounds[a][i];

            let mut hf_add = 0;
            if swth.la + half_period >= rowlim {
                hf_add = (la+half_period) - rowlim;
            }
            increment = swth.la - swth.fa - hf_add;
            if increment < 0 {continue;}

            //x_row[n:n+increment, 0 ] =  np.mean(x[fa:la-hf_add, fr:lr], axis = 1)
            //x_row[n:n+increment, 1 ] =  np.mean(x[fa + half_period:la + half_period - hf_add, fr:lr], axis=1)
            Zip::from(&mut x_row.mut_slice(s![n..n+increment, 0]))
                .and(&x.slice(s![swth.fa..swth.la-hf_add, fr..lr]).sum_axis(Axis(1)))
                .apply(|xr, xm| {
                    *xr = xm/(fr-lr as f64);
                });

            Zip::from(&mut x_row.mut_slice(s![n..n+increment, 1]))
                .and(&x.slice(s![swth.fa+half_period..swth.la+half_period-hf_add,
                                 fr..lr]).sum_axis(Axis(1)))
                .apply(|xr, xm| {
                    *xr = xm/(fr-lr as f64);
                });

                      
            n += la - fa -hf_add;
        }
        
    }
    return (x_row_eqs, eq_per_swath);
}


fn gather_interswathcol(x:Array2<f64>, swath_bounds:&[&[swath_elem]]) {
   
    let N = _num_col_equations(swath_bounds);
    let meanvals = Array2::zeros((N,2));

    let mut sn = 0;
    for a in 0 .. NUM_SUBSWATHS-1 {
        for i in 0..swath_bounds[a].len() { 
            let swth = swath_bounds[a];
            if swth.la - swth.fa < 40 { continue}
            for bnum in 0..4 {
                let step = (swth.la+1-swth.fa)/4;
                let lend:usize;
                if bnum == 3:
                    lend = swth.la+1;
                else:
                    lend = swth.fa+(bnum+1)*step;

                meanvals[(sn,0)] = x.slice(s![swth.fa + bnum*step.. lend, swth.lr-EXTANT..swth.lr]).sum()/ ( (lend-(swth.fa+bnum*step)) * (swth.lr-(swth.lr-EXTANT)) as f64);
                meanvals[(sn,1)] = x.slice(s![swth.fa + bnum*step.. lend, swth.lr+1..swth.lr+1+EXTANT]).sum() / ( (lend - (swth.fa + bnum*step)) * (swth.lr+1+EXTANT - (swth.lr+1)) as f64);
                
                sn+=1;
            }
    return meanvals;
}
