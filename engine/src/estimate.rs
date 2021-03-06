//! Estimation procedures specific to linear method.
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
use crate::parse::SwathElem;
use ndarray::prelude::*;
use ndarray::{Array1, Array2, Axis, ArrayViewMut1, ArrayViewMut2, ArrayView1, ArrayView2};
use ndarray::Zip;


// Linear algebra
use lapack::*;
use blas::*;

/// Mean along first axis
fn mean_ax0 (x:ArrayView2<f64>, sa:usize, la:usize, sb:usize, lb:usize) -> Array1<f64>{
    x.slice(s![sa..la, sb..lb])
        .sum_axis(Axis(0)) / 
        ( (la-sa) as f64)
}

///Mean of 1d array
fn mean1(x:ArrayView1<f64>, sa:usize, la:usize) -> f64{
    x.slice(s![sa..la])
        .sum() / 
        ( ( (la-sa) as f64) )
}

/*
///Mean of 2d array
fn mean2(x:ArrayView2<f64>, sa:usize, la:usize, sb:usize, lb:usize) -> f64{
    x.slice(s![sa..la, sb..lb])
        .sum() / 
        ( ( (la-sa) as f64) * ( ( lb-sb) as f64) )
}*/

///argmin of 1d array
fn argmin1(x:ArrayView1<f64>, sa:usize, la:usize) -> usize {
    x.slice(s![sa..la]).indexed_iter()
        .fold((0_usize, 9999999.0), |acc, a| if *a.1 <= acc.1 {(a.0, *a.1)} else {acc} ).0
}

///argmax of 1darray
fn argmax1(x:ArrayView1<f64>, sa:usize, la:usize) -> usize {
    x.slice(s![sa..la]).indexed_iter()
        .fold((0_usize, -9999999.0), |acc, a| if *a.1 >= acc.1 {(a.0, *a.1)} else {acc} ).0
}


const EXTENT:usize = 35;
const NUM_SUBSWATHS:usize = 5;
const NORM:f64 = 10000.0;
    

/// Estimates k values for square of image x and noise field y
/// Assumes that x and y have already been divided by 10000.0 for scaling purposes.
///
/// Arguments
///
/// - `x` square of input hv-SAR image
/// - `y` ESA calibrated noise for hv image
/// - `w` half-period in terms of row spacing
/// - `swath_bounds` boundaries of subswaths in col indices
#[allow(non_snake_case)]
pub fn estimate_k_values(x:ArrayView2<f64>,
                     y:ArrayView2<f64>,
                     w:&[usize],
                     swath_bounds:&[&[SwathElem]],//list
                     mu:f64,
                     gamma:f64,
                     lambda_:&[f64],
                     _lambda2_:f64) -> Array1<f64>{

    // get size of equations
    let n_row_eqs = _num_row_equations(x.dim(), w, swath_bounds);
    let n_inter_eqs = _num_inter_equations(swath_bounds);
    let n_reg = _num_regularization();
    let n_intrasubswath = _num_intra_equations(swath_bounds);

    let n_eqs = n_row_eqs + n_inter_eqs + n_reg + n_intrasubswath;

    // allocate matrices/vectors
    let mut C:Array2<f64> = Array2::zeros((n_eqs, 5));
    let mut m:Array1<f64> = Array1::zeros(n_eqs);

    let n_row_eqs = _num_row_equations(x.dim(), w, swath_bounds);
    // fill matrices and vectors
    _fill_row_equations(m.view_mut(),
                        C.view_mut(),
                        x,
                        y,
                        w,
                        swath_bounds,
                        n_row_eqs);

    _fill_inter_equations(m.view_mut(),
                        C.view_mut(),
                        x,
                        y,
                        swath_bounds,
                        mu,
                        n_row_eqs);

    
    _fill_regularization(m.view_mut(),
                         C.view_mut(),
                         n_row_eqs,
                         n_inter_eqs,
                         lambda_);
    
    _fill_intrasubswath_equations(m.view_mut(),
                                  C.view_mut(),
                                  x,
                                  y,
                                  swath_bounds,
                                  n_row_eqs,
                                  n_inter_eqs,
                                  n_reg,
                                  gamma);

    
    let C_vec:Vec<f64> = C.t().iter().cloned().collect();
    let mut A = vec![0.0;NUM_SUBSWATHS*NUM_SUBSWATHS];

    // Compute C^T@C
    unsafe {
        dgemm(b'T', b'N', NUM_SUBSWATHS as i32, NUM_SUBSWATHS as i32, n_eqs as i32, 1.0, &C_vec,  n_eqs as i32, &C_vec,  n_eqs as i32, 1.0, &mut A,  NUM_SUBSWATHS as i32);
    }

    // Compute C^T@b
    let mut b = vec![0.0;NUM_SUBSWATHS];
    let m_vec = m.to_vec();
    unsafe {
	dgemv(b'T', n_eqs as i32, NUM_SUBSWATHS as i32, 1.0, &C_vec, (n_eqs) as i32, &m_vec, 1_i32, 0.0, &mut b, 1_i32);
    }



    // Solve linear system.
    let mut INFO:i32 = 0;
    let mut IPIV:Vec<i32> = vec![0;NUM_SUBSWATHS];
    unsafe {
        dgesv(NUM_SUBSWATHS as i32, //num eqs
              1 as i32, //num eqs
              &mut A,
              NUM_SUBSWATHS as i32, //leading dim of A
              &mut IPIV, //pivot matrix
              &mut b, /////// right hand side
              NUM_SUBSWATHS as i32, //LDB
              &mut INFO);
    }

    return Array1::from(b);
    
}
                     

fn _num_row_equations(_shape:(usize, usize), w:&[usize], swath_bounds:&[&[SwathElem]]) -> usize{
    let mut n = 0;
    for a in 0..swath_bounds.len() {
        let half_period = w[a];
        let rowlim = swath_bounds[a].iter().fold(0, |ac, y| if y.la+1 > ac { y.la+1 } else {ac}); // find the max row

        for i in 0..swath_bounds[a].len() {
            let swth = &swath_bounds[a][i];
            let mut hf_add = 0;
            if swth.la+1 + half_period >= rowlim {
                hf_add = (swth.la+1 + half_period) - rowlim;
            }
            if swth.fa + hf_add > swth.la+1 {continue;}
            let increment = swth.la+1 - swth.fa - hf_add;
            n += increment;
        }
    }
    
    return n;
}
fn _num_inter_equations(swath_bounds:&[&[SwathElem]]) -> usize{
    let mut tot = 0;

    for a in 0..NUM_SUBSWATHS-1 {
        for swth in swath_bounds[a].iter() {
            if swth.la+1 - swth.fa >= 40 {
                tot+=4;
            }
        }
    }

    return tot;
}

fn _num_regularization() -> usize{
    return NUM_SUBSWATHS;
}
fn _num_intra_equations(swath_bounds:&[&[SwathElem]]) -> usize{
    let mut tot = 0;

    for elem in swath_bounds[0].iter() { //first subswath
        if elem.la+1 - elem.fa >= 40 { tot+=4*4; }
    }

    for a in 1..NUM_SUBSWATHS {
        for elem in swath_bounds[a].iter() {
            if elem.la+1 - elem.fa >= 40 {tot+=2*4}
        }
    }
    
    return tot;
}
#[allow(non_snake_case)]
fn _fill_row_equations(mut m:ArrayViewMut1<f64>,
                       mut C:ArrayViewMut2<f64>,
                       x:ArrayView2<f64>,
                       y:ArrayView2<f64>,
                       w:&[usize],
                       swath_bounds:&[&[SwathElem]],
                       n_row_eqs:usize
) {
    let (x_row, _) = _gather_row(x, w, swath_bounds, n_row_eqs);
    let (y_row, eq_per_swath) = _gather_row(y, w, swath_bounds, n_row_eqs);

    let mut n = 0;



    for a in 0..NUM_SUBSWATHS {
        let i = eq_per_swath[a];

        Zip::from(m.slice_mut(s![n..n+i]))
            .and(C.slice_mut(s![n..n+i, a]))
            .and(x_row.slice(s![n..n+i,0]))
            .and(x_row.slice(s![n..n+i,1]))
            .and(y_row.slice(s![n..n+i,0]))
            .and(y_row.slice(s![n..n+i,1]))
            .apply( |m, c, x0, x1, y0, y1| {
                let m_ = (x0) - (x1);
                let c_ = (y0) - (y1);
                let kratio = c_*m_/(c_*c_);
                if kratio >= 0.0 && kratio <=2.5 {
                    *m = m_;
                    *c = c_;
                }
                    
            });
        

        n+=i;
    }
    


}

#[allow(non_snake_case)]
fn _fill_inter_equations(mut m:ArrayViewMut1<f64>,
                       mut C:ArrayViewMut2<f64>,
                       x:ArrayView2<f64>,
                       y:ArrayView2<f64>,
                       swath_bounds:&[&[SwathElem]],
                       mu:f64,
                       n_row_eqs:usize)
{
    let x_col = _gather_interswathcol(x, swath_bounds);
    let y_col = _gather_interswathcol(y, swath_bounds);
    let _N = _num_inter_equations(swath_bounds);
    let mut n=0;

    for a in 0..NUM_SUBSWATHS-1 {
        let mut n_add = 0;
        for i in 0..swath_bounds[a].len() {
            let swth = &swath_bounds[a][i];
            if swth.la+1-swth.fa >= 40 {
                n_add += 4;
            }
        }


        
        Zip::from(m.slice_mut(s![n_row_eqs + n .. n_row_eqs+ n + n_add]))
            .and(x_col.slice(s![n..n+n_add,0]))
            .and(x_col.slice(s![n..n+n_add,1]))
            .apply(|m_, x0, x1| {
                *m_ = mu*(x0 - x1);
            });
        
        Zip::from(C.slice_mut(s![n_row_eqs + n .. n_row_eqs + n + n_add, a])) 
            .and(y_col.slice(s![n..n+n_add,0]))
            .apply(|ca, y0| {
                *ca = mu*y0;
            });
        
        Zip::from(C.slice_mut(s![n_row_eqs + n .. n_row_eqs + n + n_add, a+1]))
            .and(y_col.slice(s![n..n+n_add,1]))
            .apply(|ca1, y1| {
                *ca1 = -mu*y1;
            });
        
        n += n_add;
    }
}
#[allow(non_snake_case)]
fn _fill_regularization(mut m:ArrayViewMut1<f64>,
                        mut C:ArrayViewMut2<f64>,
                        n_row_eqs:usize,
                        n_col_eqs:usize,
                        lambda_:&[f64]) {
    // Bias all the K's towards 1.
    for i in 0..NUM_SUBSWATHS {
        m[n_row_eqs + n_col_eqs + i] = lambda_[i];
        C[(n_row_eqs + n_col_eqs + i,i)] = lambda_[i];
    }

}
#[allow(non_snake_case)]
fn _fill_intrasubswath_equations(mut m:ArrayViewMut1<f64>,
                                 mut C:ArrayViewMut2<f64>,
                                 x:ArrayView2<f64>,
                                 y:ArrayView2<f64>,
                                 swath_bounds:&[&[SwathElem]],
                                 n_row_eqs:usize,
                                 n_col_eqs:usize,
                                 n_reg:usize,
                                 gamma:f64) {
    let (x_vals, y_vals) = _gather_intrasubswath(x,y, swath_bounds);

    let mut n = 0;
    let _N = _num_intra_equations(swath_bounds);
    //kratio = np.zeros(N)
    for a in 0..NUM_SUBSWATHS {
        
        let mut n_add = 0;
        if a == 0 {
            for swth in swath_bounds[0] {
                    
                if swth.la+1-swth.fa >= 40 {
                    n_add += 4*4;
                }
            }
        }
        else {
            for swth in swath_bounds[a] {
                if swth.la+1-swth.fa >= 40 {
                    n_add += 2*4;
                }
            }
        }
        
        let st = n_row_eqs + n_col_eqs + n_reg + n;
        //m[st: st + n_add] = \
        //   gamma * x_vals[n:n + n_add]

        //C[st: st + n_add,a] = \
        //    gamma * y_vals[n:n + n_add]


        //kratio[n:n+n_add] = C[st:st+n_add,a]*m[st:st+n_add]/np.square(C[st:st+n_add,a])
        Zip::from(m.slice_mut(s![st..st+n_add]))
            .and(C.slice_mut(s![st..st+n_add,a]))
            .and(x_vals.slice(s![n..n+n_add]))
            .and(y_vals.slice(s![n..n+n_add]))
            .apply(|m_, c_, x_, y_| {
                let mval = gamma*x_;
                let cval = gamma*y_;
                let kratio = cval*mval/(cval*cval);
                if kratio >= 0.0 && kratio <=2.5 {
                    *m_ = mval;
                    *c_ = cval;
                }
            });

        

        n+= n_add;
    }


}


#[allow(non_snake_case)]
fn _gather_row(x:ArrayView2<f64>, w:&[usize], swath_bounds:&[&[SwathElem]], n_row_eqs:usize) -> (Array2<f64>, Vec<usize>) {
    //TODO: convert using zip/iterator functions.
    let mut n:usize = 0;
    let mut x_row = Array2::zeros((n_row_eqs, 2));
    let mut eq_per_swath = Vec::with_capacity(NUM_SUBSWATHS);
    for a in 0..NUM_SUBSWATHS {
        let rowlim = swath_bounds[a].iter().fold(0, |ac, y| if y.la+1 > ac { y.la+1 } else {ac}); // find the max row
        let half_period = w[a];
        let o = n;
        for i in 0..swath_bounds[a].len() {
            let swth = &swath_bounds[a][i];

            let mut hf_add:usize = 0;
            if swth.la+1 + half_period >= rowlim {
                hf_add = (swth.la+1+half_period) - rowlim;
            }
            if swth.fa + hf_add > swth.la+1 {continue;}
            let increment = swth.la+1  - swth.fa  - hf_add ;

            Zip::from(&mut x_row.slice_mut(s![n..n+increment, 0]))
                .and(&x.slice(s![swth.fa..swth.la+1-hf_add, swth.fr..swth.lr+1]).sum_axis(Axis(1)))
                //.par_apply(|xr, xm| {
                .apply(|xr, xm| {
                    *xr = xm/((swth.lr+1-swth.fr) as f64)/NORM;
                });

            Zip::from(&mut x_row.slice_mut(s![n..n+increment, 1]))
                .and(&x.slice(s![swth.fa+half_period..swth.la+1+half_period-hf_add,
                                 swth.fr..swth.lr+1]).sum_axis(Axis(1)))
                //.par_apply(|xr, xm|{
		.apply(|xr, xm| {
                    *xr = xm/((swth.lr+1-swth.fr) as f64)/NORM;
                });

                      
            n += swth.la+1 - swth.fa - hf_add;
        }
        eq_per_swath.push(n-o);
        
    }
    return (x_row, eq_per_swath);
}

#[allow(non_snake_case)]
fn _gather_interswathcol(x:ArrayView2<f64>, swath_bounds:&[&[SwathElem]]) -> Array2<f64> {
   
    let N = _num_inter_equations(swath_bounds);
    let mut meanvals = Array2::zeros((N,2));

    let mut sn = 0;
    for a in 0 .. NUM_SUBSWATHS-1 {
        for i in 0..swath_bounds[a].len() { 
            let swth = &swath_bounds[a][i];
            if swth.la+1 - swth.fa < 40 { continue}
            for bnum in 0..4 {
                let step = (swth.la+1-swth.fa)/4;
                let lend:usize;
                if bnum == 3 {
                    lend = swth.la+1;
                }
                else {
                    lend = swth.fa+(bnum+1)*step;
                }

                meanvals[(sn,0)] = x.slice(s![swth.fa + bnum*step.. lend, swth.lr+1-EXTENT..swth.lr+1]).sum()/
                    ( ((lend-(swth.fa+bnum*step)) as f64 )
                        * (swth.lr+1-(swth.lr+1-EXTENT)) as f64)/NORM;
                
                meanvals[(sn,1)] = x.slice(s![swth.fa + bnum*step.. lend, swth.lr+1..swth.lr+1+EXTENT]).sum() /
                    ( ((lend - (swth.fa + bnum*step)) as f64 )
                       * (swth.lr+1+EXTENT - (swth.lr+1)) as f64)/NORM;
            
                sn+=1;
            }
        }
    }
    return meanvals;
}

#[allow(non_snake_case)]
fn _gather_intrasubswath(x:ArrayView2<f64>, y:ArrayView2<f64>, swath_bounds:&[&[SwathElem]]) -> (Array1<f64>, Array1<f64>) {
    
    let pad = 60;
    let bf = 10;
    let Nelems = _num_intra_equations(swath_bounds);
    let mut xvals:Array1<f64> = Array1::zeros(Nelems);
    let mut yvals:Array1<f64> = Array1::zeros(Nelems);
    let mut sn = 0;

    for a in 0..NUM_SUBSWATHS {
        for i in 0..swath_bounds[a].len() {
            let swth = &swath_bounds[a][i];
            // Find the mins and max
            let _n = 0;
            let vy_ = mean_ax0(y.view(),swth.fa,swth.la+1, swth.fr,swth.lr+1)/NORM;
            
            //np.mean(y[fa:la+1, fr:lr+1], axis = 0);
            //let _vx_ = mean_ax0(x.view(),swth.fa,swth.la+1, swth.fr,swth.lr+1)/NORM;
            //np.mean(x[fa:la+1, fr:lr+1], axis = 0);

            if a == 0 {
                // EW has multiple peaks
                
                let Mx0:usize = pad + 20; //local max 1
                let Mx2:usize = swth.lr+1 - swth.fr - pad; //local max 3

                let mn0:usize = argmin1(vy_.view(), 0, (Mx2-Mx0)/2); //local min 1
                
                let mn1:usize = (Mx2-Mx0)/2 + argmin1(vy_.view(), (Mx2-Mx0)/2, vy_.dim() - 2*bf); //local min 2

                let Mx1:usize = mn0 + argmax1(vy_.view(), mn0, mn1); //local max 2

                if swth.la+1 - swth.fa < 40 { continue} //Skip if less than 40 samples
                for bnum in 0..NUM_SUBSWATHS-1{ // #TODO: fill
                    let step = (swth.la+1-swth.fa)/4 ;
                    let lend:usize;
                    if bnum == 3 {
                        lend = swth.la+1;
                    }
                    else {
                        lend = swth.fa+(bnum+1)*step;
                    }
                    
                    let vy = mean_ax0(y, swth.fa+bnum*step, lend, swth.fr, swth.lr+1)/NORM;//np.mean(y[fa+bnum*step:lend, fr:lr+1], axis = 0)
                    let vx = mean_ax0(x, swth.fa+bnum*step, lend, swth.fr, swth.lr+1)/NORM;//np.mean(x[fa+bnum*step:lend, fr:lr+1], axis = 0)

                    xvals[sn] = mean1(vx.view(),Mx0-bf,Mx0+bf) - mean1(vx.view(),mn0-bf,mn0+bf);
                    yvals[sn] = mean1(vy.view(),Mx0-bf,Mx0+bf) - mean1(vy.view(),mn0-bf,mn0+bf);
                    sn+=1;

                    xvals[sn] = mean1(vx.view(),mn0-bf,mn0+bf) - mean1(vx.view(),Mx1-bf,Mx1+bf);
                    yvals[sn] = mean1(vy.view(),mn0-bf,mn0+bf) - mean1(vy.view(),Mx1-bf,Mx1+bf);
                    sn+=1;

                    xvals[sn] = mean1(vx.view(),Mx1-bf,Mx1+bf) - mean1(vx.view(),mn1-bf,mn1+bf);
                    yvals[sn] = mean1(vy.view(),Mx1-bf,Mx1+bf) - mean1(vy.view(),mn1-bf,mn1+bf);
                    sn+=1;

                    xvals[sn] = mean1(vx.view(),mn1-bf,mn1+bf)- mean1(vx.view(),Mx2-bf,Mx2+bf);
                    yvals[sn] = mean1(vy.view(),mn1-bf,mn1+bf) - mean1(vy.view(),Mx2-bf,Mx2+bf);
                    sn+=1;

                }
            }
            else {
                let Mx0 = pad; //#TODO: fill
                let mut Mx1 = swth.lr+1 - swth.fr - pad;
                if a == NUM_SUBSWATHS-1 {
                    Mx1 = swth.lr+1 - swth.fr - 100;
                }
                    
                let mn0 = Mx0 + argmin1(vy_.view(),Mx0,Mx1);
                
                
                if swth.la+1 - swth.fa < 40 {continue;} // #Skip if less than 40 samples

                for bnum in 0..NUM_SUBSWATHS-1 {
                    let step = (swth.la+1-swth.fa)/4; //#TODO: fill
                    let lend:usize;
                    if bnum == 3 {
                        lend = swth.la+1;
                    }
                    else {
                        lend = swth.fa+(bnum+1)*step;
                    }
                    
                    let vy = mean_ax0(y,swth.fa+bnum*step,lend, swth.fr,swth.lr+1)/NORM;
                    let vx = mean_ax0(x,swth.fa+bnum*step,lend, swth.fr,swth.lr+1)/NORM;

                    //assert(not np.all(np.isnan(vy)) and not np.all(np.isnan(vx)))
                    
                    xvals[sn] = mean1(vx.view(),Mx0-bf,Mx0+bf) - mean1(vx.view(),mn0-bf,mn0+bf);
                    yvals[sn] = mean1(vy.view(),Mx0-bf,Mx0+bf) - mean1(vy.view(),mn0-bf,mn0+bf);
                    sn+=1;

                    xvals[sn] = mean1(vx.view(),mn0-bf,mn0+bf) - mean1(vx.view(),Mx1-bf,Mx1+bf);
                    yvals[sn] = mean1(vy.view(),mn0-bf,mn0+bf) - mean1(vy.view(),Mx1-bf,Mx1+bf);
                    sn+=1;
                }
            }
        }
    }
    return (xvals, yvals);
}

