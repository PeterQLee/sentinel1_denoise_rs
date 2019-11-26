use quick_xml::Reader;
use quick_xml::events::Event;
use std::str;
use ndarray::{Array, Array1, Array2, ArrayBase, Axis, ArrayViewMut1, ArrayViewMut2, ArrayView1, ArrayView2, Slice};
use ndarray::Zip;
use ndarray_parallel::prelude::*;
//use rayon::prelude::*;


pub struct SwathElem {
    fa:usize,
    la:usize,
    fr:usize,
    lr:usize
}

const EXTENT:usize = 35;
const NUM_SUBSWATHS:usize = 5;
    

/// Estimates k values for image x and noise field y
///
/// Arguments
///
/// - `x` input hv-SAR image
/// - `y` ESA calibrated noise for hv image
/// - `w` half-period in terms of row spacing
/// - `swath_bounds` boundaries of subswaths in col indices
pub fn estimate_k_values(x:ArrayView2<f64>,
                     y:ArrayView2<f64>,
                     w:&[usize],
                     swath_bounds:&[&[SwathElem]],//list
                     mu:f64,
                     gamma:f64,
                     lambda_:Array1<f64>,
                     lambda2_:f64) {

    // get size of equations
    let n_row_eqs = _num_row_equations(x.dim(), w, swath_bounds);
    let n_col_eqs = _num_col_equations(swath_bounds);
    let n_reg = _num_regularization();
    let n_thermal = _num_thermal(swath_bounds);

    let n_eqs = n_row_eqs + n_col_eqs + n_reg + n_thermal;

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

    _fill_col_equations(m.view_mut(),
                        C.view_mut(),
                        x,
                        y,
                        swath_bounds,
                        mu,
                        n_row_eqs);

    
    
}
                     

fn _num_row_equations(shape:(usize, usize), w:&[usize], swath_bounds:&[&[SwathElem]]) -> usize{
    let mut n = 0;
    for a in 0..swath_bounds.len() {
        let half_period = w[a];
        let rowlim = swath_bounds[a].iter().fold(0, |ac, y| if y.la > ac { y.la } else {ac}); // find the max row

        for i in 0..swath_bounds[a].len() {
            let swth = &swath_bounds[a][i];
            let mut hf_add = 0;
            if swth.la + half_period >= rowlim {
                hf_add = (swth.la + half_period) - rowlim;
            }
            if swth.fa + hf_add > swth.la {continue;}
            let increment = swth.la - swth.fa - hf_add;
            n += increment;
        }
    }
    
    return n;
}
fn _num_col_equations(swath_bounds:&[&[SwathElem]]) -> usize{
    let mut tot = 0;

    for a in 0..NUM_SUBSWATHS-1 {
        for swth in swath_bounds[a].iter() {
            if swth.la - swth.fa >= 40 {
                tot+=4;
            }
        }
    }

    return tot;
}

fn _num_regularization() -> usize{
    return NUM_SUBSWATHS;
}
fn _num_thermal(swath_bounds:&[&[SwathElem]]) -> usize{
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

    //m[:n_row_eqs] = x_row[:,0] - x_row[:,1];


    for a in 0..NUM_SUBSWATHS {
        let i = eq_per_swath[a];

        Zip::from(m.slice_mut(s![n..n+i]))
            .and(C.slice_mut(s![n..n+i, a]))
            .and(x_row.slice(s![n..n+1,0]))
            .and(x_row.slice(s![n..n+1,1]))
            .and(y_row.slice(s![n..n+1,0]))
            .and(y_row.slice(s![n..n+1,1]))
            .apply( |m, c, x0, x1, y0, y1| {
                let m_ = x0 - x1;
                let c_ = y0 - y1;
                let kratio = c_*m_/(c_*c_);
                if kratio >= 0.0 && kratio <=2.5 {
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

fn _fill_col_equations(mut m:ArrayViewMut1<f64>,
                       mut C:ArrayViewMut2<f64>,
                       x:ArrayView2<f64>,
                       y:ArrayView2<f64>,
                       swath_bounds:&[&[SwathElem]],
                       mu:f64,
                       n_row_eqs:usize)
{
    let x_col = _gather_interswathcol(x, swath_bounds);
    let y_col = _gather_interswathcol(y, swath_bounds);
    let N = _num_col_equations(swath_bounds);
    let mut n=0;

    for a in 0..NUM_SUBSWATHS-1 {
        let mut n_add = 0;
        for i in 0..swath_bounds[a].len() {
            let swth = &swath_bounds[a][i];
            if swth.la-swth.fa >= 40 {
                n_add += 4;
            }
        }
        //m[n_row_eqs + n:n_row_eqs + n + n_add] = mu*(x_col[n:n+n_add,0] - x_col[n:n+n_add,1]);
        //C[n_row_eqs + n:n_row_eqs + n + n_add, a] = mu*y_col[n:n+n_add,0];
        //C[n_row_eqs + n:n_row_eqs + n + n_add, a+1] = -mu*y_col[n:n+n_add,1];

        
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
fn _fill_regularization(mut m:ArrayViewMut1<f64>,
                        mut C:ArrayViewMut2<f64>,
                        n_row_eqs:usize,
                        n_col_eqs:usize,
                        lambda_:ArrayView1<f64>) {
    // Bias all the K's towards 1.
    for i in 0..NUM_SUBSWATHS {
        m[n_row_eqs + n_col_eqs + i] = lambda_[i];
        C[(n_row_eqs + n_col_eqs + i,i)] = lambda_[i];
    }

}
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

    let n = 0;
    let N = _num_thermal(swath_bounds);
    //kratio = np.zeros(N)
    for a in 0..NUM_SUBSWATHS {
        
        let mut n_add = 0;
        if a == 0 {
            for &swth in swath_bounds[0] {
                    
                if swth.la-swth.fa >= 40 {
                    n_add += 4*4;
                }
            }
        }
        else {
            for &swth in swath_bounds[a] {
                if swth.la-swth.fa >= 40 {
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
                if kratio >=0 && kratio <=2.5 {
                    *m_ = mval;
                    c_ = cval;
                }
            });

        

        n+= n_add;
    }


}



fn _gather_row(x:ArrayView2<f64>, w:&[usize], swath_bounds:&[&[SwathElem]], n_row_eqs:usize) -> (Array2<f64>, Vec<usize>) {
    //TODO: convert using zip/iterator functions.
    let mut n:usize = 0;
    let mut x_row = Array2::zeros((n_row_eqs, 2));
    let eq_per_swath = Vec::with_capacity(NUM_SUBSWATHS);
    for a in 0..NUM_SUBSWATHS {
        let rowlim = swath_bounds[a].iter().fold(0, |ac, y| if y.la > ac { y.la } else {ac}); // find the max row
        let half_period = w[a];

        for i in 0..swath_bounds[a].len() {
            let swth = &swath_bounds[a][i];

            let mut hf_add:usize = 0;
            if swth.la + half_period >= rowlim {
                hf_add = (swth.la+half_period) - rowlim;
            }
            if swth.fa + hf_add > swth.la {continue;}
            let increment = swth.la  - swth.fa  - hf_add ;

            //x_row[n:n+increment, 0 ] =  np.mean(x[fa:la-hf_add, fr:lr], axis = 1)
            //x_row[n:n+increment, 1 ] =  np.mean(x[fa + half_period:la + half_period - hf_add, fr:lr], axis=1)
            Zip::from(&mut x_row.slice_mut(s![n..n+increment, 0]))
                .and(&x.slice(s![swth.fa..swth.la-hf_add, swth.fr..swth.lr]).sum_axis(Axis(1)))
                .apply(|xr, xm| {
                    *xr = xm/((swth.fr-swth.lr) as f64);
                });

            Zip::from(&mut x_row.slice_mut(s![n..n+increment, 1]))
                .and(&x.slice(s![swth.fa+half_period..swth.la+half_period-hf_add,
                                 swth.fr..swth.lr]).sum_axis(Axis(1)))
                .apply(|xr, xm| {
                    *xr = xm/((swth.fr-swth.lr) as f64);
                });

                      
            n += swth.la - swth.fa - hf_add;
        }
        
    }
    return (x_row, eq_per_swath);
}


fn _gather_interswathcol(x:ArrayView2<f64>, swath_bounds:&[&[SwathElem]]) -> Array2<f64> {
   
    let N = _num_col_equations(swath_bounds);
    let mut meanvals = Array2::zeros((N,2));

    let mut sn = 0;
    for a in 0 .. NUM_SUBSWATHS-1 {
        for i in 0..swath_bounds[a].len() { 
            let swth = &swath_bounds[a][i];
            if swth.la - swth.fa < 40 { continue}
            for bnum in 0..4 {
                let step = (swth.la+1-swth.fa)/4;
                let lend:usize;
                if bnum == 3 {
                    lend = swth.la+1;
                }
                else {
                    lend = swth.fa+(bnum+1)*step;
                }

                meanvals[(sn,0)] = x.slice(s![swth.fa + bnum*step.. lend, swth.lr-EXTENT..swth.lr]).sum()/
                    ( ((lend-(swth.fa+bnum*step)) as f64 )
                        * (swth.lr-(swth.lr-EXTENT)) as f64);
                
                meanvals[(sn,1)] = x.slice(s![swth.fa + bnum*step.. lend, swth.lr+1..swth.lr+1+EXTENT]).sum() /
                    ( ((lend - (swth.fa + bnum*step)) as f64 )
                       * (swth.lr+1+EXTENT - (swth.lr+1)) as f64);
            
                sn+=1;
            }
        }
    }
    return meanvals;
}


fn _gather_intrasubswath(x:ArrayView2<f64>, y:ArrayView2<f64>, swath_bounds:&[&[SwathElem]]) -> (Array1<f64>, Array1<f64>) {
    
    let pad = 60;
    let bf = 10;
    let Nelems = _num_intrasubswath(swath_bounds);
    let xvals:Array1<f64> = Array1::zeros(Nelems);
    let yvals:Array1<f64> = Array1::zeros(Nelems);
    let mut sn = 0;

    for a in 0..NUM_SUBSWATHS {
        for i in 0..swath_bounds[a].len() {
            let swth = &swath_bounds[a][i];
            // Find the mins and max
            let mut n = 0;
            let vy_ = np.mean(y[fa:la+1, fr:lr+1], axis = 0);
            let vx_ = np.mean(x[fa:la+1, fr:lr+1], axis = 0);

            if a == 0 {
                // EW has multiple peaks
                
                let Mx0 = pad + 20;
                let Mx2 = lr - fr - pad;

                let mn0 = np.argmin(vy_[0:(Mx2-Mx0)/2]);
                let mn1 = (Mx2-Mx0)/2 + np.argmin(vy_[(Mx2-Mx0)/2:])

                Mx1 = mn0 + np.argmax(vy_[mn0:mn1])

                if la - fa < 40: continue #Skip if less than 40 samples
                for bnum in range(4): #TODO: fill
                    step = (la+1-fa)//4
                    if bnum == 3:
                        lend = la+1
                    else:
                        lend = fa+(bnum+1)*step
                    
                    vy = np.mean(y[fa+bnum*step:lend, fr:lr+1], axis = 0)
                    vx = np.mean(x[fa+bnum*step:lend, fr:lr+1], axis = 0)

                    xvals[sn] = np.mean(vx[Mx0-bf:Mx0+bf] - vx[mn0-bf:mn0+bf])
                    yvals[sn] = np.mean(vy[Mx0-bf:Mx0+bf] - vy[mn0-bf:mn0+bf])
                    sn+=1

                    xvals[sn] = np.mean(vx[mn0-bf:mn0+bf] - vx[Mx1-bf:Mx1+bf])
                    yvals[sn] = np.mean(vy[mn0-bf:mn0+bf] - vy[Mx1-bf:Mx1+bf])
                    sn+=1

                    xvals[sn] = np.mean(vx[Mx1-bf:Mx1+bf] - vx[mn1-bf:mn1+bf])
                    yvals[sn] = np.mean(vy[Mx1-bf:Mx1+bf] - vy[mn1-bf:mn1+bf])
                    sn+=1

                    xvals[sn] = np.mean(vx[mn1-bf:mn1+bf] - vx[Mx2-bf:Mx2+bf])
                    yvals[sn] = np.mean(vy[mn1-bf:mn1+bf] - vy[Mx2-bf:Mx2+bf])
                    sn+=1



            else:
                Mx0 = pad#TODO: fill
                Mx1 = lr - fr - pad
                if a == 4:
                    Mx1 = lr - fr - 100
                    
                mn0 = Mx0 + np.argmin(vy_[Mx0:Mx1])
                
                
                if la - fa < 40: continue #Skip if less than 40 samples

                for bnum in range(4):
                    step = (la+1-fa)//4#TODO: fill
                    if bnum == 3:
                        lend = la+1
                    else:
                        lend = fa+(bnum+1)*step
                    vy = np.mean(y[fa+bnum*step:lend, fr:lr+1], axis = 0)
                    vx = np.mean(x[fa+bnum*step:lend, fr:lr+1], axis = 0)

                    assert(not np.all(np.isnan(vy)) and not np.all(np.isnan(vx)))
                    xvals[sn] = np.mean((vx[Mx0-bf:Mx0+bf] - vx[mn0-bf:mn0+bf]))
                    yvals[sn] = np.mean((vy[Mx0-bf:Mx0+bf] - vy[mn0-bf:mn0+bf]))
                    sn+=1

                    xvals[sn] = np.mean((vx[mn0-bf:mn0+bf] - vx[Mx1-bf:Mx1+bf]))
                    yvals[sn] = np.mean((vy[mn0-bf:mn0+bf] - vy[Mx1-bf:Mx1+bf]))
                    sn+=1

    return xvals, yvals
                }
