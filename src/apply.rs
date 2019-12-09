use crate::parse::SwathElem;
use ndarray::prelude::*;

use ndarray::{ArrayViewMut2, ArrayView1, ArrayView2, Slice};
use ndarray::Zip;
use ndarray_parallel::prelude::*;
use rayon::prelude::*;


const NUM_SUBSWATHS:usize = 5;


/// Applies subswath noise subtraction scaling
///
/// Arguments:
///
/// x: Squared signal array
/// y: noise array
/// k: array of scaling factors
pub fn apply_swath_scale(mut x:ArrayViewMut2<f64>,
                         y:ArrayView2<f64>,
                         k:ArrayView1<f64>,
                         swath_bounds:&[&[SwathElem]]
) {

    //Zip::from(&mut x)
    //    .apply(|x_| {*x_ = (*x_) * (*x_);});

    for a in 0..NUM_SUBSWATHS {
        let ks = k[a];
        for swth in swath_bounds[a].iter() {
            // apply the swath scaling
            Zip::from(&mut x.slice_mut(s![swth.fa..swth.la+1, swth.fr..swth.lr+1]))
                .and(y.slice(s![swth.fa..swth.la+1, swth.fr..swth.lr+1]))
                .apply(|x_, y_| {
                    *x_ = *x_ - ks*y_;
                });
        }
    }
}



pub fn prep_measurement(x:ArrayView2<u16>, mut y:ArrayViewMut2<f64>) -> Array2<f64>{
    let mut result:Array2<f64> = Array2::zeros(x.dim());

    Zip::from(&mut result)
        .and(x)
        .par_apply(|a,b| *a = ((*b as f64) * (*b as f64)));
               
    return result;
}

pub fn restore_scale(mut x:ArrayViewMut2<f64>, mut y:ArrayViewMut2<f64>)  {
    //Zip::from(&mut x)
    //    .apply(|a| *a = (*a)*10000.0);

    //Zip::from(&mut y)
     //   .apply(|a| *a = (*a)*10000.0);
    
}



