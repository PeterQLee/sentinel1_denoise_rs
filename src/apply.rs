
use crate::parse::{SwathElem,BurstEntry, HyperParams};
use crate::read_from_archive::SentinelFormatId;
use crate::prep_lp::{ArrToArr, TwoDArray, MidPoint, };
use ndarray::prelude::*;

use ndarray::{ArrayViewMut2, ArrayView1, ArrayView2, Slice};
use ndarray::Zip;
//use ndarray_parallel::prelude::*;
use lapack::*;
use blas::*;

use std::sync::Arc;
use std::thread;
const NUM_SUBSWATHS:usize = 5;

macro_rules! get_num_subswath {
    ($sw:expr) => {
	match $sw.sentmode.as_str() {
	    "EW" => {5},
	    "IW" => {3},
	    _ => panic!("Incorrect sentmode")
	}
    }
}

macro_rules! unpack_bound {
    ($sw:expr) => {
	($sw.fa, $sw.la, $sw.fr, $sw.lr)
    }
}

/// Get number of looks for specified image type
/// These values were taken from https://sentinel.esa.int/web/sentinel/user-guides/sentinel-1-sar/resolutions/level-1-ground-range-detected
macro_rules! get_look {
    ($sw:expr) => {
	match $sw.sentmode.as_str() {
	    "EW" => {
		match $sw.quality.as_str() {
		    "H" => {2.7},
		    "M" => {10.7},
		    e => {panic!("Unrecognized image quality {}",e)}
		}
	    }
	    "IW" => {
		match $sw.quality.as_str() {
		    "H" => {4.4},
		    "M" => {81.8},
		    e => {panic!("Unrecognized image quality {}",e)}
		}
	    }
	    e => {panic!("Unrecognized imaging type {}",e)}
	}
    }
}

/// Reduce mean along both axes
fn reduce_mean(x:&TwoDArray,
		fa:usize, la:usize, fr:usize, lr:usize) -> f64
{
    let mut total = 0.0;
    for i in fa..la {
	for j in fr..lr {
	    total += x[(i,j)]/(((la-fa)*(lr-fr)) as f64);
	}
    }
    return total;
}

fn reduce_col_mean(x:&TwoDArray,
		fa:usize, la:usize, fr:usize, lr:usize) -> Vec<f64>  {
    let mut total = vec![0.0;(la-fa)];
    for i in fa..la {
	for j in fr..lr {
	    total[i-fa] += x[(i,j)]/((lr-fr) as f64);
	}
    }
    return total;
}


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
            //.par_apply(|x_, y_| {
		.apply(|x_, y_| {
                    *x_ = (*x_).max(0.0) - ks*y_.max(0.0);
                });
        }
    }
}



pub fn prep_measurement(x:ArrayView2<u16>, mut y:ArrayViewMut2<f64>) -> Array2<f64>{
    let mut result:Array2<f64> = Array2::zeros(x.dim());

    Zip::from(&mut result)
        .and(x)
        //.par_apply(|a,b| *a = ((*b as f64) * (*b as f64)));
	.apply(|a,b| *a = ((*b as f64) * (*b as f64)));
               
    return result;
}

pub fn restore_scale(mut x:ArrayViewMut2<f64>) {//, mut y:ArrayViewMut2<f64>)  {
    Zip::from(&mut x)
        .apply(|a| {
            if (*a).is_sign_negative() {
                *a = (*a).abs().sqrt() * (-1.0);
            }
            else {
                *a = (*a).sqrt();
            }
        });
              

    //Zip::from(&mut y)
     //   .apply(|a| *a = (*a)*10000.0);
    
}



pub fn convert_to_u16_f32(x:ArrayView2<u16>) -> Array2<f32>{
    let mut outarr:Array2<f32> = Array2::zeros(x.dim());

    Zip::from(&mut outarr)
        .and(x)
        .apply(|a,b| *a = *b as f32);
    return outarr;
}

pub fn convert_to_f64_f32(x:ArrayView2<f64>) -> Array2<f32>{
    let mut outarr:Array2<f32> = Array2::zeros(x.dim());

    Zip::from(&mut outarr)
        .and(x)
        .apply(|a,b| *a = *b as f32);
    return outarr;
}


pub struct LpApply{
}


impl LpApply {
    /// Gets gapless antenna pattern.
    fn get_gapless_ant(fr:usize, lr:usize,
		       mp:ArrToArr,
		       midpoint:&Vec<(f64,f64)>,
		       lin_params:&[crate::est_lp::lin_params]) -> Vec<f64> {
	let k:Vec<usize> = (fr..lr).collect();
	let ant = mp(&k);
	let mut p_ant = vec![0.0;lr-fr];
	let num_mp = midpoint.len();

	for (i,m) in midpoint.iter().enumerate() {
	    /* determine start and end */
	    let prev = match i {
		0 => fr,
		_ => fr.max(m.0.round() as usize)
	    };
	    let nex = lr.min(m.1.round() as usize);

	    /* compute p_ant in this undisputed region */
	    p_ant[prev-fr..nex-fr].iter_mut()
		.zip(ant[prev-fr..nex-fr].iter())
		.for_each(|x| *x.0 = lin_params[i].b.exp() * x.1.powf(lin_params[i].m));

	    /* Determine in between region*/
	    if i < num_mp - 1 {
		let s_nex = lr.min(midpoint[i+1].0.round() as usize);
		let slope_prog = (nex..s_nex)
		    .map(|x| {
			let tot = (s_nex - nex) as f64;
			let n = (x - nex) as f64;
			lin_params[i].m*(1.0-n/tot) +
			    lin_params[i+1].m*n/tot});
		let int_prog = (nex..s_nex)
		    .map(|x| {
			let tot = (s_nex - nex) as f64;
			let n = (x - nex) as f64;
			lin_params[i].b*(1.0-n/tot) +
			    lin_params[i+1].b*n/tot});
		p_ant[prev-fr..nex-fr].iter_mut()
		    .zip(ant[prev-fr..nex-fr].iter()
			 .zip(slope_prog.zip(int_prog)))
		    .for_each(|x| *x.0 = ((x.1).1).1.exp() * (x.1).0.powf(((x.1).1).0));
	    }
	}
	/* Interpolate to lr if appropriate*/
	let nex = lr.min(midpoint[num_mp-1].1.round() as usize);
	if nex < lr {
	    p_ant[nex-fr..lr-fr].iter_mut()
		.zip(ant[nex-fr..lr-fr].iter())
		.for_each(|x| *x.0 = lin_params[num_mp-1].b.exp() * x.1.powf(lin_params[num_mp-1].m));
	}
	return p_ant;
    }

    fn subtract_along_burst(x: &mut TwoDArray,
			    p_ant:&Vec<f64>,
			    azimuth_noise:&TwoDArray,
			    fa:usize,
			    la:usize,
			    fr:usize,
			    lr:usize,
			    
    ) {
	for i in fa..la {
	    for j in fr..lr {
		x[(i,j)] = x[(i,j)] - p_ant[j-fr] * azimuth_noise[(i,j)];
	    }
	}
    }


    /// Applies interior noise removal within threads.
    fn apply_interior(x_v:&mut Arc<TwoDArray>,
		      x_m:&mut Arc<TwoDArray>,
		      burst:Arc<BurstEntry>,
		      next_burst:Option<Arc<BurstEntry>>,
		      mp:ArrToArr,
		      sind:Arc<Vec<(f64,f64)>>,
		      lin_par:Vec<crate::est_lp::lin_params>,
		      aznoise:Arc<TwoDArray>,
		      swath_bounds_:Arc<Vec<Vec<SwathElem>>>,
		      cur_burst:usize,
		      num_burst:usize,
		      swath:usize,
    )
		      
    {
	let (b_fa, b_la, _b_fr, _b_lr) = unpack_bound!(burst);
	// TODO: ensure that burst_coords are sorted by fa.
	
	let mut apply_subtract = |fa, la, fr, lr| -> Vec<f64>{
	    let p_ant = LpApply::get_gapless_ant(fr, lr, mp.clone(),
						 &sind,
						 &lin_par);
	    LpApply::subtract_along_burst(unsafe{Arc::get_mut_unchecked(x_v)}, //unsafe is bad, but I guaruntee no race conditions.
					  &p_ant,
					  &aznoise,
					  fa, la, fr, lr);
	    
	    //return pattern used.
	    p_ant
		
	};
	
	// First burst and swath.
	if cur_burst == 0 && swath_bounds_[swath][0].fa < b_fa {
	    let (fa, _la_, fr, lr_) = unpack_bound!(swath_bounds_[swath][0]);
	    let lr = lr_ + 1;
	    apply_subtract(fa, b_fa, fr, lr);
	    // Apply azimuth noise to the p_ant.
	}
	
	// Last subswath.
	if cur_burst == num_burst-1 && swath_bounds_[swath].last().unwrap().fa > b_fa{
	    let (fa, la, fr, lr_) = unpack_bound!(swath_bounds_[swath].last().unwrap());
	    let lr = lr_ + 1;
	    apply_subtract(fa, la, fr, lr);
	}
	
	// Iterate over every swath entry and insert in segments where it is fit.
	for e in 0..swath_bounds_[swath].len() {
	    let (fa, la_, fr, lr_)  = unpack_bound!(swath_bounds_[swath][e]);
	    let la = la_ + 1;
	    let lr = lr_ + 1;
	    
	    if fa > b_la || b_fa > la {continue;}
	    
	    let s_fa = b_fa.max(fa);
	    let s_la = b_la.min(la);
				    
	    let p_ant = apply_subtract(s_fa, s_la, fr, lr);
	    
	    // Check missing row
	    if cur_burst == num_burst-1 && la > b_la { // TODOThis literally doesn't make sense
		/* appply */
		LpApply::subtract_along_burst(unsafe{Arc::get_mut_unchecked(x_m)},
					      &p_ant,
					      &aznoise,
					      b_la, la, fr, lr);
	    }
	    else if cur_burst < num_burst-1 { // Apply between burst coordinates
		let (n_fa, _n_la, _n_fr, _n_lr) = unpack_bound!(next_burst.clone().unwrap());
		if n_fa > b_la && n_fa <= la {
		    /* apply */
		    LpApply::subtract_along_burst(unsafe{Arc::get_mut_unchecked(x_m)},
						  &p_ant,
						  &aznoise,
						  b_la, n_fa, fr, lr);
		    
		}
	    }
	}
    }
    
    /// Applies the power function noise floor
    /// without affine offsets.
    pub fn apply_lp_noisefield(x:Arc<TwoDArray>,
			       azimuth_noise:Arc<TwoDArray>,
			       mp_dict:Vec<Vec<ArrToArr>>,
			       burst_coords:&Vec<Vec<Arc<BurstEntry>>>,
			       split_indices:&MidPoint,
			       swath_bounds:Arc<Vec<Vec<SwathElem>>>,
			       lin_params:&Vec<Vec<crate::est_lp::lin_params>>,
			       hyper:Arc<HyperParams>,
			       id:&SentinelFormatId) -> ()
    {
	let num_subswaths:usize = get_num_subswath!(id);

	
	match split_indices {
	    MidPoint::Est(_) => {
		panic!("Wrong Midpoint. Test segments needed");
	    },
	    MidPoint::Test(splitinds) => {
		// Note that I'm reversing order from the python implementation because it will make
		// it easier to parallelize.
		for swath in 0..num_subswaths {
		    let num_burst = burst_coords[swath].len();

		    // Paralelize this.
		    let mut thread_handles:Vec<_> = (0..num_burst).map(
			|cur_burst| {
			    let mut x_v = x.clone();
			    let mut x_m = x.clone();
			    let burst = burst_coords[swath][cur_burst].clone();
			    
			    let next_burst:Option<Arc<BurstEntry>> =
				if cur_burst == num_burst-1 {None}
			    else {Some(burst_coords[swath][cur_burst+1].clone())};
			    
			    let mp = mp_dict[swath][cur_burst].clone();
			    let sind = splitinds[swath][cur_burst].clone();
			    let lin_par = lin_params[swath].clone();
			    let aznoise = azimuth_noise.clone();
			    let swath_bounds_ = swath_bounds.clone();
			    thread::spawn(move || {
				LpApply::apply_interior(
				    &mut x_v,
				    &mut x_m,
				    burst,
				    next_burst,
				    mp,
				    sind,
				    lin_par,
				    aznoise,
				    swath_bounds_,
				    cur_burst,
				    num_burst,
				    swath);
			    
			    })}).collect();
		    let b = thread_handles.len();
		    for _i in 0..b {thread_handles.pop().unwrap().join().unwrap();}
		}
	    }
	}
    }
    
    /// Premptively compute weights using the original array.
    pub fn compute_weights_for_affine(
	original:&TwoDArray, // probably change to array ref
	burst_coords:&Vec<Vec<Arc<BurstEntry>>>,
	hyper:Arc<HyperParams>,
	id:&SentinelFormatId) -> Vec<Vec<f64>> {
	let mut w_vec = Vec::new();

	let num_subswaths:usize = get_num_subswath!(id);
	let base_l:f64 = get_look!(id); // looks in the sentinel image type

	let lowpad = hyper.affine_lowpad;
	let highpad = hyper.affine_highpad;
	let look:f64 = ((lowpad - highpad) as f64)*base_l;
	
	let var_norm = hyper.affine_var_norm; // Normalization for variance computation

	for swath in 0..num_subswaths-1 {
	    let mut tmp = Vec::new();
	    for st in burst_coords[swath].iter() {
		let (fa, la, _fr, lr) = unpack_bound!(st);
		
		let ref_left = reduce_col_mean(original, fa, la+1, lr-lowpad,lr-highpad);
		let ref_right = reduce_col_mean(original, fa, la+1, lr+1+highpad, lr+1+lowpad);

		let var_left:Vec<_> = ref_left.iter().map(|z| z*z/look).collect();//
		let var_right:Vec<_> = ref_right.iter().map(|z| z*z/look).collect();//

		let combined_var = var_left.iter().zip(var_right.iter()).map(|z| z.0+z.1);
		let valid = var_left.iter().zip(var_right.iter()).map(|z| (*z.0!=0.0) && (*z.1!=0.0));

		let w_:(f64, usize) = combined_var.zip(valid).fold((0.0,0), |acc,z| {
		    if z.1 {return (acc.0+z.0, acc.1+1);}
		    acc});
		let w = var_norm/(w_.0/(w_.1 as f64));//

		tmp.push(w);
	    }
	    w_vec.push(tmp);
	}
	
	w_vec
    }
	
	
	
    /*def square_gapless_segmethod_denoise(x, y, mp_dict, slope_list, intercept_list, burst_coords, swath_bounds,  azimuth_noise, split_indices, sent_mode = 'EW', one_flag = False):*/
    pub fn apply_affine(x:Arc<TwoDArray>,
			w_vec:Vec<Vec<f64>>,
			burst_coords:&Vec<Vec<Arc<BurstEntry>>>,
			swath_bounds:Arc<Vec<Vec<SwathElem>>>,
			hyper:Arc<HyperParams>,
			id:&SentinelFormatId) {
	let num_subswaths:usize = get_num_subswath!(id);
	let o = LpApply::compute_affine(x.clone(),
					w_vec,
					burst_coords,
					hyper.clone(),
					&id);

	assert!(o.len() == num_subswaths);

	let mut handles = Vec::new();
	for swath in 0..num_subswaths {
	    for m in 0..swath_bounds[swath].len() {
		let st = swath_bounds.clone();
		let mut xv = x.clone();
		let oval = o[swath];

		handles.push(thread::spawn(move|| {
		    LpApply::sub_offset(unsafe{Arc::get_mut_unchecked(&mut xv)},
			       st,
			       swath,
			       m,
			       oval);
			       
		}));
	    }
	}
	let b = handles.len();
	for _i in 0..b {handles.pop().unwrap().join().unwrap();}
    }
    /// Adds the computed offset
    fn sub_offset(xv:&mut TwoDArray,
		  swath_bounds:Arc<Vec<Vec<SwathElem>>>,
		  swath:usize,
		  m:usize,
		  oval:f64) {

	let (fa, la_, fr, lr_) = unpack_bound!(swath_bounds[swath][m]);
	let la = la_ + 1;
	let lr = lr_ + 1;
	for i in fa..la {
	    for j in fr..lr {
		xv[(i,j)] = xv[(i,j)] - oval;
	    }
	}
	
    }
    
    /// Least squares compute offset.
    #[allow(non_snake_case)]
    fn compute_affine(x:Arc<TwoDArray>,
		      vec_w:Vec<Vec<f64>>,
		      burst_coords:&Vec<Vec<Arc<BurstEntry>>>,
		      hyper:Arc<HyperParams>,
		      id:&SentinelFormatId) -> Vec<f64>{
			  
	
	let num_subswaths:usize = get_num_subswath!(id);

	let lowpad = hyper.affine_lowpad;
	let highpad = hyper.affine_highpad;

	let LIM = if id.sentmode.as_str() == "EW" {2000000.0}
	else {5000000.0};

	let num_entries:usize = (0..num_subswaths-1).fold(0,|acc1, z1| {
	    acc1 + burst_coords[z1].len()});
	
	let mut m:Vec<f64> = vec![0.0;num_entries+1]; // num swath bounds plus 1 for regularization
	let mut C:Vec<f64> = vec![0.0;(num_entries+1)*num_subswaths]; //(num_entries+1) by num_subswaths matrix

	let mut count:usize = 0;
	for swath in 0..num_subswaths-1 {
	    for (e,st) in burst_coords[swath].iter().enumerate() {
		let (fa, la, _fr, lr) = unpack_bound!(st);
		let left_a = reduce_col_mean(&x, fa, la+1, lr-lowpad, lr-highpad);
		let right_a = reduce_col_mean(&x, fa, la+1, lr+1+highpad, lr+1+lowpad);
		
		let w = vec_w[swath][e];
		
		let num_mas:(f64,usize) = left_a.iter().zip(right_a.iter()).fold((0.0,0),|acc,z| {
		    if (z.0-z.1).abs() < LIM {return (acc.0+z.0-z.1,acc.1+1)}
		    acc});
		
		m[count] = w * num_mas.0/(num_mas.1 as f64);//

		C[swath*(num_entries+1) + count] = w;
		C[(swath+1)*(num_entries+1) + count] = -w;
		count += 1;
	    }
	}


	// Add regularization
	for swath in 0..num_subswaths {C[swath*(num_entries+1) + num_entries] = 1.0;}
	m[num_entries] = 0.0;
	//Solve least squares problem..

	let mut A = vec![0.0;num_subswaths*num_subswaths];
	// Compute C^T@C
	unsafe {
	    dgemm(b'T', b'N', 
		  num_subswaths as i32,
		  num_subswaths as i32,
		  (num_entries+1) as i32,
		  1.0,
		  &C,
		  (num_entries+1) as i32,
		  &C,
		  (num_entries+1) as i32,
		  1.0,
		  &mut A,
		  num_subswaths as i32);
	}
	// Compute C^T@b
	let mut b = vec![0.0;num_subswaths];
	unsafe {
	    dgemv(b'T',
		  (num_entries+1) as i32,
		  num_subswaths as i32,
		  1.0,
		  &C,
		  (num_entries+1) as i32,
		  &m,
		  1_i32,
		  0.0,
		  &mut b,
		  1_i32);
	}


	let mut INFO:i32 = 0;
	let mut IPIV:Vec<i32> = vec![0;num_subswaths];

	unsafe {
            dgesv(num_subswaths as i32, //num eqs
		  1 as i32, //num eqs
		  &mut A,
		  num_subswaths as i32, //leading dim of A
		  &mut IPIV, //pivot matrix
		  &mut b, /////// right hand side
		  num_subswaths as i32, //LDB
		  &mut INFO);
	}
	assert!(INFO==0);
	return b;
    }
    /*
def compute_variance_weighted_subswath_offsets(x, reference, swath_bounds, sent_mode = 'EW'):
*/
}


