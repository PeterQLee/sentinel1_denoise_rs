use crate::parse::{SwathElem,BurstEntry, HyperParams};
use crate::read_from_archive::SentinelFormatId;
use crate::prep_lp::{ArrToArr, TwoDArray, MidPoint, };
use ndarray::prelude::*;

use ndarray::{ArrayViewMut2, ArrayView1, ArrayView2, Slice};
use ndarray::Zip;
//use ndarray_parallel::prelude::*;
use rayon::prelude::*;

use std::sync::Arc;
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
                    *x_ = *x_ - ks*y_;
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

macro_rules! unpack_bound {
    ($sw:expr) => {
	($sw.fa, $sw.la, $sw.fr, $sw.lr)
    }
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
	    println!("pantsize = {}, max={}",p_ant.len(), nex-fr);
	    println!("antsize = {}, max={}",ant.len(), nex-fr);
	    println!("plen = {} {}",lin_params.len(),i);
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
		x[(i,j)] = x[(i,j)] - p_ant[i] * azimuth_noise[(i,j)];
	    }
	}
    }




    /// Applies the power function noise floor
    /// without affine offsets.
    pub fn apply_lp_noisefield(x:Arc<TwoDArray>,
			       mp_dict:Vec<Vec<ArrToArr>>,
			       burst_coords:&Vec<Vec<Arc<BurstEntry>>>,
			       split_indices:&MidPoint,
			       swath_bounds:&[Vec<SwathElem>],
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
		    for cur_burst in 0..num_burst {
			let (b_fa, b_la, _b_fr, _b_lr) = unpack_bound!(burst_coords[swath][cur_burst]);

		    
			// TODO: ensure that burst_coords are sorted by fa.

			let apply_subtract = |fa, la, fr, lr| {
			    let p_ant = LpApply::get_gapless_ant(fr, lr, mp_dict[swath][cur_burst].clone(),
							&splitinds[swath][cur_burst],
							&lin_params[swath]);
			};

			// First burst and swath.
			if cur_burst == 0 && swath_bounds[swath][0].fa < b_fa {
			    let (fa, la, fr, lr) = unpack_bound!(swath_bounds[swath][cur_burst]);
			    apply_subtract(fa, b_fa, fr, lr);
			    // Apply azimuth noise to the p_ant.
			}

			// Last subswath.
			if cur_burst == num_burst-1 && swath_bounds[swath].last().unwrap().fa > b_fa{
			    let (fa, la, fr, lr) = unpack_bound!(swath_bounds[swath].last().unwrap());
			    apply_subtract(fa, la, fr, lr);
			}

			// Iterate over every swath entry and insert in segments where it is fit.
			for e in 0..swath_bounds[swath].len() {
			    let (fa, la, fr, lr)  = unpack_bound!(swath_bounds[swath][e]);
			    
			    if fa > b_la || b_fa > la {continue;}
			    
			    let s_fa = b_fa.max(fa);
			    let s_la = b_la.min(la);

			    apply_subtract(s_fa, s_la, fr, lr);

			    // Check missing row
			    if e == num_burst-1 && la > s_la {
				/* appply */
			    }
			    else if e < num_burst-1 { // Apply between burst coordinates
				let (n_fa, n_la, mn_fr, n_lr) = unpack_bound!(burst_coords[swath][cur_burst+1]);
				if n_fa > b_la && n_fa <= la {
				    /* apply */
				}
			    }
			}
		    }
		}
	    }
	}
    }


    /*
def square_gapless_segmethod_denoise(x, y, mp_dict, slope_list, intercept_list, burst_coords, swath_bounds,  azimuth_noise, split_indices, sent_mode = 'EW', one_flag = False):
    """
    Interpolates radiation pattern to match the full extent of the swathbounds.
    """
    result = np.zeros(x.shape)
    N = NUMSUBSWATHS[sent_mode]
    
    if one_flag:
        # on the fly replace the function if flag.
        _get_gapless_ant = _get_one_ant
    else:
        _get_gapless_ant = _get_gapless_ant_

    ###swath_mask = np.zeros(x.shape, dtype=np.int)
    #check_mask = np.zeros(x.shape, dtype=np.int)
    for e in range(N):
        swath = SUBSWATH_NAMES[sent_mode][e]

        bursts = burst_coords[swath]
        al = [a for a, _, __, ___ in burst_coords[swath]]
        burst_inds = np.argsort(al) # sort according to azimuth (row)
        slope = slope_list[e]
        intercept = intercept_list[e]
        
        cur_burst = burst_inds[0]
        k = cur_burst
        swath_count = 0
        for fa, la, fr, lr in swath_bounds[e]:

            
            la += 1
            lr += 1
            #swath_mask[fa:la, fr:lr] += 1
            
            # cur_burst is the last used burst, it is prserved from 
            if swath_count == 0 and fa < bursts[cur_burst][0]:
                # swath bound is outside the ascribed burst
                b_fa, _b_la, _fr, _lr = bursts[cur_burst]

                p_ant = _get_gapless_ant(fr, lr, mp_dict, swath, cur_burst, split_indices, intercept, slope)
                assert(np.all(p_ant!=0))
                result[fa: b_fa, fr:lr] = x[fa: b_fa, fr:lr] - p_ant[np.newaxis,:] * azimuth_noise[fa: b_fa, fr:lr]
                #print('outside')
                #check_mask[fa: b_fa, fr:lr] += 1

            # This implies that there are no bursts that fit the last subswath
            if fa > bursts[burst_inds[-1]][0]:
                #print(fa,la,fr,lr,k,'trigger')
                cur_burst = burst_inds[-1]
                b_fa, _b_la, _fr, _lr = bursts[cur_burst]

                p_ant = _get_gapless_ant(fr, lr, mp_dict, swath, cur_burst, split_indices, intercept, slope)
                result[fa: la, fr:lr] = x[fa: la, fr:lr] - p_ant[np.newaxis,:] * azimuth_noise[fa: la, fr:lr]
                #check_mask[fa: la, fr:lr] += 1


            # iterate over every burst entry to see if the burst fits in this subswath split
            for q, cur_burst in enumerate(burst_inds):
                k = cur_burst
                # account for everything within the burst.
                b_fa, b_la, _b_fr, _b_lr = bursts[cur_burst]
                
                #b_la += 1
                if fa > b_la or b_fa > la: continue
                b_fa = max(b_fa, fa)
                b_la = min(b_la, la) # Asserts changes remain between subswath divisions

                ant = mp_dict[swath][cur_burst](np.arange(fr, lr)) # compute pattern for the initial one.
                # Get adjusted antenna denoising pattern
                p_ant = _get_gapless_ant(fr, lr, mp_dict, swath, cur_burst, split_indices, intercept, slope)
                
                result[b_fa:b_la, fr:lr] = x[b_fa:b_la, fr:lr] - p_ant[np.newaxis,:] * azimuth_noise[b_fa:b_la, fr:lr]
                #check_mask[b_fa:b_la, fr:lr] += 1
                assert(np.all(p_ant!=0))

                #print('q={} fa={} la={} b_fa={} b_la={} subswath = {}'.format(q, fa, la, b_fa, b_la, swath))


                # Fill in any missing rows
                if la > b_la and q == len(burst_inds)-1: # Last row, apply to the rest
                    result[b_la:la, fr:lr] = x[b_la:la, fr:lr] - p_ant[np.newaxis,:] * azimuth_noise[b_la:la, fr:lr]
                    #print('remain')
                    #check_mask[b_la:la, fr:lr] += 1

                    
                elif q < len(bursts)-1: # apply between burst coordinates.
                    n_fa, n_la, n_fr, n_lr = bursts[burst_inds[q+1]]

                    if n_fa > b_la and n_fa <= la: #there is a gap
                        #print('between q={} n_fa={} b_la={} lr={}'.format(q, n_fa, b_la, lr))
                        result[b_la:n_fa, fr:lr] = x[b_la:n_fa, fr:lr] - p_ant[np.newaxis,:] * azimuth_noise[b_la:n_fa, fr:lr]
                        #check_mask[b_la:n_fa, fr:lr] += 1
            swath_count+=1
    #assert(np.all(check_mask == swath_mask))
    return result
     */
    
    pub fn compute_affine() {
    }
}


