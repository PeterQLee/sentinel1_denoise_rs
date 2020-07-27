/// Top level module for preparing measurements for the linear programming implementation.
/// This is everything up-to the actual program estimation

use crate::parse::{BurstEntry, TimeRowLut, RawPattern, HyperParams};
use crate::read_from_archive::SentinelFormatId;
use lapack::*;
use std::sync::{Arc, Mutex};
use std::ops::{Index, IndexMut};

// FFT libs for fast convolution.
use rustfft::FFTplanner;
use rustfft::num_complex::Complex64;
use ndarray::{ArrayView2, Array2};
use std::thread;

pub type ArrToArr = Arc<Fn(&[usize]) -> Vec<f64> + Sync + Send> ;

pub struct TwoDArray{
    pub rows:usize,
    pub cols:usize,
    data:Vec<f64>
}

impl Index<(usize, usize)> for TwoDArray {
    type Output = f64;
    ///column major order
    fn index(&self, x:(usize, usize)) -> &f64{
	//&self.data[x.0 + x.1*self.rows]
	&self.data[x.1 + x.0*self.cols]
    }
}

impl IndexMut<(usize, usize)> for TwoDArray {
    fn index_mut(&mut self, x:(usize, usize)) -> &mut Self::Output {
	//&mut self.data[x.0 + x.1*self.rows]
	&mut self.data[x.1 + x.0*self.cols]
    }
}

impl TwoDArray {
    pub fn from_ndarray(ndarr:Array2<f64>) -> TwoDArray {
	let row = ndarr.nrows();
	let col = ndarr.ncols();
	
	if ndarr.is_standard_layout() {
	    return TwoDArray{
		rows:row,
		cols:col,
		data:ndarr.into_raw_vec()
	    };
	}
	else {
	    return TwoDArray{
		rows:row,
		cols:col,
		data:ndarr.reversed_axes().into_raw_vec()
	    };
	}
    }
}

macro_rules! get_num_subswath {
    ($sw:expr) => {
	match $sw.sentmode.as_str() {
	    "EW" => {5},
	    "IW" => {3},
	    _ => panic!("Incorrect sentmode")
	}
    }
}

#[macro_use]
macro_rules! argmax_row {
    ($sw:expr) => {
	$sw.iter().enumerate().fold((0, 0.0, false), |acc, x| {
	    if !acc.2 || *x.1 > acc.1 {return (x.0,*x.1,true);}
	    acc
	})
    };
}
#[macro_use]
macro_rules! argmin_row {
    ($sw:expr) => {
	$sw.iter().enumerate().fold((0, 0.0, false), |acc, x| {
	    if !acc.2 || *x.1 < acc.1 {return (x.0,*x.1,true);}
	    acc
	})
    };
}

macro_rules! unpack_bound {
    ($sw:expr) => {
	($sw.fa, $sw.la, $sw.fr, $sw.lr)
    }
}


/// Gets the clear midpoints.
/// eval: determines whether to do padded eval midpoints
/// Or just return the midpoints as is.
fn compute_midpoint(pat_angle:&Vec<f64>, pat_vals:&Vec<f64>, swath:usize, id:&SentinelFormatId,
		    eval:bool) -> Vec<f64>{ 
    	    
    match id.sentmode.as_str() {
	"EW" => {
	    match eval {
		true => {
		    let offsets = &[(-0.6,0.6), (-0.7, 1.0), (-0.8, 0.8), (-0.7, 0.8), (-0.7,0.9)];
		    if swath == 0 {
			let (a, _am, _fa) = argmax_row!(pat_vals);
			let (b, _bm, _fb) = argmax_row!(pat_vals[0..a/2]);
			let (c, _cm, _fc) = argmin_row!(pat_vals[b..a]);
			
			/*
			
			#midlist = [angles_list[b], angles_list[b+c], angles_list[a]]
			midlist = [angles_list[b], angles_list[b+c], angles_list[a]]*/
			return vec![pat_angle[0], pat_angle[b] + offsets[0].0,
				    pat_angle[b] + offsets[0].1, pat_angle[b+c] + offsets[0].0,
				    pat_angle[b+c] + offsets[0].1, pat_angle[a] + offsets[0].0,
				    pat_angle[a] + offsets[0].1, *pat_angle.last().unwrap()
			];
		    }
		    else {
			let (a, _am, _fa) = argmax_row!(pat_vals);
			return vec![pat_angle[0], pat_angle[a] + offsets[swath].0,
				    pat_angle[a] + offsets[swath].1, *pat_angle.last().unwrap()
			];
		    }

		}
		false => {
		    if swath == 0 {
			let (a, _am, _fa) = argmax_row!(pat_vals);
			let (b, _bm, _fb) = argmax_row!(pat_vals[0..a/2]);
			let (c, _cm, _fc) = argmin_row!(pat_vals[b..a]);

			return vec![(pat_angle[b]),
				    (pat_angle[b+c]),
				    (pat_angle[a])];
		    }
		    else {
			let (a, _am, _fa) = argmax_row!(pat_vals);
			return vec![(pat_angle[a])];
		    }
		}
	    }
	},
	"IW" => {
	    match eval {
		true => {
		    let offsets = &[(-1.0,1.0), (-1.0,1.0), (-1.0,1.0)];
		    let (a, _am, _fa) = argmax_row!(pat_vals);
		    return vec![pat_angle[0], pat_angle[a] + offsets[swath].0,
				pat_angle[a] + offsets[swath].1, *pat_angle.last().unwrap()];
		}
		false => {
		    let (a, _am, _fa) = argmax_row!(pat_vals);
		    return vec![(pat_angle[a])];
		}
	    }
	},
	_ => panic!("Incorrect sentmode")
    };
    
}

/// Builds a quadratic spline to map angles to the index/intensity lookup table.
/// http://nmbooks.eng.usf.edu/ebooks/05inp_spline/inp_05_spline_300_quadratic_example.html
fn quadratic_spline(x:&[f64], y:&[f64]) -> ArrToArr {
    // We assume that x is sorted
    let x_c:Vec<f64> = x.to_vec();
    let k = 2;

    let mut coef_list:Vec<(f64, f64, f64)> = Vec::new();
    
    for e in 0..x.len()-2 {
	let x0 = x[e];
	let x1 = x[e+1];
	let x2 = x[e+2];
	let y0 = y[e];
	let y1 = y[e+1];
	let y2 = y[e+2];
	
	let mut A = vec![0.0;(k+1)*(k+1)];
	let mut b = vec![0.0;k+1];
	
	for i in 0..k+1 {
	    b[i] = x0.powi(i as i32) * y0 + x1.powi(i as i32) * y1 + x2.powi(i as i32) * y2;
	    for j in 0..k+1 {
		A[j*(k+1) + i] = x0.powi((i+j) as i32) + x1.powi((i+j) as i32) + x2.powi((i+j) as i32);
	    }
	}
	// Solve linear system.
	let mut INFO:i32 = 0;
	let mut IPIV:Vec<i32> = vec![0;k+1];

	// dgesv (square system)
	unsafe {
	    dgesv((k+1) as i32,//num eqs
		  1 as i32,//num eqs
		  &mut A,
		  (k+1) as i32,//leading dim of A
		  &mut IPIV,//pivot matrix
		  &mut b,/////// right hand side
		  (k+1) as i32,//LDB
		  &mut INFO);
	}
	assert!(INFO == 0);
	coef_list.push((b[0], b[1], b[2]));
    }

    fn apply_quad(s:f64, coef:(f64, f64, f64)) -> f64 {
	coef.0 + s*coef.1 + s*s*coef.2
    }
    let len_x = x_c.len();
    Arc::new(move |cols:&[usize]| -> Vec<f64> {
	cols.iter().map(|s_| {
	    let s = (*s_) as f64;
	    if s <= x_c[0] {apply_quad(s, coef_list[0])}
	    else if s >= x_c[len_x-3] {apply_quad(s, *coef_list.last().unwrap())}
	    else {
		match x_c.binary_search_by(|sd| sd.partial_cmp(&s).expect("Nan encountered") ) {
		    Ok(t) => apply_quad(s,coef_list[t]),
		    Err(e) => apply_quad(s, coef_list[e-1])
		}
	    }

	}).collect()
    })
    
    
}

pub enum MidPoint {
    Est(Vec<Vec<Arc<Vec<f64>>>>),
    Test(Vec<Vec<Arc<Vec<(f64,f64)>>>>),
}
/// Returns a dictionary that provides the antenna values.
/// Along with indices of where to "split"
pub fn get_interpolation_pattern (buffer:&str,
				  burst_coords:&Vec<Vec<Arc<BurstEntry>>>,
				  time_row_lut:&TimeRowLut,
				  pat_raw:&RawPattern, id:&SentinelFormatId, eval:bool)
				  -> (Vec<Vec<ArrToArr>>, MidPoint) {
    let num_subswaths:usize = get_num_subswath!(id);

    let mut mp_dict:Vec<Vec<ArrToArr>> = Vec::new(); // vector of closures for pattern value derivation
    let mut test_indices:Vec<Vec<Arc<Vec<(f64,f64)>>>> = Vec::new();
    let mut est_indices:Vec<Vec<Arc<Vec<f64>>>> = Vec::new();

    for (swath, bu_co) in burst_coords.iter().enumerate() {
	mp_dict.push(Vec::new());
	match eval {
	    true => test_indices.push(Vec::new()),
	    false => est_indices.push(Vec::new())
	}
	for (n, burst) in bu_co.iter().enumerate() {

	    // Get raw pattern values.
	    let pat_angle = &pat_raw.angle[swath][n];
	    let pat_vals = &pat_raw.pattern[swath][n];


	    // get interp_model
	    let multi = time_row_lut.interp_angle_to_col();
	    
	    // take the midrow as the one to use for interpolation
	    let mid_angles = compute_midpoint(pat_angle, pat_vals, swath, id, eval);
	    
	    let row = (burst.fa+burst.la)/2;

	    let pat_cols = multi(&pat_angle, row);

	    let mp = quadratic_spline(&pat_cols, &pat_vals);

	    
	    mp_dict[swath].push(mp);

	    // Match the midpoint formula according to evaluation or estimation.
	    match eval {
		true => {
		    let t = multi(&mid_angles, row);
		    test_indices[swath].push(Arc::new(
					     t.iter().enumerate().step_by(2).map(|x| (*x.1,t[x.0+1])).collect()));
		}
		false => {est_indices[swath].push(Arc::new(multi(&mid_angles, row)));}
	    }
	}
    }
    let split_indices:MidPoint = match eval {
	true => MidPoint::Test(test_indices),
	false => MidPoint::Est(est_indices)
    };
    (mp_dict, split_indices)
    
}



fn compute_mean_slice(x:Arc<TwoDArray>, swath_bounds:BurstEntry)  -> Vec<f64> {
    let mut res = vec![0.0; swath_bounds.lr - swath_bounds.fr];
    for i in swath_bounds.fr..swath_bounds.lr {
	let mut cum_val = 0.0;
	let norm = (swath_bounds.la-swath_bounds.fa) as f64;
	let mut nonzero = true;
	for j in swath_bounds.fa..swath_bounds.la {
	    nonzero = nonzero & (x[(j,i)] != 0.0);
	    cum_val += x[(j,i)].max(0.0)/norm;
	}
	if nonzero {
	    res[i-swath_bounds.fr] = cum_val;
	}
    }
    res
}

/// Performs boxcar smoothing on the entry via FFT.
fn boxcar(sl:&[f64], hyper:&HyperParams) -> Vec<f64> {

    // get input data prepared
    let mut sl_r:Vec<Complex64> = sl.iter().map(|x| Complex64::new(*x/(sl.len() as f64), 0.0)).collect();
    let mut f_sl_r:Vec<Complex64> = vec![Complex64::new(0.0,0.0); sl.len()];
    
    // Construct box filter
    let mut box_f:Vec<Complex64> = vec![Complex64::new(0.0,0.0); sl.len()];
    let mut f_box_f:Vec<Complex64> = vec![Complex64::new(0.0,0.0); sl.len()];
    //let a = &mut box_f[0..hyper.box_l/2+1];
    box_f[0..hyper.box_l/2+1].iter_mut().for_each(|x| *x = Complex64::new(1.0/(hyper.box_l as f64), 0.0));

    box_f[sl.len()-hyper.box_l/2..sl.len()].iter_mut().for_each(|x| *x = Complex64::new(1.0/(hyper.box_l as f64),0.0));

    // Get the FFT stuff ready
    let mut f_planner = FFTplanner::new(false);
    let mut i_planner = FFTplanner::new(true);

    let fft = f_planner.plan_fft(sl.len());
    let ifft = i_planner.plan_fft(sl.len());

    // Apply fft
    fft.process(&mut box_f, &mut f_box_f);
    fft.process(&mut sl_r, &mut f_sl_r);

    //Apply convolution
    let mut f_res:Vec<Complex64> = f_box_f.iter().zip(f_sl_r.iter()).map(
	|x| x.0*x.1).collect();

    let mut res:Vec<Complex64> = vec![Complex64::new(0.0,0.0); sl.len()];
    // apply inverse fft
    ifft.process(&mut f_res, &mut res);
    
    // Finally convert values back to f64.
    res.iter().map(|x| x.re).collect()

    
}

type EstSegment = Vec<Vec<f64>>;
fn process_segment(x:Arc<TwoDArray>, burst_coords:Arc<BurstEntry>, swath:usize, ant:ArrToArr,
		   split_index:Arc<Vec<f64>>, o_value:f64, hyper:Arc<HyperParams> ) -> (EstSegment, EstSegment){
    let padding:usize = hyper.burst_padding;
    let num_subswaths = 5; // TODO: fix

    // Get mean value
    let mut sl = compute_mean_slice(x, BurstEntry{
	fa:burst_coords.fa + padding,
	la:burst_coords.la - padding + 1 ,
	fr:burst_coords.fr + padding,
	lr:burst_coords.lr - padding + 1});
    // ensure that values are above 0
    sl.iter_mut().for_each(|x| if *x < 0.0 { *x = 0.0;});

    // get zero mask to remove the low values that were convolved.
    let zero_mask = sl.iter().map(|x| *x <= 0.0);
    // get filtered slice
    let mut filt_sl = boxcar(&sl, &hyper);
    // rezero the appropriate entries
    // (this ensures we don't get weird drops in output that inhibit LP estimation)
    filt_sl.iter_mut().zip(zero_mask).for_each(|x| if x.1 {*x.0 = 0.0});

    let len_filt = filt_sl.len();
    // next, zero out the weird values on the edges
    if swath == 0 { // note, this is redundant, but I'm keeping it bc it was in the python code
	filt_sl[0..25].iter_mut().for_each(|x| *x = 0.0);
    }

    else if swath == num_subswaths - 1 {
	filt_sl[len_filt-70..].iter_mut().for_each(|x| *x = 0.0);
    }

    // Now find out where to amputate
    let mut mn = 10000;
    let mut mnfound = false;
    let mut mxfound = false;
    let mut mx = 0;
    for (e,x) in filt_sl.iter().enumerate() {
	if *x != 0.0 { mn = e; mnfound = true; break;}
    }
    for (e,x) in filt_sl.iter().enumerate() {
	if *x != 0.0 { mx = e; mxfound = true; }
    }
    if !mnfound  || !mxfound { mn = 0; mx = 50}

    filt_sl = filt_sl[mn+25..mx-25].to_vec();

    // Finally, return the segments divided by the bursts
    let mut ln_real_list:EstSegment = Vec::new();
    let mut ln_ant_list:EstSegment = Vec::new();
    
    let mut prev:usize;
    let mut nex:usize = 0;
    let fr = burst_coords.fr + padding + mn + 25;

    // PRocess log segment
    fn get_seg(prev:usize, nex:usize, filt_sl:&[f64], antvals:&[f64], o_value:f64) -> (Vec<f64>, Vec<f64>){
	let mut lreal = Vec::new();
	let mut lant = Vec::new();
	
	for i in prev..nex {
	    if filt_sl[i]-o_value <= 0.0 {continue;}
	    lreal.push((filt_sl[i]-o_value).ln());
	    lant.push(antvals[i].ln());
	}

	(lreal, lant)
    }
    let inds:Vec<usize> = (fr..fr+filt_sl.len()).collect();
    let ant_vals = ant(&inds);
    for i in 0..split_index.len() {
	prev = nex + 2*hyper.add_pad;
	nex = (split_index[i].round() as usize - fr) as usize - hyper.add_pad;
	
	let (real, ant) = get_seg(prev, nex, &filt_sl, &ant_vals, o_value);
	ln_real_list.push(real);
	ln_ant_list.push(ant);
    }
    // final segment.
    prev = nex + 2*hyper.add_pad;
    nex = filt_sl.len();
    let (real, ant) = get_seg(prev, nex, &filt_sl, &ant_vals, o_value);
    ln_real_list.push(real);
    ln_ant_list.push(ant);

    (ln_real_list, ln_ant_list)
    
    
}

/// Processes the segments in the splits/ bursts and
/// returns the slope/intercept parameters for each slope.
pub fn select_and_estimate_segments(x:Arc<TwoDArray>, mp_dict:Vec<Vec<ArrToArr>>,
				    burst_coords:&Vec<Vec<Arc<BurstEntry>>>,
				    split_indices:&MidPoint,
				    o_list:Vec<Vec<f64>>,
				    hyper:Arc<HyperParams>,
				    id:&SentinelFormatId) -> Vec<Vec<crate::est_lp::lin_params>>{
    let num_subswaths:usize = get_num_subswath!(id);
    let mut ret:Vec<Vec<crate::est_lp::lin_params>> = Vec::with_capacity(num_subswaths);
    match split_indices {
	MidPoint::Est(splitind) => {
	    for swath in 0..num_subswaths {
		let n_splits:usize = splitind[swath][0].len() + 1;
		let mut lreal:Arc<Mutex<EstSegment>> = Arc::new(Mutex::new(vec![Vec::new();n_splits]));
		let mut lant:Arc<Mutex<EstSegment>> = Arc::new(Mutex::new(vec![Vec::new();n_splits]));

		/* Prepare threads for gathering segments. */
		let mut thread_handles:Vec<_> = (2..burst_coords[swath].len()).map(
		    |e| { /* Clone reference counters.*/
			let x_ = x.clone();
			let hyper_ = hyper.clone();
			let burst = burst_coords[swath][e].clone();
			let o = o_list[swath][e];
			let _lreal =  lreal.clone();
			let _lant = lant.clone();
			let mp = mp_dict[swath][e].clone();
			let sind = splitind[swath][e].clone();
			thread::spawn(move || {
			    let (mut real, mut ant) = process_segment(x_, burst, swath, mp, sind, o, hyper_);
			    for i in (0..n_splits).rev() {
				{ // This will ensure both aspects are locked until mutation is finished.
				    let mut zd = _lreal.lock().unwrap();
				    zd[i].extend(real.pop().unwrap());
				    _lant.lock().unwrap()[i].extend(ant.pop().unwrap());
				}
			    }
			})
		    }).collect();

		/* Run the threads */
		let b = thread_handles.len();
		for i in 0..b {thread_handles.pop().unwrap().join();}

		/* Calculate the parameters via linear programming. */
		ret.push(
		    lreal.lock().unwrap().iter()
			.zip(lant.lock().unwrap().iter())
			.map(|x| crate::est_lp::solve_lp(&x.0, &x.1)).collect());
		
	    }
	}
	MidPoint::Test(_) => panic!("Wrong midpoints for estimation.")
    }
    return ret;
}
//https://cvxopt.org/userguide/coneprog.html#cvxopt.solvers.lp



fn determine_mino_value(base:Arc<TwoDArray>,
			fa:usize,
			la:usize,
			fr:usize,
			lr:usize) -> f64
{
    let mut minval:f64 = 99999999999.0;
    for j in fr..lr {
	let mut cur = 0.0;
	let mut count = 0;
	let mut nonzero = true;
	for i in fa..la {
	    nonzero = nonzero & (base[(i,j)] != 0.0);
	    cur += base[(i,j)];
	    count+=1;
	}
	if nonzero && cur > 0.0 {
	    minval = minval.min(cur/(count as f64));
	}
    }
    minval
    //minval.max(0.0)
}
pub fn compute_mino_list(base:Arc<TwoDArray>,
			 burst_coords:&Vec<Vec<Arc<BurstEntry>>>,
			 hyper:Arc<HyperParams>,
			 id:&SentinelFormatId) -> Vec<Vec<f64>> {
    let padding = hyper.burst_padding;
    let num_subswaths:usize = get_num_subswath!(id);
    let mut mino_list:Vec<Vec<f64>> = Vec::new();


    for swath in 0..num_subswaths {
	let tmp:Arc<Mutex<Vec<f64>>> = Arc::new(Mutex::new(vec![9999.0;burst_coords[swath].len()]));
	let mut handles = Vec::new();
	for c in 0..burst_coords[swath].len() {
	    let bt = burst_coords[swath][c].clone();
	    let tmp_ = tmp.clone();
	    let base__ = base.clone();

	    handles.push(thread::spawn(move || {
	    	let (fa_, la_, fr_, lr_ ) = unpack_bound!(bt);
	    	let fa = fa_ + padding;
	    	let la = la_ - padding + 1;
	    	let fr = fr_ + padding;
	    	let lr = lr_ - padding + 1;

	    	let mino = determine_mino_value(base__,
	    					fa,la,fr,lr);
	    	tmp_.lock().unwrap()[c] = mino;
	    }));

	}
	let b = handles.len();
	for _i in 0..b {handles.pop().unwrap().join().unwrap();}

	mino_list.push(Arc::try_unwrap(tmp).unwrap()
		       .into_inner().unwrap());
	
    }

    mino_list
    
}
/*
def square_yoff(x, y, burst_coords, mp_dict, sent_mode = 'EW'):
    """
    Selects o to be the minimum denoised value from using x-y standard denoising.
    Should have more success on arctic scenes since it is better calibrated there.
    """
    padding = 40
    ynoised = denoise._denoise_sim(np.sqrt(np.clip(x,0, x.max())), np.sqrt(np.clip(y, 0, y.max())))
    N = NUMSUBSWATHS[sent_mode]
    o_list = [ [] for i in range(N)]
    

    for e in range(N):
        swath = SUBSWATH_NAMES[sent_mode][e]
        c=0
        for fa_, la_, fr_, lr_ in burst_coords[swath]:
            la_+=1
            lr_+=1
            fa = fa_ + padding
            la = la_ - padding
            fr = fr_ + padding
            lr = lr_ - padding

            real_ = ynoised[fa:la,fr:lr] #?

            real_sum = np.sum(real_, axis=0)
            real_mean = np.zeros(real_.shape[1])

            rmask = np.all(real_ != 0, axis=0)

            #real_mean[rmask ] = real_sum[rmask] / np.sum(real_ != 0, axis=0)[rmask]
            real_mean[rmask] = np.mean(real_, axis=0)[rmask]
            #np.min(real_mean[real_mean>0])
            if np.all(~(real_mean>0)):
                s = 0
            else:
                s = np.min(real_mean[real_mean>0])

            o_list[e].append(s)

            c+=1
    return o_list
*/
