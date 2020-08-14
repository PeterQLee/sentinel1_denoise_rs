//! Estimation method specific to linear programming method
use std::sync::Arc;
use crate::parse::HyperParams;
    
# [ repr ( C ) ] # [ derive ( Debug , Copy , Clone ) ]
pub struct lin_params { pub m : f64 , pub b : f64 , }

#[repr(C)]
#[derive ( Debug , Copy , Clone )]
pub struct lp_scs_settings {
    pub normalize : ::std::os::raw::c_int,
    pub scale : f64 ,
    pub rho_x : f64 ,
    pub max_iters : ::std::os::raw::c_int,
    pub eps : f64 ,
    pub cg_rate : f64 ,
    pub verbose : ::std::os::raw::c_int , } 

extern "C" {
    pub fn scs_solve_lp ( n : ::std::os::raw::c_uint,
			  m : ::std::os::raw::c_uint ,
			  A : * mut f64 ,
			  b : * mut f64 ,
			  c : * mut f64,
			  settings : * mut lp_scs_settings) -> lin_params ; }

#[allow(non_snake_case)]
/// Sets up linear program to be solved in SCS.
/// * `lreal`- slice of log measured values
/// * `lant`- slice of log antenna values
/// * `hyper`- hyper params to send to SCS.
pub fn solve_lp(lreal:&[f64], lant:&[f64], hyper:Arc<HyperParams>) -> lin_params {
    let m = lreal.len() + 2; // num constraints
    let n  = 2; // num variables

    let mut A:Vec<f64> = vec![0.0;m*n];
    let mut b:Vec<f64> = vec![0.0;m];
    let mut c:Vec<f64> = vec![0.0;n];

    // Prepare constraints
    // Note A is in column major order.
    for i in 0..lreal.len() {
	A[i] = lant[i];
	A[m + i] = 1.0;
	b[i] = lreal[i];
    }

    A[lreal.len()] = 1.0;
    b[lreal.len()] = -0.75;
    A[lreal.len()+1] = -1.0;
    b[lreal.len()+1] = 1.25;

    // compute gamma (use low_lp)
    let (minind, minval, __) = argmin_row!(lant);
    let (maxind, maxval, __) = argmax_row!(lant);
    let vs = ((lreal[maxind].exp() + lreal[minind].exp())/2.0).ln();
    let percent = (vs - lreal[minind])/(lreal[maxind]-lreal[minind]);
    let gamma = (1.0-percent)*(maxval - minval) + minval;

    c[0] = -gamma;
    c[1] = -1.0;

    let mut scs_variables = hyper.get_scs();

    // Solve the problem
    unsafe {scs_solve_lp(n as std::os::raw::c_uint,
			 m as std::os::raw::c_uint,
			 A.as_mut_ptr(),
			 b.as_mut_ptr(),
			 c.as_mut_ptr(),
			 &mut scs_variables
    )}
    
}
