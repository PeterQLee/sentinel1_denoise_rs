// Linear algebra

# [ repr ( C ) ] # [ derive ( Debug , Copy , Clone ) ]
pub struct lin_params { pub m : f64 , pub b : f64 , }

extern "C" {
    pub fn scs_solve_lp ( n : :: std :: os :: raw :: c_uint , m : :: std :: os :: raw :: c_uint , A : * mut f64 , b : * mut f64 , c : * mut f64 ) -> lin_params ; }

#[allow(non_snake_case)]
pub fn solve_lp(lreal:&[f64], lant:&[f64]) -> lin_params {
    let m = lreal.len() + 2; // num constraints
    let n  = 2; // num variables

    let mut A:Vec<f64> = vec![0.0;m*n];
    let mut b:Vec<f64> = vec![0.0;m];
    let mut c:Vec<f64> = vec![0.0;n];

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

    c[0] = gamma;
    c[1] = 1.0;

    unsafe {scs_solve_lp(n as std::os::raw::c_uint, m as std::os::raw::c_uint, A.as_mut_ptr(), b.as_mut_ptr(), c.as_mut_ptr())}
    
}
