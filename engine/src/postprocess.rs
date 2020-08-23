//! Routines for post processing the scene (e.g. multilooking, handling negative)
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

use std::thread;
use std::sync::{Arc};
use crate::prep_lp::TwoDArray;
/// Applies multilook processing and flooring negative input values
///
/// Parameters:
/// * `x`- array
/// * `row_factor`- amount to reduce row by: newrow = floor(oldrow/row_factor)
/// * `col_factor`- amount to reduce col by: newcol = floor(oldrow/col_factor)
/// * `num_cores`- number of cpus to parallize over.
pub fn multilook_and_floor(x:Arc<TwoDArray>, row_factor:usize, col_factor:usize, num_cores:usize) -> TwoDArray {

    let newrow = x.rows/row_factor;
    let newcol = x.cols/col_factor;

    // TODO: Check that the output array is of valid dimensions

    // Allocate output array.
    let output:Arc<TwoDArray> = Arc::new(TwoDArray::zeros(newrow, newcol));


    let subsum = move |r:&Arc<TwoDArray>,i:usize,j:usize| -> f64 {
	let mut total = 0.0;
        let a_end = row_factor;
        let b_end = col_factor;
        let norm = 1.0/((a_end*b_end) as f64);
        for a in 0..a_end {
            for b in 0..b_end {
                total += r[(row_factor*i+a,col_factor*j+b)] * norm;
            }
        }
	return total;
    };
    
    // Since the array is row major, we divide jobs by number of rows
    let mut thread_handles:Vec<_> = (0..num_cores).map( {
	|t| {
	    let start_dest_row = t*newrow/num_cores;
	    let end_dest_row = if t < num_cores -1 {(t+1)*newrow/num_cores} else {newrow};
	    let x_ = x.clone();
	    let mut o_ = output.clone();
	    thread::spawn(move ||{
		let op = unsafe{Arc::get_mut_unchecked(&mut o_)};
		for i in start_dest_row..end_dest_row {
		    for j in 0..newcol {
			op[(i,j)] = subsum(&x_, i, j).max(0.0).sqrt();
		    }
		}
	    })
	}
    }).collect();

    for _i in 0..num_cores {thread_handles.pop().unwrap().join().unwrap();}

    Arc::try_unwrap(output).unwrap()
}

