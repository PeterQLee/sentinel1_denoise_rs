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

///Entry point for denoising zip file from C api
/// Under failure, assume that all elements are null.
#include <stdarg.h>
#include <stdint.h>
#include <stdlib.h>
//#include <new>

/// Struct for holding results from Linear results
/// crosspol: pointer to crosspol in f64 (square)
/// copol: pointer to u16
/// k: pointer to linear scales
/// rows: # rows in image
/// cols: # cols in image
typedef struct {
  double *crosspol;
  uint16_t *copol;
  double *k;
  uintptr_t rows;
  uintptr_t cols;
  uintptr_t n_subswaths;
}LinearResult;

void destroy_linearresult(LinearResult *a) {
  free(a->crosspol);
  free(a->copol);
  free(a->k);
}
						  

/// struct for holding raw results.
/// crosspol: pointer to crosspol in f64 (square)
/// copol: pointer to u16
/// y: pointer to noise field.
/// rows: # rows in image
/// cols: # cols in image
typedef struct  {
  uint16_t *crosspol;
  uint16_t *copol;
  double *y;
  uintptr_t rows;
  uintptr_t cols;
} RawResult;

void destroy_rawresult(RawResult *a) {
  free(a->crosspol);
  free(a->copol);
  free(a->y);
}

/// Struct for holding results from Linear results
/// crosspol: pointer to crosspol in f64 (square)
/// copol: pointer to u16
/// k: pointer to noise field.
/// m: pointer to slope params
/// b: pointer to intercept params
/// rows: # rows in image
/// cols: # cols in image
/// n_subswaths: number of subswaths in image.
/// plen: number of elements in m.
typedef struct  {
  double *crosspol;
  uint16_t *copol;
  double *k;
  double *m;
  double *b;
  uintptr_t rows;
  uintptr_t cols;
  uintptr_t n_subswaths;
  uintptr_t plen;
} LPResult;

void destroy_lpresult(LPResult *a) {
  free(a->crosspol);
  free(a->copol);
  free(a->k);
  free(a->m);
  free(a->b);
}


/// Simple array structure for multilook processing
/// Assumes Row major order, contiguous.
/// data: pointer to data
/// rows: number of rows in image
/// cols: number of columns in image
typedef struct  {
  double *data;
  uintptr_t rows;
  uintptr_t cols;
}SimpleArray2D;

typedef struct {
  SimpleArray2D cross;
  SimpleArray2D co;
}TwoStruct;

typedef struct {
  SimpleArray2D cross;
  SimpleArray2D co;
  SimpleArray2D y;
}ThreeStruct;

void destroy_simplearray2d(SimpleArray2D *a) {
  free(a->data);
}



/// Applies the linear scaling method using custom user provided scales, k.
/// Returns the values back in square intensity units.
///
/// Parameters:
///
/// path: char*
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// pathlen: size_t
///
/// kp: *double
///     One dimensional array in 64-bit float that indicates the linear scaling parameters
///     to apply to each subswath. Length of array must equal five.
///     To apply the standard ESA noise removal
/// klen: size_t
///
///
/// Returns:
/// LinearResult
LinearResult linear_get_customscale_data(const uint8_t *path,
                                         uintptr_t pathlen,
                                         double *kp,
                                         uintptr_t klen);

/// Applies the lstsquares estimation method to retrieve scaling parameters, k,
/// rescales the noise floor, y, and subtracts it from the image, x.
/// Returns the values back in square intensity units.
/// Parameters:
///
/// path: char*
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// pathlen: size_t
/// configpath: char*
///     Path to the config ini file.
/// configpathlen: size_t
///
/// Returns:
/// Struct of LinearResult. Upon error, all entries will be null.
LinearResult linear_get_dualpol_data(const uint8_t *path,
                                     uintptr_t pathlen,
                                     const uint8_t *configpath,
                                     uintptr_t configpathlen);

/// Returns the original cross pol, co pol, and noise field from the archive.
/// Note that arrays are in linear units.
///
/// Parameters:
///
/// path: char *
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// pathlen: size_t
///
/// Returns:
/// RawResult
RawResult linear_get_raw_data(const uint8_t *path, uintptr_t pathlen);

///  Applies the power function noise floor obtained from linear programming,
///   with parameters given by the user.
/// Returns the values back in square intensity units.
///
/// Parameters:
///
/// archpath: str
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// m: list
///     List of slope parameters
/// b: list
///     List of intercept parameters.
/// config_path: str or None
///     Optional path to config file. If None (or non-string) will use default configuration.
LPResult lp_get_customscale_data(const uint8_t *path,
                                 uintptr_t pathlen,
                                 double *m,
                                 double *b,
                                 uintptr_t plen,
                                 const uint8_t *configpath,
                                 uintptr_t configpathlen);

/// Applies the linear programming method to restimate a noise floor based
/// on the characteristics of the original image and the
/// Returns the values back in square intensity units.
///
/// Parameters:
///
/// path: char*
///     Path to the zip or directory unpacked from the Sentinel-1 zip archive
/// pathlen: size_t
///
/// lstsq_rescale: unsigned char
///     Indicate whether you want to apply the least squares method from
///     linear_get_dualpol_data to get the baseline minimum offset values for
///     the method. Ignored if product type is IW
///     >0 for applying the method
///     =0 to just use the default ESA noise floor for this.
/// configpath: char *
/// configpathlen: size_t
///     set to 0 if there is no config file and you want to use default settings.
///
/// Returns:
/// LPResult
LPResult lp_get_dualpol_data(const uint8_t *path,
                             uintptr_t pathlen,
                             uint8_t lstsq_rescale,
                             const uint8_t *configpath,
                             uintptr_t configpathlen);

/// Applies multilooking to the input image, sets negative values to 0, and
/// square roots the output values.
///
/// Parameters:
/// x: Input array for multilooking (square units) - SimpleArray2D
/// row_factor: integer amount to mean reduce row by.
/// col_factor: integer amount to mean reduce col by.
/// num_cores: number of multithreading cores to use.
///
/// Output:
/// x : multilooked image (linear units) - SimpleArray2D
SimpleArray2D post_multilook_and_floor(SimpleArray2D *x,
                                       uintptr_t row_factor,
                                       uintptr_t col_factor,
                                       uintptr_t num_cores);



/// All of the below are functions that conver the result types to SimpleArray2D types.
double* copy_u16_to_f64 (uint16_t *d, uintptr_t rows, uintptr_t cols) {
  double *o =(double*) malloc(sizeof(double)*rows*cols);
  for (uintptr_t i=0;i<rows;i++) {
    for (uintptr_t j=0;j<cols;j++) {
      o[i*cols+j] = (double) d[i*cols+j];
    }
  }
  return o;
}

TwoStruct linear_to_simplearray(LinearResult *x) {
  SimpleArray2D crosspol;

  crosspol.data = x->crosspol;
  crosspol.rows = x->rows;
  crosspol.cols = x->cols;
  
  SimpleArray2D copol;
  double *d = copy_u16_to_f64(x->copol, x->rows, x->cols);
  copol.data = d;
  copol.rows = x->rows;
  copol.cols = x->cols;

  free(x->copol);
  free(x->k);
  
  x->crosspol = NULL;
  x->copol = NULL;
  x->k = NULL;

  TwoStruct a = {crosspol, copol};
  return a;
}

TwoStruct lpresult_to_simplearray(LPResult *x) {
  SimpleArray2D crosspol;

  crosspol.data = x->crosspol;
  crosspol.rows = x->rows;
  crosspol.cols = x->cols;

  SimpleArray2D copol;
  double *d = copy_u16_to_f64(x->copol, x->rows, x->cols);
  copol.data = d;
  copol.rows = x->rows;
  copol.cols = x->cols;

  free(x->copol);
  free(x->k);
  free(x->m);
  free(x->b);

  x->crosspol = NULL;
  x->copol = NULL;
  x->k = NULL;
  x->m = NULL;
  x->b = NULL;
  TwoStruct a = {crosspol, copol};
  return a;
}



ThreeStruct rawresult_to_simplearray(RawResult *x) {
  SimpleArray2D crosspol;


  double *o = copy_u16_to_f64(x->crosspol, x->rows, x->cols);
  crosspol.data = o;
  crosspol.rows = x->rows;
  crosspol.cols = x->cols;

  SimpleArray2D copol;
  double *d = copy_u16_to_f64(x->copol, x->rows, x->cols);
  copol.data = d;
  copol.rows = x->rows;
  copol.cols = x->cols;

  SimpleArray2D y;
  y.data = x->y;
  y.rows = x->rows;
  y.cols = x->cols;

  free(x->crosspol);
  free(x->copol);
  x->crosspol = NULL;
  x->copol = NULL;
  x->y = NULL;

  ThreeStruct a = {crosspol, copol, y};
  return a;
}
