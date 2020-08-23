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

#include "scs/glbopts.h"
#include "scs/cones.h"
#include "scs/linalg.h"
#include "scs/scs.h"
#include "scs/util.h"
#include <stdlib.h>

/* Output results */
typedef struct {
  double m;
  double b;
} lin_params;

/* this struct defines the data matrix A */
struct SCS_A_DATA_MATRIX {
  /* A is supplied in column compressed format */
  scs_float *x; /* A values, size: NNZ A */
  scs_int *i;   /* A row index, size: NNZ A */
  scs_int *p;   /* A column pointer, size: n+1 */
  scs_int m, n; /* m rows, n cols */
};

typedef struct {
  int normalize;
  double scale;
  double rho_x;
  int max_iters;
  double eps;
  double cg_rate;
  int verbose;
}lp_scs_settings;
lin_params scs_solve_lp(unsigned int n,
			unsigned int m,
			double *A,
			double *b,
			double *c,
			lp_scs_settings *settings
			);
