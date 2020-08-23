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

#include "scs_solve.h"
#include <string.h>
#include <stdio.h>

void copy_settings(ScsSettings *x,
		   lp_scs_settings *y) {
  x->normalize = y->normalize;
  x->scale = y->scale;
  x->rho_x = y->rho_x;
  x->max_iters = y->max_iters;
  x->eps = y->eps;
  x->cg_rate = y->cg_rate;
  x->verbose = y->verbose;
  
}
lin_params scs_solve_lp(unsigned int n,
			unsigned int m,
			double *A_,
			double *b_,
			double *c_,
			lp_scs_settings *settings) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));

  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  d->stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));

  k-> l = m; // set only to inequality constraints
  d->m = m;
  d->n = n; //number of variables.

  ScsMatrix *A = d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  scs_float *b = d->b = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *c = d->c = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *z = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_int i, j, r, rn, rm;

  
  A->i = (scs_int *)scs_calloc(n*m, sizeof(scs_int));
  A->p = (scs_int *)scs_calloc((n + 1), sizeof(scs_int));
  A->x = (scs_float *)scs_calloc(n*m, sizeof(scs_float));
  A->n = d->n;
  A->m = d->m;

  memcpy(A->x, A_, n*m*sizeof(double));
  memcpy(b, b_, m*sizeof(double));
  memcpy(c, c_, n*sizeof(double));

  // setup index arrays for the matrix
  for (int i=0; i < n+1; i++) {
    A->p[i] = i*m;
  }

  for (int i=0; i < m; i++) {
    A->i[i] = i;
    A->i[m+i] = i;
  }
  SCS(set_default_settings)(d);
  copy_settings(d->stgs, settings);

  // solve
  scs(d, k, sol, &info);
  // extract result
  lin_params result = {sol->x[0], sol->x[1]};

  // free data.
  SCS(free_data(d,k));
  SCS(free_sol(sol));
  
  return result;
}
