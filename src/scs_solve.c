#include "scs_solve.h"
#include <string.h>
lin_params scs_solve_lp(unsigned int n,
		    unsigned int m,
		    double *A_,
		    double *b_,
		    double *c_) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  d->stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));

  k-> l = m; // set only to inequality constraints
  d->m = m;
  d->n = n;

  ScsMatrix *A = d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  scs_float *b = d->b = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_float *c = d->c = (scs_float *)scs_calloc(n, sizeof(scs_float));
  scs_float *z = (scs_float *)scs_calloc(m, sizeof(scs_float));
  scs_int i, j, r, rn, rm;

  
  A->i = (scs_int *)scs_calloc(n*m, sizeof(scs_int));
  A->p = (scs_int *)scs_calloc((n + 1), sizeof(scs_int));
  A->x = (scs_float *)scs_calloc(n*m, sizeof(scs_float));

  memcpy(A->x, A_, n*m);
  memcpy(b, b_, m);
  memcpy(c, c_, n);

  // setup index arrays for the matrix
  A->p[0] = 0;
  A->p[1] = m;
  A->p[2] = m*2;

  for (int i=0; i < m; i++) {
    A->i[i] = i;
    A->i[m+i] = i;
  }
  SCS(set_default_settings)(d);
  
  scs(d, k, sol, &info);
  return lin_params{sol->x[0], sol->x[1]};


}
