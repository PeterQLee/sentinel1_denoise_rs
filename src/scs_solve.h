#include "glbopts.h"
#include "amatrix.h"
#include "cones.h"
#include "linalg.h"
#include "scs.h"
#include "util.h"
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

lin_params scs_solve_lp(unsigned int n,
		    unsigned int m,
		    double *A,
		    double *b,
		    double *c);
