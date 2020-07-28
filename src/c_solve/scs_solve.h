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
