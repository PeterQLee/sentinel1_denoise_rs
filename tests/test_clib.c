#include "../include/interface.h"
#include <string.h>
int main () {
  const char * path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
  size_t pathlen = strlen(path);
  LinearResult res = linear_get_dualpol_data(path,
			  pathlen,
			  "",
			  0);

  destroy_linearresult(&res);

  double *k = malloc(sizeof(double)*5);
  for (int i=0;i<5;i++) k[i] = 1.0;

  

  LinearResult res1 = linear_get_customscale_data(path,
						 pathlen,
						 k,
						 5);
  destroy_linearresult(&res1);

  RawResult res2 = linear_get_raw_data(path,
				       pathlen);

  destroy_rawresult(&res2);


  LPResult res3 = lp_get_dualpol_data(path,
					  pathlen,
					  1,
					  "",
					  0
					  );


  double *m = malloc(sizeof(double)*res3.plen);
  double *b = malloc(sizeof(double)*res3.plen);
  for (int i=0;i<res3.plen;i++) m[i] = res3.m[i];
  for (int i=0;i<res3.plen;i++) b[i] = res3.b[i];


  LPResult res4 = lp_get_customscale_data(path,
					  pathlen,
					  m,
					  b,
					  res3.plen,
					  "",
					  0
					  );
						 
						 
  destroy_lpresult(&res3);

  destroy_lpresult(&res4);

  
}
