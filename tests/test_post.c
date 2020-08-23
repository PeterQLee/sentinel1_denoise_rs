#include "../include/interface.h"
#include <string.h>
int main () {
  const char * path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
  size_t pathlen = strlen(path);
  LinearResult res = linear_get_dualpol_data(path,
			  pathlen,
			  "",
			  0);


  TwoStruct t = linear_to_simplearray(&res);
  
  SimpleArray2D x = post_multilook_and_floor(&(t.co),
					     16,16,16);
  
  destroy_simplearray2d(&x);
  destroy_simplearray2d(&(t.cross));
  destroy_simplearray2d(&(t.co));

  
}
