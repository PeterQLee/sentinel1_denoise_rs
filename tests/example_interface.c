#include "interface.h"
#include <string.h>

int main() {
  char *path = "/mnt/D2/Data/Sentinel/beaufortredown/S1A_EW_GRDM_1SDH_20180902T164932_20180902T165032_023522_028FAA_5A8B.zip";
  OutArr res = denoise_zip(path, strlen(path));
}
