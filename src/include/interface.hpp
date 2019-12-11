typedef struct {
  float *crosspol;
  float *copol;
  int rows;
  int cols;
} OutArr;

extern "C" OutArr denoise_zip (char *path, int pathlen);
