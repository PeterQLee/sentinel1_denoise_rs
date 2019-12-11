typedef struct {
  float *crosspol;
  float *copol;
  int rows;
  int cols;
} OutArr;

extern OutArr denoise_zip (char *path, int pathlen);
