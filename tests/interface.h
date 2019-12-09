typedef struct {
  double data;
  int rows;
  int cols;
} OutArr;

extern OutArr denoise_zip (char *path, int pathlen);
