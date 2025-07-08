#include <stdio.h>
#include <stdlib.h>

#define pow2(x) ((x)*(x))

void cuad_min(double *x, double *y, int n, double *a, double *b)
{

  // Fits a line to n data points given in x, y.

  int i;
  double norm_x, norm_y, norm_x2, norm_xy;

  // Calculate norms:
  //=================
  norm_y = 0.0;
  norm_x = 0.0;
  norm_x2 = 0.0;
  norm_xy = 0.0;
  for(i=0; i<n; i++)
    {
      norm_y += y[i];
      norm_x += x[i];
      norm_x2 += pow2(x[i]);
      norm_xy += x[i] * y[i];
    }

  // Calculate the parameters:
  //==========================
  *a = (norm_y * norm_x2 - norm_xy * norm_x) / 
    (n * norm_x2 - pow2(norm_x));

  *b = (norm_x * norm_y - n * norm_xy) / 
    (pow2(norm_x) - n * norm_x2);

  // Output a test:
  //===============
  FILE *fp;
  if( (fp = fopen("test_cuad_min.dat", "w")) == NULL)
    {
      fprintf(stderr, "Problems opening the output file for cuad_min test\n");
      exit(0);
    }
  for(i=0; i<n; i++)
    fprintf(fp, "%lf  %lf\n", x[i], y[i]);
  
  fclose(fp);
  
}
