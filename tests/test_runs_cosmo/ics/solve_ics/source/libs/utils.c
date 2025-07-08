#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solve.h"
#include "utils.h"

//=======================
// modulus
//=======================
double modulus(double x, double y)
{

  return x-(floor(x/y)*y);

}

//=======================
// fopen_check
//=======================
FILE *fopen_check(char *file, char *type)
{
 
  FILE *fp;

  if( (fp = fopen(file, type)) == NULL)
    {
      printf("Problems opening the file %s.  Aborting.\n", file);
      exit(0);
    }

  return fp;

}


//=========================
// euclides2
//=========================
double euclides2(double x, double y, double z)
{

  return pow2(x) + pow2(y) + pow2(z);


}

//=========================
// euclides
//=========================
double euclides(double x, double y, double z)
{

  return sqrt(pow2(x) + pow2(y) + pow2(z));


}

//=========================
// euclides2
//=========================
double dist(double x, double y, double z)
{

  return sqrt(pow2(x) + pow2(y) + pow2(z));


}

//=========================
// change_system
//=========================
void change_system(double vec[3][3], double *x, double *y, double *z)
{

  // Changes coordinates (x,y,z) to a system defined by vec.
  //
  // - This routine overwrites (x,y,z).
  //
  
  double x_new, y_new, z_new;

  x_new = vec[0][0]*(*x) + vec[1][0]*(*y) + vec[2][0]*(*z);
  y_new = vec[0][1]*(*x) + vec[1][1]*(*y) + vec[2][1]*(*z);
  z_new = vec[0][2]*(*x) + vec[1][2]*(*y) + vec[2][2]*(*z);
  
  // Prepares output
  *x = x_new;
  *y = y_new;
  *z = z_new;

  
}

//========================
// sort_3_values
//========================
void sort_3_values(double *a)
{
  // Sort the three first elements of a in decreasing order.
  // 
  // The result is given in a.

  int i, j, k;
  double p;

  for (i=0;i<3;i++) 
    {
      k = i;
      p = a[i];
      for (j=i+1;j<3;j++)
	{
	  if (a[j] >= p) 
	    {
	      k = j;
	      p = a[j];
	    }
	}
      if (k != i) 
	{
	  a[k] = a[i];
	  a[i] = p;
	}
    }
  
}

//=========================
// indexes_1d_to_3d
//=========================
void index_1d_to_3d(double q, int *i, int *j, int *k, int n)
{
  // Determines 3D indexes from a 1D description of a 3D array.
  // Arguments:
  //             - n = number of nodes per dimenssion
  //

  float u;

  // Get k:
  *k = (int)floor(q/pow2((double)n));
  
  // Get j:
  u = fmod(q, pow2((double)n));
  *j = (int)floor(u/(double)n);
  
  // Get i:
  *i = (int)fmod(u, (double)n);
    
}
