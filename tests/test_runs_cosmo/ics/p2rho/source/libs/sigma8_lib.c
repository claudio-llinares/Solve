#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow6(x) ((x)*(x)*(x)*(x)*(x)*(x))

double f(double phase, double p, double r);


double sigma8(double *phase, double *p, int n_points, double r)
{
  // Calculate sigma_r for a given r and a given power spectrum
  // making the integral.
  //
  // Arguments:
  //          - phase, p = power spectrum
  //          - n_points = number of points in these arrays.
  //          - r = radius you want for sigma (standard = 8).
  // Comments:
  //          - The size of teh interval is not uniform.
  //          - Assumes radius is given in Mpc/h (user should make the 
  //            conversion before).
  //          - The power spectrum should be given in (h/Mpc)^3.  In fact, the 
  //            units of the power spectrum should be the units of r.
  //          - Formula taken from Tegmark 2002.

  double aa, delta_a;
  double res;  // result
  int i;

  double pi;

  pi = 4.0*atan(1.0);
  
  n_points--;  // we need one more point at the end

  if(n_points/3.0 == (double)(n_points/3))
    n_points--;
  
  // Makes the integral:
  //====================
  aa = 0.0;
  res = 0.0;
  for(i=0; i<=n_points-1; i++)
    {
      delta_a = phase[i+1] - phase[i];
      res += ( (f(phase[i+1], p[i+1], r) + f(phase[i], p[i], r)) /
	       2.0 ) * (delta_a);
      aa += delta_a;
      //fprintf(stderr, "%lf\n", res, phase[i], );
    }
  
  res *= (4.0*pi) / pow3(2.0*pi);
  res = sqrt(res);

  return res;

}

//========================
// f
//========================
double f(double phase, double p, double r)
{
  // Integrand in the definition of sigma8.
    
  double window;
  double kr;
  
  kr = phase * r;

  window = 1/phase * 9.0/pow6(kr) * pow2(sin(kr) - kr*cos(kr));

  return p * pow3(phase) * window;
  

}


