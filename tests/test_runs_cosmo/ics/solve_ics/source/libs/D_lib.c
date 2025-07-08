#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

double integrand(double omegam, double omegal, double a);

//=======================
// H_dot
//=======================
double H_D_lib(double a, double omegam, double omegar, double omegal)
{

  // Since D is independient of H0, doesn't matter the value of H0 used.
  
  return sqrt(omegam/pow3(a) + omegal);

}

//======================
// D
//======================
double D_a(double a, double omegam, double omegar, double omegal, double H0, double *norm_D)
{

  // Calculates D(a) as integral, normalized at redshift 0.
  // 
  // Comments:
  //           - Uses trapecium to make the integral.
  
  double d;  // growth factor to calculate

  double ap;  // a prime
  double delta = 1.0e-7;
  double integ;   // integral.
  
  // Calculates the normalization:
  //==============================
  // Makes the integral:
  integ = 0.0;
  for(ap=delta; ap<=1.0; ap+=delta)
    {
      //fprintf(stderr, "%lg\n", ap);
      integ += ( integrand(omegam, omegal, ap) +
		 integrand(omegam, omegal, ap+delta) ) * delta / 2.0;
    }
  
  // Multiplies by H(a):
  *norm_D = H_D_lib(1.0, omegam, omegar, omegal) * integ;

  // Calculates the unnormalized function:
  //======================================
  // Makes the integral:
  integ = 0.0;
  for(ap=delta; ap<=a; ap+=delta)
    integ += ( integrand(omegam, omegal, ap) +
	       integrand(omegam, omegal, ap+delta) ) * delta / 2.0;
  
  // Multiplies by H(a):
  d = H_D_lib(a, omegam, omegar, omegal) * integ;

  // Normalizes:
  //============
  // Here, all the H0 cancel.
  d /= *norm_D;

  // Put units in the norm:
  //=======================
  *norm_D /= pow2(H0);

  return d;

}


//======================
// integrand
//======================
double integrand(double omegam, double omegal, double a)
{

  // Calculates the integrand of the integral (1/a_dot^3) in a numerically
  // friendly way.
  // The H0 factor that should appear because of the expresion for a_dot
  // is not here because it will cancel when the normalization is made.

  // In the case of D_dot, the H0 factor must be but is taken into account
  // in other place.

  return pow(a / (omegam + omegal*pow3(a)), 1.5);


}

