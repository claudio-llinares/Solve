#include "D_lib.h"
#include <math.h>
#include <stdio.h>

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))

// Routines to calculate D_dot.
// You should link D_lib to use this.

//=========================
// eta
//=========================
double eta(double a, double omegam, double omegar, double omegal)
{

  // Warning.  This eta is not the one used by Bertsinger.
  //           (it is that one squared times some a factor).
  //
  
  return omegam * (1.0/pow3(a) - 1.0/pow2(a)) + 
    omegar * (1.0/pow4(a) - 1.0/pow2(a)) +
    omegal * (1.0 - 1.0/pow2(a)) +
    1.0/pow2(a);
  
}

//=========================
// eta_prime
//=========================
double eta_prime(double a, double omegam, double omegar, double omegal)
{

  // Warning.  This eta is not the one used by Bertsinger.
  //           (it is that one squared times some a factor).
  //
  
  return omegam * (1.0/a - 1.0) + 
    omegar * (1.0/pow2(a) - 1.0) +
    omegal * (pow2(a) - 1.0) +
    1.0;

}

//===========================
// eta_3
//===========================
double eta_3(double a, double omegam, double omegar, double omegal)
{

  // Warning.  This eta is not the one used by Bertsinger.
  //           (it is that one squared times some a factor).
  //
  
  return 1.5*omegam / pow3(a) + 
    2.0*omegar / pow4(a) +
    (1.0 - omegam - omegar - omegal) / pow2(a);

}

//================================
//  D_dot_nonorm
//================================
double D_dot_nonorm(double a, double omegam, double omegar, double omegal, double H0)
{

  // WARNING:  Omega radiation is not coded in the calculation of D.
  // Calculates D_dot.
  //
  // Arguments:
  //             - expansion factor were you want the thing.
  //             - cosmology.
  //
  double delta = 1.0e-7;
  double integ;   // integral.
  double ddot;
  double ap;

  // Makes the integral:
  integ = 0.0;
  for(ap=delta; ap<=a; ap+=delta)
    integ += ( integrand(omegam, omegal, ap) +
	       integrand(omegam, omegal, ap+delta) ) * delta / 2.0;

  ddot = -eta_3(a, omegam, omegar, omegal) * integ + 
    1.0/(sqrt(eta(a, omegam, omegar, omegal))*pow2(a));
  ddot /= H0;

  return ddot;

}

//================================
//  D_dot_nonorm
//================================
double D_dot_1000(double a, double omegam, double omegar, double omegal, double H0)
{

  // WARNING:  Omega radiation is not coded in the calculation of D.
  // Calculates D_dot assuming D normalize to 1 at redshift 1000.
  //
  // Arguments:
  //             - expansion factor were you want the thing.
  //             - cosmology.
  //
  double delta = 1.0e-7;
  double integ;   // integral.
  double ddot;
  double a_1000 = 1.0/(1.0+1000.0);
  double ap;

  double c; // normalization.


  // Calculates the normalization (D_nonorm(z=1000)):
  //=================================================
  // Makes the integral:
  // integrand = 1/(eta_prime)^(3/2)  
  integ = 0.0;
  for(ap=delta; ap<=a_1000; ap+=delta)
    integ += ( integrand(omegam, omegal, ap) +
	       integrand(omegam, omegal, ap+delta) ) * delta / 2.0;

  // Multiplies by H(a) (the H0 is needed explicitly):
  //c = H0 * (a_dot_D_lib(a_1000, omegam, omegar, omegal) / a_1000) * integ;
  c = sqrt(eta(a_1000, omegam, omegar, omegal))/pow2(H0) * integ;

  // Calculates the unnormalized function:
  //======================================
  ddot = D_dot_nonorm(a, omegam, omegar, omegal, H0);

  // Normalizes:
  //============
  ddot /= c;

  return ddot;

}

//================================
//  D_dot_nonorm
//================================
double D_dot_a(double a, double a_norm, double omegam, double omegar, double omegal, double H0)
{

  // WARNING:  Omega radiation is not coded in the calculation of D.
  // Calculates D_dot assuming D normalize to 1 at redshift 1000.
  //
  // Arguments:
  //             - expansion factor were you want the thing.
  //             - cosmology.
  //
  double delta = 1.0e-7;
  double integ;   // integral.
  double ddot;
  double ap;

  double c; // normalization.


  // Calculates the normalization (D_nonorm(z=1000)):
  //=================================================
  // Makes the integral:
  // integra = 1/(eta_prime)^(3/2)  
  integ = 0.0;
  for(ap=delta; ap<=a_norm; ap+=delta)
    integ += ( integrand(omegam, omegal, ap) +
	       integrand(omegam, omegal, ap+delta) ) * delta / 2.0;

  // Multiplies by H(a) (the H0 is needed explicitly):
  c = sqrt(eta(a_norm, omegam, omegar, omegal))/pow2(H0) * integ;

  // Calculates the unnormalized function:
  //======================================
  ddot = D_dot_nonorm(a, omegam, omegar, omegal, H0);

  // Normalizes:
  //============
  ddot /= c;

  return ddot;

}
