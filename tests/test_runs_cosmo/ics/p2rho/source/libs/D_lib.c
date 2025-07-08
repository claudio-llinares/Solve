#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

double integrand(double omega0, double omegal, double a);

//=======================
// a_dot
//=======================
double a_dot(double a, double omega0, double omegar, double omegal)
{

  // Since D is independient of H0, doesn't matter the value of H0 used.
  
  double H0 = 1.0;

  //return unit.H0 * sqrt(par.omega0 * (1.0/a - 1.0) + 
  //                    par.omegal * (pow2(a) - 1.0) +
  //                    1.0);

  return H0 * sqrt(omega0 * (1.0/a - 1.0) + 
		   omegar * (1.0/pow2(a) - 1.0) +
		   omegal * (pow2(a) - 1.0) +
		   1.0);

}


//======================
// D
//======================
double D_1000(double a, double omega0, double omegar, double omegal)
{

  // Calculates D(a) as integral formula in thesis.
  // 
  // Comments:
  //           - Normalized to one at redshift 1000.
  //           - Uses trapecium to make the integral.
  
  double d;  // growth factor to calculate
  double c;  // normalization = D(z=1000)

  double ap;  // a prime
  double delta = 1.0e-7;
  double integ;   // integral.
  double a_1000 = 1.0/(1.0+1000.0);
  
  //double b1, b2;  // for the trapecium stuff.

  // Calculates the normalization:
  //==============================
  // Makes the integral:
  integ = 0.0;
  for(ap=delta; ap<=a_1000; ap+=delta)
    integ += ( integrand(omega0, omegal, ap) +
	       integrand(omega0, omegal, ap+delta) ) * delta / 2.0;
  
  //integ += 1.0/pow3(a_dot(ap, omega0, omegar, omegal)) * delta;
  
  // Multiplies by H(a):
  c = (a_dot(a_1000, omega0, omegar, omegal) / a_1000) * integ;

  // Calculates the unnormalized function:
  //======================================
  // Makes the integral:
  integ = 0.0;
  for(ap=delta; ap<=a; ap+=delta)
    integ += ( integrand(omega0, omegal, ap) +
	       integrand(omega0, omegal, ap+delta) ) * delta / 2.0;
  
  // Multiplies by H(a):
  d = (a_dot(a, omega0, omegar, omegal) / a) * integ;

  // Normalizes:
  //============
  d /= c;

  return d;

}

//======================
// D
//======================
double D_0(double a, double omega0, double omegar, double omegal)
{

  // Calculates D(a) as integral formula in thesis.
  // 
  // Comments:
  //           - Normalized to one at redshift 0.
  //           - Uses trapecium to make the integral.
  
  double d;  // growth factor to calculate
  double c;  // normalization = D(z=1000)

  double ap;  // a prime
  double delta = 1.0e-7;
  double integ;   // integral.
  
  // Calculates the normalization:
  //==============================
  // Makes the integral:
  integ = 0.0;
  for(ap=delta; ap<=1.0; ap+=delta)
    integ += ( integrand(omega0, omegal, ap) +
	       integrand(omega0, omegal, ap+delta) ) * delta / 2.0;
  //integ += integrand(omega0, omegal, ap) * delta;
  //integ += 1.0/pow3(a_dot(ap, omega0, omegar, omegal)) * delta;
  
  // Multiplies by H(a):
  c = (a_dot(1.0, omega0, omegar, omegal) / 1.0) * integ;

  // Calculates the unnormalized function:
  //======================================
  // Makes the integral:
  integ = 0.0;
  for(ap=delta; ap<=a; ap+=delta)
    integ += ( integrand(omega0, omegal, ap) +
	       integrand(omega0, omegal, ap+delta) ) * delta / 2.0;
  //integ += integrand(omega0, omegal, ap) * delta;
  //integ += 1.0/pow3(a_dot(ap, omega0, omegar, omegal)) * delta;
  
  // Multiplies by H(a):
  d = (a_dot(a, omega0, omegar, omegal) / a) * integ;

  // Normalizes:
  //============
  d /= c;

  return d;

}

//======================
// D
//======================
double D_a(double a, double a_norm, double omega0, double omegar, double omegal)
{

  // Calculates D(a) as integral formula in thesis.
  // 
  // Comments:
  //           - Normalized to one at redshift given by a.
  //           - Uses trapecium to make the integral.
  
  double d;  // growth factor to calculate
  double c;  // normalization = D(z=given by a_norm)

  double ap;  // a prime
  double delta = 1.0e-7;
  double integ;   // integral.
  
  //double b1, b2;  // for the trapecium stuff.

  // Calculates the normalization:
  //==============================
  // Makes the integral:
  integ = 0.0;
  for(ap=delta; ap<=a_norm; ap+=delta)
    integ += ( integrand(omega0, omegal, ap) +
	       integrand(omega0, omegal, ap+delta) ) * delta / 2.0;
  
  //integ += 1.0/pow3(a_dot(ap, omega0, omegar, omegal)) * delta;
  
  // Multiplies by H(a):
  c = (a_dot(a_norm, omega0, omegar, omegal) / a_norm) * integ;

  // Calculates the unnormalized function:
  //======================================
  // Makes the integral:
  integ = 0.0;
  for(ap=delta; ap<=a; ap+=delta)
    integ += ( integrand(omega0, omegal, ap) +
	       integrand(omega0, omegal, ap+delta) ) * delta / 2.0;
  
  // Multiplies by H(a):
  d = (a_dot(a, omega0, omegar, omegal) / a) * integ;

  // Normalizes:
  //============
  d /= c;

  return d;

}


//======================
// integrand
//======================
double integrand(double omega0, double omegal, double a)
{

  // Calculates the integrand of the integral (1/a_dot^3) in a numerically
  // friendly way.
  // The H0 factor that should appear because of the expresion for a_dot
  // is not here because it will cancel when the normalization is made.

  return pow(a / (omega0 + omegal*pow3(a)), 1.5);


}

