#include <stdio.h>
#include <math.h>

#include <gsl/gsl_integration.h>

#include "../solve.h"
#include "cosmology.h"

struct params_int
{
  struct params *par;
};


//=======================
// hubble
//=======================
double hubble(struct params *par, struct units *unit, double a)
{

  return unit->H0 * sqrt(par->omegam/pow3(a) + par->omegal);

}

//=======================
// hubble2
//=======================
double hubble2(struct params *par, struct units *unit, double a)
{

  return pow2(unit->H0) * (par->omegam/pow3(a) + par->omegal);

}

//=======================
// hubble_dot
//=======================
double hubble_dot(struct params *par, struct units *unit, double a)
{

  //return -pow2(unit->H0)/2.0 * ( 3.0*par->omegam/pow3(a) );
  return -3.0/2.0 * par->omegam/pow4(a) / sqrt(par->omegam/pow3(a) + par->omegal);

}

//=======================
// a_dot
//=======================
double a_dot(struct params *par, struct units *unit, double a)
{

  return unit->H0 * sqrt(par->omegam/a +
			 par->omegal * pow2(a));
  
}

//=========================
// integrate_chi_of_a
//=========================
double integrate_chi_of_a(struct params *par, double a)
{

  gsl_integration_workspace *int_work;
  gsl_function int_f;
  double result, error;

  struct params_int par_int;

  int status;

  int_work = gsl_integration_workspace_alloc(8000);
  int_f.function = integrand_chi;
  int_f.params = &par_int;
  par_int.par = par;

  //status = gsl_integration_qags(&int_f, 1.0e-3, par->l_max, 0, 1e-5, 8000, int_work, &result, &error);
  status = gsl_integration_qag(&int_f, 1.0, a, 0, 1e-7, 8000, 1, int_work, &result, &error);
  
  if(status)
    {
      fprintf(stderr, "integrate_chi_of_a: there is a problem.\n");
    }

  gsl_integration_workspace_free(int_work);

  return -par->clight/par->H0*result;
  
  
}

//===============================
// integrand_chi_s
//===============================
double integrand_chi(double a, void *parameters)
{
  
  struct params *par_int;

  par_int = ((struct params_int *)parameters)->par;
  
  //return 1.0/(pow2(a)*sqrt(par->omega_m/pow3(a) + par->omega_l));
  return 1.0/(pow2(a)*sqrt(par_int->omegam/pow3(a) + par_int->omegal));
      
}
