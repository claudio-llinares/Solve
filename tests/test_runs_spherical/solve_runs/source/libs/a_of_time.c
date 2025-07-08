#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
//#include <gsl/gsl_odeiv.h>

#include "../solve.h"

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))

//=================================
// a_dot_supercom_for_rk
//=================================
int a_dot_supercom_for_rk(double t, const double y[], double f[], void *params)
{
  
  double *mu = (double *)params;
  
  f[0] = mu[0] * pow(y[0], 1.5) * sqrt(mu[1] + mu[2]*pow3(y[0]));

  //fprintf(stderr, "%g  %g  %g  %g  %g\n", 
  //	  mu[0], mu[1], mu[2], y[0], f[0]);

  return GSL_SUCCESS;
 
}


//=================================
// jac_a_dot_supercom_for_rk
//=================================
int jac_a_dot_supercom_for_rk(double t, const double y[], double *dfdy, double dfdt[], void *params)
{

  double *mu = (double *)params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 1, 1);
  gsl_matrix * m = &dfdy_mat.matrix;

  gsl_matrix_set (m, 0, 0,
  		  1.5 * sqrt(y[0]) * mu[0] *
  		  ( (mu[1]+2.0*pow3(y[0])*mu[2]) /
  		    sqrt( mu[1] + pow3(y[0])*mu[2])  )   );
  
  dfdt[0] = 0.0;
  
  return GSL_SUCCESS;
  
}

//=============================
// a_of_supercom_time
//=============================
double a_of_supercom_time(struct params *par, struct units *unit, double tau)
{
  
  //
  // Calculates expansion factor that corresponds to a given
  // supercomoving time t.
  //
  // Assumes 
  //          a(today) = 0
  //          tau(today) = 0
  // Times are negative for the past.
  //

  double mu[3];

  gsl_odeiv2_system sys = {a_dot_supercom_for_rk, 
			   jac_a_dot_supercom_for_rk, 
			   1, 
			   &mu};
  gsl_odeiv2_driver *d; 

  double t;
  int status;

  double y[1] = {1.0};  //1.0e-6};

  d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, -1e-6, 1e-6, 0.0);
  t = 0.0;   //1.0e-6;

  mu[0] = unit->H0;
  mu[1] = par->omegam;
  mu[2] = par->omegal;

  status = gsl_odeiv2_driver_apply(d, &t, tau, y);
  if (status != GSL_SUCCESS)
    {
      fprintf (stderr, "a_of_supercom_time: error, return value=%d\n", status);
      exit(0);
    }
  fprintf(stderr, "a_of_supercom_time: tau = %g  %g  a(tau) = %g\n", t, tau, y[0]);
  
  gsl_odeiv2_driver_free(d);

  return y[0];

}

//===========================
// time_of_a_supercom
//===========================
double time_of_a_supercom(struct params *par, struct units *unit, double a_init)
{

  //
  // Calculates supercomoving time that corresponds to a 
  // given expansion factor.
  //
  // Assumes 
  //          a(today) = 0
  //          tau(today) = 0
  // Times are negative for the past.
  //

  double mu[3];

  gsl_odeiv2_system sys = {a_dot_supercom_for_rk, 
			   jac_a_dot_supercom_for_rk, 
			   1, 
			   &mu};
  gsl_odeiv2_driver *d; 

  double t;
  int status;

  double tau, delta;
  //double tau_old;
  
  double y[1] = {1.0};  //1.0e-6};


  fprintf(stderr, "Calculating initial time.  a_init = %g\n", a_init);

  d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, -1e-6, 1e-6, 0.0);

  mu[0] = unit->H0;
  mu[1] = par->omegam;
  mu[2] = par->omegal;

  // Find a time to the left
  delta = 1.0e-4;
  t = 0.0;
  tau = -2.0*delta;
  //  tau_old = -2.0*delta;
  while(1)
    {
      // Calculate a corresponding to tau.
      
      status = gsl_odeiv2_driver_apply(d, &t, tau, y);
      if (status != GSL_SUCCESS)
	{
	  fprintf (stderr, "time_of_a_supercom:  error, return value=%d\n", status);
	  exit(0);
	}
      //fprintf(stderr, "tau = %g   a_init = %g    a(tau) = %g\n", 
      //	      tau, a_init, y[0]);
      if(y[0] < a_init)
	{
	  fprintf(stderr, "time_of_a_supercom: tau = %g   a_init = %g    a(tau) = %g\n", 
		  tau, a_init, y[0]);
	  break;
	}
      else
	{
	  //y[0] = 1.0;
	  //t = 0.0;
	  //tau_old = tau;
	  tau -= delta;
	}
    }

  //fprintf(stderr, "tau = %g  %g  a(tau) = %g\n", t, tau, y[0]);
  
  gsl_odeiv2_driver_free(d);

  return tau;



}



//===========================
// test_a_of_time
//===========================
void test_a_of_time(struct params *par, struct units *unit)
{
  
  double tau;

  for(tau=0.01; tau<13.0; tau++)
    {
      //fprintf(stderr, "test: %g  %g\n", 
      //      tau,  
      //      a_of_supercom_time(par, unit, tau));

      a_of_supercom_time(par, unit, tau);
    }


}
