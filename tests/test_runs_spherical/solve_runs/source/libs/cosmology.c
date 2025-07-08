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
;//=======================
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

  return -pow2(unit->H0)/2.0 * ( 3.0*par->omegam/pow3(a) );

}

//=======================
// a_dot
//=======================
double a_dot(struct params *par, struct units *unit, double a)
{

  return unit->H0 * sqrt(par->omegam/a + par->omegal * pow2(a));
  
}

//=========================
double F(struct params *par, struct units *unit, double a)
{

  // A dimensionless version of a_dot

  //return sqrt(par->omegam/a + pow2(a)*par->omegal);
  //return pow(a, 1.0+par->p) * sqrt(par->omegam/pow3(a) + par->omegal);
  return sqrt(par->omegam/pow3(a) + par->omegal);


}

//=========================
double A(struct params *par, struct units *unit, double a)
{

  // A dimensionless version of (a_dot/a)^2

  return (par->omegam/pow3(a) + par->omegal);


}

//=========================
// integrate_time_of_a
//=========================
double integrate_time_of_a(struct params *par, double a)
{

  gsl_integration_workspace *int_work;
  gsl_function int_f;
  double result, error;

  struct params_int par_int;

  int status;

  int_work = gsl_integration_workspace_alloc(8000);
  int_f.function = integrand_time;
  int_f.params = &par_int;
  par_int.par = par;

  status = gsl_integration_qag(&int_f, 1.0e-6, a, 0, 1e-7, 8000, 1, int_work, &result, &error);
  
  if(status)
    {
      fprintf(stderr, "integrate_chi_of_a: there is a problem.\n");
    }

  gsl_integration_workspace_free(int_work);

  return result/par->H0;
  
  
}



//===============================
// integrand_chi_s
//===============================
double integrand_time(double a, void *parameters)
{
  
  // This is conformal time!!!

  struct params *par_int;

  par_int = ((struct params_int *)parameters)->par;
  
  return 1.0/(pow2(a)*sqrt(par_int->omegam/pow3(a) + par_int->omegal));  // Conformal
  //return 1.0/(a*sqrt(par_int->omegam/pow3(a) + par_int->omegal));  // Newtonian
      
}

//===================================================================
//===================================================================
//===================================================================


//===============================
// tabulate_chi_of_a
//===============================
void tabulate_t_of_a(struct params *par)
{

  int i;
  double a;

  double a_min = 1.0e-6;
  double n_t = 1000.0;
  
  //for(i=0; i<n_t-1; i++)
  for(i=0; i<n_t; i++)
  {
    a = a_min + (1.1-a_min)/(n_t-1) * i;  //(model->n_chi-i-1);
    
    par->t_of_a[i] = integrate_time_of_a(par, a);
    
    //fprintf(stderr, "%lg  %lg\n", a, par->t_of_a[i]);
    
  }
  
}

//===============================
// a_of_t
//===============================
double a_of_t(struct params *par, double t)
{

  int i, ii, i_up, i_down;
  double a, a_down, a_up;
  double a_min = 1.0e-6;
  double n_t = 1000.0;
  
  // Find interval
  i_up = n_t-1;  // a=1
  i_down = 0;             // a=a_s
  ii = (int)(n_t-0)/2;
  for(i=0; i<n_t-1; i++)
    {
      a_down = a_min + (1.1-a_min)/(n_t-1) * i_down;
      a_up   = a_min + (1.1-a_min)/(n_t-1) * i_up;
      //fprintf(stderr, "%4d  %4d  %4d | %7.5lg  %7.5lg  %7.5lg  %7.5lg | %lg  %lg\n", 
      //    i_down, ii, i_up, 
      //    par->t_of_a[i_down], t, par->t_of_a[ii], par->t_of_a[i_up], 
      //    a_down, a_up);
      
      if(t==par->t_of_a[i_down])
	{
	  i_up = i_down+1;
	  //fprintf(stderr, "break 1  %d  %d\n", i_down, i_up);
	  break;
	}
      if(t==par->t_of_a[i_up])
	{
	  i_down = i_up-1;
	  //fprintf(stderr, "break 2  %d  %d\n", i_down, i_up);
	  break;
	}
      
      if(t>=par->t_of_a[ii] && t<=par->t_of_a[ii+1])
	{
	  if(t==par->t_of_a[ii])
	    {
	      i_down = ii;
	      i_up = ii+1;
	    }
	  if(t==par->t_of_a[ii+1])
	    {
	      i_down = ii;
	      i_up = ii+1;
	    }
	  
	  i_down = ii;
	  i_up = ii+1;
	  //fprintf(stderr, "break 3  %d  %d  %lg  %lg\n", i_down, i_up, par->t_of_a[ii], par->t_of_a[ii+1]);
	  break;
	}
      
      if(t<par->t_of_a[ii])
	i_up = ii;
      else
	i_down = ii;
      
      ii = (int)((i_up+i_down)/2.0+0.5);
	
    }

  // Interpolate
  a_down = a_min + (1.1-a_min)/(n_t-1) * i_down;
  a_up   = a_min + (1.1-a_min)/(n_t-1) * i_up;
  //fprintf(stderr, "DONE %lg  %lg  %lg\n", t, a_down, a_up);
  a = (t - par->t_of_a[i_down])/(par->t_of_a[i_up] - par->t_of_a[i_down]) * (a_up - a_down) + a_down;

  //fprintf(stderr, "%lg  %lg  %lg\n", a_down, a, a_up);
    
  return a;
  
}
