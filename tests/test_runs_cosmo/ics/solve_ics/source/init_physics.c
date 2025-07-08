#include <stdio.h>
#include <math.h>

#include "solve.h"

//====================
// init_physics
//====================
void init_physics(struct params *par, struct units *unit)
{

  unit->pi      = 4.0*atan(1.0);
  unit->twopi   = 2.0*unit->pi;
  //unit->G_mks   = 6.67300e-11;     // m+3 kg-1 s-2
  unit->G_mks   = 6.91125e-11;   // grafic
  unit->m2mpc   = 1.0/(3.08568025e22);
  unit->mpc2m   = 3.08568025e22;
  unit->kg2msun = 1.0/(1.9198e+30);
  unit->sec2gyr = 1.0 / (1.0e+9 * 365.25 * 23.56 * 3600.0);
  unit->gyr2sec = 1.0e+9 * 365.25 * 23.56 * 3600.0;
  unit->km2mpc  = 1.0/3.08568025e19;
  unit->mpc2km  = 3.08568025e19;
    
  unit->H0 = 100.0 * par->hsmall * unit->km2mpc / unit->sec2gyr;    // gyr^-1
  par->H0 = unit->H0;
  unit->fourpiG = 4.0 * unit->pi * unit->G_mks * pow3(unit->m2mpc) / unit->kg2msun / pow2(unit->sec2gyr);
  unit->clight = 2.99792458e+8 * unit->m2mpc / unit->sec2gyr;
  par->clight = unit->clight;
  unit->cte_poisson = 3.0/2.0 * par->omegam * pow2(unit->H0);  // gyr^-2 
  fprintf(stderr, "unit->cte_poisson = %lf\n", unit->cte_poisson);

  // Get rid of hsmall in all teh quantities:
  par->box /= par->hsmall;
  par->nparts = pow3(par->grid_p);
  
  // My version.
  //unit->rho_crit = (3.0*pow2(unit->H0)) / 
  //  (8.0*unit->pi*unit->G_mks*pow3(unit->m2mpc)/unit->kg2msun/pow2(unit->sec2gyr));
  //unit->rho_mean = par->omegam * unit->rho_crit;
  // grafic1
  unit->rho_crit = 2.776e+7 * pow2(100.0*par->hsmall);
  unit->rho_mean = par->omegam * unit->rho_crit;

  fprintf(stderr, "cte_poisson = %g \n", unit->cte_poisson);
  fprintf(stderr, "H0 = %lg gyr^-1\n", unit->H0);
  fprintf(stderr, "rho_crit = %lg msun mpc^-3 = %lg msun Mpc^-3 h^2\n", 
	  unit->rho_crit, unit->rho_crit/pow2(par->hsmall));
  fprintf(stderr, "rho_mean = %lg msun mpc^-3 = %lg msun Mpc^-3 h^2\n", 
	  unit->rho_mean, unit->rho_mean/pow2(par->hsmall));
  fprintf(stderr, "G = %lg mpc^3 msun^-1 sec^-2\n", unit->G_mks * pow3(unit->m2mpc)/unit->kg2msun/pow2(unit->sec2gyr));
  fprintf(stderr, "clight = %g mpc/gyr\n", unit->clight);


  // Determine times:
  //=================
  /* a = 1.0/(1.0+par->z_init); */
  /* a_end = 1.0/(1.0+par->z_end); */
  /* da = (a_end-a)/par->nsteps; */
  /* delta_a_out = (a_end-a) / par->n_outputs; */
  /* a_out = a; */
	 
  /* l_out = 0; */
  
  /* fprintf(stderr, "Time parameters:\n"); */
  /* fprintf(stderr, "  a           = %le\n", a); */
  /* fprintf(stderr, "  da          = %le\n", da); */
  /* fprintf(stderr, "  a_end       = %le\n", a_end); */
  /* fprintf(stderr, "  a_out       = %le\n", a_out); */
  /* fprintf(stderr, "  delta_a_out = %le\n", delta_a_out); */
  /* fprintf(stderr, "  l_out       = %d\n", l_out); */
  /* fprintf(stderr, "\n"); */

}
