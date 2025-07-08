#include <stdio.h>
#include <math.h>

#include "solve.h"

#include "libs/cosmology.h"

/**
 * Initializes a lot of physical constants including the factor in front of 
 * the density in the Poisson's equation (which somehow contains the effective 
 * gravitational constant).  So if you want to change $G$, this is probably 
 * the place to do it.   Although, keep in mind that additional changes 
 * are made later on to include the dependence with expansion factor.
 *
 * Parameters:
 *     par: parameter structure.
 *     unit: units structure.
 *
 * Returns:
 *     Changes are included in the two structures that you input.
 */
void init_physics(struct params *par, struct units *unit)
{

  unit->pi      = 4.0*atan(1.0);
  unit->twopi   = 2.0*unit->pi;
  unit->cm2mpc  = 1.0/3.08568025e24;
  unit->m2mpc   = 1.0/3.08568025e22;
  unit->mpc2m   = 3.08568025e22;
  unit->kg2msun = 1.0/(1.9198e+30);
  unit->sec2gyr = 1.0 / (1.0e+9 * 365.25 * 23.56 * 3600.0);
  unit->gyr2sec = 1.0e+9 * 365.25 * 23.56 * 3600.0;
  unit->km2mpc  = 1.0/3.08568025e19;
  unit->mpc2km  = 3.08568025e19;
    
  unit->mpc2km_over_gyr2sec = unit->mpc2km / unit->gyr2sec;

  unit->H0 = 100.0 * par->hsmall * unit->km2mpc / unit->sec2gyr;    // gyr^-1
  par->H0  = unit->H0;

  // Get rid of hsmall in all the quantities:
  par->box  /= par->hsmall;
  par->mm_M *= par->hsmall;

  par->nparts = pow3(par->grid_p);

  par->a = 1.0/(1.0+par->z_init);

  unit->cte_poisson = 1.5 * pow2(par->H0) * par->omegam;
  par->cte_poisson  = unit->cte_poisson;

  if (par->gravity==1)
    {
      fprintf(stderr, "Warning!  init_physics: changing cte_poisson!\n");
      par->cte_poisson = par->cte_poisson/15.0;
    }

  par->mm_a0 = par->mm_a0*unit->cm2mpc/pow2(unit->sec2gyr);
  
  // Tabulate t_of_a:
  //=================
  // We need this so we can inverted later on to get a(t).
  tabulate_t_of_a(par);

  fprintf(stderr, "H0 = %lg gyr^-1\n", unit->H0);
  fprintf(stderr, "a0 = %lg Mpc / Gyr^2\n", par->mm_a0);
  fprintf(stderr, "\n");

}
