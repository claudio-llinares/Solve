#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solve.h"
#include "../grid/init_grid.h"
#include "../libs/utils.h"

#include "mond_mass/solve_mond_mass_multi.h"
#include "mond_mass/operators_mond_mass.h"

#include <omp.h>


void calc_mond_operator(struct params *par, double ***phi, double ***lap_phi, double a);


//======================================
// calc_mond_operator
//======================================
void calc_mond_operator(struct params *par, double ***phi, double ***lap_phi, double a)
{
  // Evaluation of the diference operator using the fi in g.
  // for tunc_err and residual.
  // g = place where you have the potential.
  //
  // The variable used for storing the thing is called lap_phi, but in fact it has the complete non-linear operator.

  double ev_mu[3][2];
  double sum;

  // optimization
  double aa0;
  double fourh;
  
  // Calculates MOND operator

  int ip1, im1, jp1, jm1, kp1, km1;  // i plus 1, etc.
  double h;   // spacing
  int i, j, k;

  int gm1 = par->grid - 1;  // optimization/openmp
  int n = par->grid;

  h = par->box / par->grid;
  fourh = 4.0 * h;
  aa0 = a * par->mm_a0;

  // for each point in the grid.
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(lap_phi, phi, gm1, n, h2)
  {
#pragma omp for schedule(static)
#endif
  for(i=0; i<n; i++)
    {
      ip1 = (i==(gm1) ? 0     : (i+1));
      im1 = (i==(0)   ? (gm1) : (i-1));
      for(j=0; j<n; j++)
	{
	  jp1 = (j==(gm1) ? 0     : (j+1));
	  jm1 = (j==(0)   ? (gm1) : (j-1));
	  for(k=0; k<n; k++)
	    {
	      // periodic stuff
	      kp1 = (k==(gm1) ? 0     : (k+1));
	      km1 = (k==(0)   ? (gm1) : (k-1));
	      
	      // Evaluate the operator:
	      ev_mu[X][MINUS] = mu(dist( ( phi[i][j][k] - phi[im1][j][k] ) / h,
					 ( phi[i][jp1][k] + phi[im1][jp1][k] - phi[i][jm1][k] - phi[im1][jm1][k] ) / (fourh),
					 ( phi[i][j][kp1] + phi[im1][j][kp1] - phi[i][j][km1] - phi[im1][j][km1] ) / (fourh) )/aa0 );
	      
	      ev_mu[X][PLUS]  = mu(dist( ( phi[ip1][j][k] - phi[i][j][k] ) / h,
					 ( phi[ip1][jp1][k] + phi[i][jp1][k] - phi[ip1][jm1][k] - phi[i][jm1][k] ) / (fourh),
					 ( phi[ip1][j][kp1] + phi[i][j][kp1] - phi[ip1][j][km1] - phi[i][j][km1] ) / (fourh) )/aa0 );
	      
	      ev_mu[Y][MINUS] = mu(dist( ( phi[ip1][j][k] + phi[ip1][jm1][k] - phi[im1][j][k] - phi[im1][jm1][k] ) / (fourh),
					 ( phi[i][j][k] - phi[i][jm1][k] ) / h,
					 ( phi[i][j][kp1] + phi[i][jm1][kp1] - phi[i][j][km1] - phi[i][jm1][km1] ) / (fourh) )/aa0 );
	      
	      ev_mu[Y][PLUS]  = mu(dist( ( phi[ip1][jp1][k] + phi[ip1][j][k] - phi[im1][jp1][k] - phi[im1][j][k] ) / (fourh),
					 ( phi[i][jp1][k] - phi[i][j][k] ) / h,
					 ( phi[i][jp1][kp1] + phi[i][j][kp1] - phi[i][jp1][km1] - phi[i][j][km1] ) / (fourh) )/aa0 );
	      
	      ev_mu[Z][MINUS] = mu(dist( ( phi[ip1][j][k] + phi[ip1][j][km1] - phi[im1][j][k] - phi[im1][j][km1] ) / (fourh),
					 ( phi[i][jp1][k] + phi[i][jp1][km1] - phi[i][jm1][k] - phi[i][jm1][km1] ) / (fourh),
					 ( phi[i][j][k] - phi[i][j][km1] ) / h )/aa0 );
	      
	      ev_mu[Z][PLUS]  = mu(dist( ( phi[ip1][j][kp1] + phi[ip1][j][k] - phi[im1][j][kp1] - phi[im1][j][k] ) / (fourh),
					 ( phi[i][jp1][kp1] + phi[i][jp1][k] - phi[i][jm1][kp1] - phi[i][jm1][k] ) / (fourh),
					 ( phi[i][j][kp1] - phi[i][j][k] ) / h )/aa0 );
	      
	      sum = ev_mu[X][MINUS] * (phi[im1][j][k] - phi[i][j][k]) +
		ev_mu[X][PLUS]  * (phi[ip1][j][k] - phi[i][j][k]) +
		ev_mu[Y][MINUS] * (phi[i][jm1][k] - phi[i][j][k]) +
		ev_mu[Y][PLUS]  * (phi[i][jp1][k] - phi[i][j][k]) +
		ev_mu[Z][MINUS] * (phi[i][j][km1] - phi[i][j][k]) +
		ev_mu[Z][PLUS]  * (phi[i][j][kp1] - phi[i][j][k]);

	      lap_phi[i][j][k] = sum/pow2(h);
	    }
	}
    }
#ifdef OPENMP
  }
#endif
  
}
 

