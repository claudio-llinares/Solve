#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solve.h"
#include "../grid/init_grid.h"

#include "../libs/allocator.h"
#include "../libs/releaser.h"

#include <omp.h>

// Routines to calculate laplacian of potential for many gravities.


//=======================
// calc_lap
//=======================
void calc_lap(struct params *par, double ***phi, double ***lap_phi)
{
  // Calculates laplacian of phi.

  int ip1, im1, jp1, jm1, kp1, km1;  // i plus 1, etc.
  double h, h2;   // spacing
  int i, j, k;

  int gm1 = par->grid - 1;  // optimization/openmp
  int n = par->grid;

  h = par->box / par->grid;
  h2 = pow2(h);

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
	      // do the derivatives
	      lap_phi[i][j][k] = (phi[ip1][j][k] + phi[im1][j][k] +
				  phi[i][jp1][k] + phi[i][jm1][k] +
				  phi[i][j][kp1] + phi[i][j][km1] -
				  6.0*phi[i][j][k]) / h2;
	    }
	}
    }
#ifdef OPENMP
  }
#endif
}


