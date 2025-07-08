/***
 * Calculates force in the grid by deriving a potential.
 ***/
#include<stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solve.h"

#include "get_forces.h"

#define pow2(x) ((x)*(x))

/** 
 * Derives a potential.
 * 
 *  Arguments:
 * 
 *    - par = pointer to ``params`` structure (we need it to know grid size, etc).
 *    - phi = pointer to the potential.
 *    - dphi_dx... = derivatives.
 **/
void forces(struct params *par, double ***phi, double ***dphi_dx, double ***dphi_dy, double ***dphi_dz)
{
  
  // Calculates forces of whatever thing you put in phi

  int i, j, k;
  double h;

  int ip1, im1, jp1, jm1, kp1, km1;

  h = par->box/par->grid;

#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(phi, par, h, dphi_dx, dphi_dy, dphi_dz)
  {
#pragma omp for schedule(static)
#endif
  for(i=0; i<par->grid; i++)
    {
      for(j=0; j<par->grid; j++)
	{
	  for(k=0; k<par->grid; k++)
	    {
	      // periodic stuff:
	      ip1 = (i==(par->grid-1) ? 0            : (i+1));
	      im1 = (i==(0)          ? (par->grid-1) : (i-1));
	      jp1 = (j==(par->grid-1) ? 0            : (j+1));
	      jm1 = (j==(0)          ? (par->grid-1) : (j-1));
	      kp1 = (k==(par->grid-1) ? 0            : (k+1));
	      km1 = (k==(0)          ? (par->grid-1) : (k-1));

	      // do the derivatives:
	      dphi_dx[i][j][k] = (phi[ip1][j][k] - phi[im1][j][k]) / (2.0*h); 
	      dphi_dy[i][j][k] = (phi[i][jp1][k] - phi[i][jm1][k]) / (2.0*h); 
	      dphi_dz[i][j][k] = (phi[i][j][kp1] - phi[i][j][km1]) / (2.0*h); 
	      
	    }
	}    
    }
#ifdef OPENMP
  }
#endif

}
