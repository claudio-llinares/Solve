#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <fftw3.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"
#include "../../libs/D_lib.h"

//===============================
// analytic_rho_hernquist
//===============================
void calc_newt_source(struct params *par, struct grids *grid, struct units *unit, double a)
{
  // Calculates Newtonian source given an overdensity.
  //
  // Arguments:
  //              rd = flag to put rho or delta.
  // Comments:
  //           - delta should be allocated outside.

  int i, j, k;
  double move;
  //double D_initial;
  
  fprintf(stderr, "Entering calc_newt_source...\n");

  //D_initial = D_a(a, a, par->omegam, par->omegar, par->omegal);
  //fprintf(stderr, "calc_newt_source: D for source = %lf\n", D_initial);
  //fprintf(stderr, "calc_newt_source: normalizing growth functions to 1 at z=%lf\n", 1.0/(a)-1.0);
  move = unit->cte_poisson/a;  // optimization
#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(par, grid, move)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<par->grid; i++)  // x directon
      for(j=0; j<par->grid; j++)  // y direction
	for(k=0; k<par->grid; k++)  // z direction
	  //grid->source[i][j][k] = move * D_initial * grid->delta[i][j][k];
	  grid->source[i][j][k] = move * grid->delta[i][j][k];
#ifdef OPENMP
  }
#endif
  
}
