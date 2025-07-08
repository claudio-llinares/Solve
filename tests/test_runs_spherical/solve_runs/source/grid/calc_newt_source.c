#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <fftw3.h>

#include "../solve.h"
#include "init_grid.h"

//===============================
// analytic_rho_hernquist
//===============================
void calc_newt_source(struct params *par, struct grids *grid, struct units *unit)
{
  // Calculates Newtonian source given an overdensity.
  //
  // Arguments:
  //              rd = flag to put rho or delta.
  // Comments:
  //           - delta should be allocated outside.

  int i, j, k;
  double move;
  
  fprintf(stderr, "Entering calc_newt_source...\n");

  move = unit->cte_poisson/par->a;  // optimization
#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(par, grid, move)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<par->grid; i++)  // x directon
      for(j=0; j<par->grid; j++)  // y direction
	for(k=0; k<par->grid; k++)  // z direction
	  grid->source[i][j][k] = move * (grid->delta[i][j][k]);
#ifdef OPENMP
  }
#endif
  
}
