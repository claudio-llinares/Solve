#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"

#include "../../libs/allocator.h"
#include "../../libs/releaser.h"

#include "../commons_multigrid/commons_multigrid.h"
#include "../fifth/init_phi_fifth.h"
#include "../fifth/solve_fifth_multi.h"

#include "../../grid/write_field.h"

void write_symm(struct params *par, struct grids *grid, int l_out)
{

  // Calculates and writes static symmetron field
  // given an overdensity.

  double ***rho, ***static_phi;

  int i, j, k;

  double rho_crit, rho_mean;
  rho_crit = 2.776e+7 * pow2(100.0*par->hsmall);
  rho_mean = par->omegam * rho_crit;

  fprintf(stderr, "Entering write_symm...\n");

  // Allocates density and solution:
  //================================
  alloc_double_3(&rho,      par->grid, par->grid, par->grid);
  alloc_double_3(&static_phi, par->grid, par->grid, par->grid);

  // Fills density:
  //===============
  for(i=0; i<par->grid; i++)
	  for(j=0; j<par->grid; j++)
	    for(k=0; k<par->grid; k++)  // z direction
	      rho[i][j][k] = (grid->delta[i][j][k] + 1.0)*rho_mean;

  // Solve equation:
  //================
  init_phi_fifth(par, static_phi, rho);
  solve_fifth_multigrid(par, rho, static_phi, par->box, par->grid, par->a);

  // Output:
  //========
  write_only_field(par,
		   static_phi,
		   "static_fifth_phi", 0.0, l_out);

  free_double_3(rho,         par->grid, par->grid, par->grid);
  free_double_3(static_phi, par->grid, par->grid, par->grid);

}
