/***
 * Calculates forces on the grid from scratch (including calculation of the 
 * density and potential).
 * 
 * Comments:
 *     - The result is returned in the corresponding variable in the grid 
 *       structure.
 *     - The density and potential will also end up in the coresponding 
 *       variables in the structure grid.
 ***/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../solve.h"
#include "../pm/init_parts.h"
#include "../grid/init_grid.h"

#include "../grid/add_two_fields.h"

#include "../libs/cosmology.h"

#include "../pm/calc_delta_cic.h"

#include "../gravity/newton/calc_newt_source.h"
#include "../gravity/newton/fft_solver/solve_poisson_fftw.h"
#include "../gravity/get_forces.h"

#include "../gravity/mond_mass/initial_guess_mm.h"
#include "../gravity/mond_mass/solve_multigrid_mm.h"



//===========================
// get_gravity
//===========================
void get_gravity(struct params *par, struct grids *grid, struct particles *part, struct units *unit, double t)
{

  static int first = 1;

  // Calculate overdensity:
  //=======================
  calc_delta_cic(par, unit, part->pos_1, part->pos_2, part->pos_3, grid->delta);

  switch(par->gravity)
    {
    case 0:
      if(first)
	{
	  first = 0;
	  fprintf(stderr, "Using Newtonian gravity...\n");
	}
      // Newtonian gravity:
      //===================
      // Add some constants to the density to get the source of Poisson's equation
      calc_newt_source(par, grid, unit, a_of_t(par, t));  // puts source in grid->source
      // Solve Poisson's equation:
      solve_poisson_fftw(par, grid, unit);
      break;
    case 1:
      if(first)
	{
	  first = 0;
	  fprintf(stderr, "Using MOND with mass gravity...\n");
	}
      // MOND with mass:
      //================
      // Get initial guess for multigrid solver:
      initial_guess_mm(par, grid->chi_bar_mm, grid->phi_bar_mm);
      // Calculate the two potentials:
      solve_multigrid_mm(par, grid->delta, grid->phi_bar_mm, grid->chi_bar_mm, par->box, par->grid, par->a);
      // Add the two potentials and put it where the metric pertubation (i.e. the newtonian potential) goes:
      add_two_fields(grid->phi_bar_mm, grid->chi_bar_mm, par->grid, grid->phi);
      break;
    default:
      fprintf(stderr, "Wrong gravity flag (%d).  Aborting.\n", par->gravity);
      exit(0);
    }
      
  // Calculate the force from the potential:
  //========================================
  forces(par, grid->phi, grid->force_1, grid->force_2, grid->force_3);

  return;

}
