#include 
#include "../pm/calc_rho.h"
#include "../gravity/commons_multigrid/commons_multigrid.h"
#include "../gravity/fifth/init_phi_fifth.h"
#include "../gravity/fifth/solve_fifth_multi.h"

//===================
// symm_force
//===================
double symm_force(struct params *par, struct grids *grid, struct units *unit)
{

  
  // Store Newtonian fields
  store_fields_newton(par, grid);
  
    // Recover MONDian values from previous grids
  // For the boundaries of the refinements, etc.
  recover_fields_fifth(par, grid);

  // Calculate density
  calc_rho(par, unit, 
	   part->pos_1, part->pos_2, part->pos_3, 
	   grid->delta, par->mass, part->kernel);
  
  // Calculate scalar field
  init_phi_fifth(par, grid);
  solve_fifth_multigrid(par, grid->delta, grid->phi, par->box, par->grid, 1.0);

  // Calculate fifth force
  forces(par, grid->phi, grid->force_1, grid->force_2, grid->force_3, par->a);
  
  // Store symmetron quantities
  store_fields_fifth(par, grid);
  
  // Recover fields Newton (ONLY DURING DEVELOPMENT!)
  recover_fields_newton(par, grid);
  

}
