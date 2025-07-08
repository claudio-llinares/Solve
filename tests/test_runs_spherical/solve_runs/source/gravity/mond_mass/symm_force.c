#include 
#include "../pm/calc_rho.h"
#include "../gravity/commons_multigrid/commons_multigrid.h"
#include "../gravity/mond_mass/init_phi_mm.h"
#include "../gravity/mond_mass/solve_mm_multi.h"

//===================
// mm_force
//===================
double mm_force(struct params *par, struct grids *grid, struct units *unit)
{

  
  // Store Newtonian fields
  store_fields_newton(par, grid);
  
    // Recover MONDian values from previous grids
  // For the boundaries of the refinements, etc.
  recover_fields_mm(par, grid);

  // Calculate density
  calc_rho(par, unit, 
	   part->pos_1, part->pos_2, part->pos_3, 
	   grid->delta, par->mass, part->kernel);
  
  // Calculate scalar field
  init_phi_mm(par, grid);
  solve_mm_multigrid(par, grid->delta, grid->phi, par->box, par->grid, 1.0);

  // Calculate mm force
  forces(par, grid->phi, grid->force_1, grid->force_2, grid->force_3, par->a);
  
  // Store mmetron quantities
  store_fields_mm(par, grid);
  
  // Recover fields Newton (ONLY DURING DEVELOPMENT!)
  recover_fields_newton(par, grid);
  

}
