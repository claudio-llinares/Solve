#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "../solve.h"
#include "init_parts.h"
#include "../grid/init_grid.h"

#include "../gravity/newton/calc_newt_source.h"
#include "../gravity/newton/solve_poisson_fftw.h"
#include "../gravity/get_forces.h"

#include "../libs/D_lib.h"
#include "../libs/utils.h"
#include "../libs/cosmology.h"

#include "../pm/interp_force.h"

void zeldovich(struct params *par, struct units *unit, struct grids *grid, struct particles *part, double a, double (*kernel)(double, double, double, double, double, double, struct params *))
{


  // Move particles and velocities according to Zeldovich approximation
  // from z_init to z_zeldovich.
  // Actualizes also the value of a to be z_zeldovich.
  //

  double move, move_vel;

  int i, j, k;
  int i_part;
  double force_i1, force_i2, force_i3;
  
  double D, D_vel, norm_D;

  double da_over_two;
  
  // Calculate potential and force:
  //===============================
  calc_newt_source(par, grid, unit, 1.0);
  solve_poisson_fftw(par, grid, unit);
  forces(par, grid->phi, grid->force_1, grid->force_2, grid->force_3);
      
  // Do Zeldovich desplacements:
  //============================
  da_over_two = (1.0-par->a_end)/par->nsteps_nbody/2.0;
  D = D_a(par->a_end, par->omegam, par->omegar, par->omegal, par->H0, &norm_D);
  D_vel = D_a(par->a_end+da_over_two, par->omegam, par->omegar, par->omegal, par->H0, &norm_D);
  fprintf(stderr, "zeldovich: D      = %lg\n", D);
  fprintf(stderr, "zeldovich: norm_D = %lg\n", norm_D);
  move = D/unit->cte_poisson;
  move_vel = 1.0/(pow2(par->a_end+da_over_two)*hubble(par, unit, par->a_end+da_over_two)) * 
    (D_vel/(par->a_end+da_over_two) - 1.0/(unit->cte_poisson*(norm_D)));
  fprintf(stderr, "move     = %le\n", move);
  fprintf(stderr, "move_vel = %le\n", move_vel);

  i_part = 0;
  for(i=0; i<par->grid_p; i++)
    for(j=0; j<par->grid_p; j++)
      for(k=0; k<par->grid_p; k++)  // z direction
	{
	  // Grid:
	  part->pos_1[i_part] = i*par->box/par->grid_p + par->box/par->grid_p/2.0;
	  part->pos_2[i_part] = j*par->box/par->grid_p + par->box/par->grid_p/2.0;
	  part->pos_3[i_part] = k*par->box/par->grid_p + par->box/par->grid_p/2.0;

	  // Move particles:
	  interp_force(par,
		       grid->force_1, grid->force_2, grid->force_3,
		       part->pos_1[i_part], part->pos_2[i_part], part->pos_3[i_part],
		       &force_i1, &force_i2, &force_i3, kernel);
	  
	  part->pos_1[i_part] = modulus( part->pos_1[i_part] - move*force_i1, par->box );
	  part->pos_2[i_part] = modulus( part->pos_2[i_part] - move*force_i2, par->box );
	  part->pos_3[i_part] = modulus( part->pos_3[i_part] - move*force_i3, par->box );
	  	  
	  part->vel_1[i_part] = move_vel * force_i1;
	  part->vel_2[i_part] = move_vel * force_i2;
	  part->vel_3[i_part] = move_vel * force_i3;
	  
	  // Construct momentum:
	  // DON'T FORGET, THE vel VARIABLE IS A MOMENTUM!!!!
	  part->vel_1[i_part] *= pow2(par->a_end);
	  part->vel_2[i_part] *= pow2(par->a_end);
	  part->vel_3[i_part] *= pow2(par->a_end);

	  i_part++;
	  
	}


}
