#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include "../../solve.h"
#include "../../pm/init_parts.h"
#include "../../grid/init_grid.h"

#include "../../libs/allocator.h"

#include "../../pm/calc_delta.h"

#include "../commons_multigrid/commons_multigrid.h"
#include "init_phi_fifth.h"
#include "solve_fifth_multi.h"
//#include "../chamaleon/solve_cham_multi.h"
#include "../commons_multigrid/prolonge.h"


//============================
// test_convergence_fifth
//============================
void test_convergence_fifth(struct params *par, struct grids *grid, struct particles *part, struct units *unit)
{

  double ***delta, ***phi_fifth;
  double h;
  
  h = par->box/par->grid;

  // Allocate space for fine grid:
  //==============================
  alloc_double_3(&delta, 2*par->grid, 2*par->grid, 2*par->grid);
  alloc_double_3(&phi_fifth, 2*par->grid, 2*par->grid, 2*par->grid);

  // Calculate density and interpolate to fine grid:
  //================================================
  calc_delta(par, unit, 
	   part->pos_1, part->pos_2, part->pos_3, 
	   delta, par->mass, part->kernel);

  //prolonge(grid->delta, delta, par->grid, h);

  // Call the solver:
  //=================
  //par->grid = 2*par->grid;
  init_phi_fifth(par, phi_fifth, delta);
  /* if(0) */
  /*   { */
  /*     // Put something different in phi: */
  /*     double rho_crit, rho_mean; */
  /*     rho_crit = 2.776e+7 * pow2(100.0*par->hsmall); */
  /*     rho_mean = par->omegam * rho_crit; */
  /*     for(i=0; i<par->grid; i++)  // x directon */
  /* 	for(j=0; j<par->grid; j++)  // y direction */
  /* 	  for(k=0; k<par->grid; k++)  // z direction */
  /* 	    {  */
  /* 	      //phi[i][j][k] = phi_b(par, par->a);  // symmetron */
  /* 	      phi_fifth[i][j][k] = -pow3(par->symm_assb/par->a)*delta[i][j][k]/rho_mean/2.0; */
  /* 	      phi_fifth[i][j][k] += 1.0;  // for full solution */
  /* 	      if(phi_fifth[i][j][k]<0) */
  /* 		phi_fifth[i][j][k] = 0.0;	     */
  /* 	      phi_fifth[i][j][k] = inverse_change_cham(par->cham_change, phi_fifth[i][j][k]); */
  /* 	    } */
  /*     solve_cham_multigrid(par, delta, phi_fifth, par->box, par->grid, par->a); */
  /*   } */
  /* else */
    {
      solve_fifth_multigrid(par, delta, phi_fifth, par->box, par->grid, par->a);
    }

  fprintf(stderr, "test_convergence_fifth finish.  Exit now.\n");
  exit(0);

}

