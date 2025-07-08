#include "../../solve.h"
#include "../../grid/init_grid.h"

#include "solve_mond_mass_multi.h"
#include "../commons_multigrid/restrict.h"
#include "../commons_multigrid/prolonge.h"   // for truncation error

#include "solve_mond_mass_multi_fine.h"
#include "operators_mond_mass.h"

//======================================
// go_up_fine_mm
//======================================
void go_up_mm(struct params *par, struct grids_gravity_mm *gg, double a)
{
  //
  // Moves everything to the next coarser grid.
  //      - Calculates source in coarse grid.
  //      - Transfer density and solution.
  //
  // The source is calculated with:
  //                                   S = residual(R(phi), R(rho), R(S)) - R(residual(phi, rho, S))
  // 
  // For first term:  
  //                    - Restrict three fields (including the source in this level, which is zero in the fine grid).
  //                    - Calculate residuals with these fields.
  //                    - The restrictions has to the done outside the loop below, but the calculation of the residual can be calculated on the fly.
  // For the second term:
  //                    - Calculate the residual in the level where we are.
  //                    - Restrct that residual.
  //  
    
  int i, j, k;

  // optimization:
  int nm1_p1 = gg->n[gg->g+1]-1;
  //  double val_phi, val_chi;

  int grid, grid_p1;

  if((gg->par)->verbose_mm) fprintf(stderr, "Going up %d -> %d\n", gg->n[gg->g], gg->n[gg->g+1]);
 

  // Copy levels to local variables:
  //================================
  grid    = gg->g;
  grid_p1 = gg->g+1;

  // Restrict fields for first term:
  //================================
  restrict_grid(gg->rho[grid],        gg->rho[grid_p1],        gg->n[grid]);
  restrict_grid(gg->source_phi[grid], gg->source_phi[grid_p1], gg->n[grid]);
  restrict_grid(gg->source_chi[grid], gg->source_chi[grid_p1], gg->n[grid]);
  restrict_grid(gg->phi[grid],        gg->phi[grid_p1],        gg->n[grid]);
  restrict_grid(gg->chi[grid],        gg->chi[grid_p1],        gg->n[grid]);

  // Calculate residual in this level (and restrict it):
  //====================================================
  residual_mm(gg, grid, gg->res_phi[grid], gg->res_chi[grid]     );
  restrict_grid(gg->res_phi[grid], gg->res_phi[grid_p1], gg->n[grid]);
  restrict_grid(gg->res_chi[grid], gg->res_chi[grid_p1], gg->n[grid]);

  // Calculate residual in coarse level:
  //====================================
  residual_mm(gg, grid_p1, gg->temp_phi[grid_p1], gg->temp_chi[grid_p1] );

  // Substract both residuals:
  //==========================
#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(gg, grid, grid_p1)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<= nm1_p1; i++)   // Fills all the values in the grid.
      for(j=0; j<= nm1_p1; j++)
	for(k=0; k<= nm1_p1; k++)
	  {
	    gg->source_phi[grid_p1][i][j][k] = gg->temp_phi[grid_p1][i][j][k] - gg->res_phi[grid_p1][i][j][k];
	    gg->source_chi[grid_p1][i][j][k] = gg->temp_chi[grid_p1][i][j][k] - gg->res_chi[grid_p1][i][j][k];

	    // Delete the density (we do not need it):
	    // but we did use it to calculate the residuals for moving things up.
	    gg->rho[grid_p1][i][j][k] = 0.0;
	  }
#ifdef OPENMP
  }
#endif

}


//======================================
// go_down
//======================================
void go_down_mm(struct grids_gravity_mm *gg, double a)
{

  //
  // Moves everything to the next finer grid
  //
  // Comments:
  //            - gg->g is the coarse.
  //            - gg->g-1 fine
  
  int i, j, k;

  // optimization
  int n     = gg->n[gg->g];
  int n_m1  = gg->n[gg->g-1];

  if((gg->par)->verbose_mm) fprintf(stderr, "Going down %d -> %d\n", n, n_m1);

  // Restrict solution:
  //===================
  restrict_grid(gg->phi[gg->g-1], gg->temp_phi[gg->g], gg->n[gg->g-1]);
  restrict_grid(gg->chi[gg->g-1], gg->temp_chi[gg->g], gg->n[gg->g-1]);
    
  // Substracts temp to coarse solution:
  //====================================
  // Other good boundaries:  1, -2
#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(gg, n)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i< n; i++)
      for(j=0; j< n; j++)
	for(k=0; k< n; k++)
	  {
	    gg->temp_phi[gg->g][i][j][k] = gg->phi[gg->g][i][j][k] - gg->temp_phi[gg->g][i][j][k];
	    gg->temp_chi[gg->g][i][j][k] = gg->chi[gg->g][i][j][k] - gg->temp_chi[gg->g][i][j][k];
	  }
#ifdef OPENMP
  }
#endif

  // Prolongue the fine temp:
  //=========================
  prolonge(gg->temp_phi[gg->g], gg->temp_phi[gg->g-1], gg->n[gg->g], gg->h[gg->g]);
  prolonge(gg->temp_chi[gg->g], gg->temp_chi[gg->g-1], gg->n[gg->g], gg->h[gg->g]);

  // Add this stuff to the fine solution:
  //=====================================
#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(gg, n_m1)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i< n_m1; i++)
      for(j=0; j< n_m1; j++)
	for(k=0; k< n_m1; k++)
	  {
	    gg->phi[gg->g-1][i][j][k] += gg->temp_phi[gg->g-1][i][j][k];
	    gg->chi[gg->g-1][i][j][k] += gg->temp_chi[gg->g-1][i][j][k];
	  }
#ifdef OPENMP
  }
#endif
 
}

