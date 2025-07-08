/***
 * Multigrid solver for the non-linear Bekenstein-Milgrom MM equation 
 * in a 3D uniform grid given a source.
 * See Brent 1977, Trottenberg 2000 or Wesseling 1992.
 *
 * Comments:
 *     - User should give and initial condition in phi for the iterations 
 *       (zero is wrong).
 *     - Grids start from zero, meaning that the V cycle looks more like an 
 *       inverted V (i.e. when we go up, we go to coarse grids). 
**/

#include<stdio.h>
#include <math.h>
#include <stdlib.h>


#include "../../libs/allocator.h"
#include "../../libs/releaser.h"

#include "../../solve.h"

#include "solve_multigrid_mm.h"
#include "cycles_mm.h"
#include "differential_operators_mm.h"

#include <omp.h>

void allocate_grids_mm(struct grids_gravity_mm *gg, double ***dens, double ***chi, double ***phi);
void deallocate_grids_mm(struct grids_gravity_mm *gg, double ***dens, double ***chi, double ***phi);


/**
 * Main routine of the multigrid solver for MOND-MM potentials */
int solve_multigrid_mm(struct params *par, double ***dens, double ***phi, double ***chi, double box, double grid, double a)
{
        
  // This is the structure that has everything:
  //===========================================
  static struct grids_gravity_mm gg;

  int conv;

  fprintf(stderr, "Entering solve_multigrid_mm...\n");

  conv = 1;
  
  // Initialize some constants:
  //===========================
  // Transfer things from par to gg so we do not need to access too many pointers.
  gg.par = par;  // maybe we do not need this
  gg.verbose            = par->verbose_mm;
  gg.verbose_coarse     = par->verbose_mm_coarse;
  gg.cte_poisson_over_a = par->cte_poisson/a;
  gg.aa0                = a*par->mm_a0;
  gg.a2M2               = pow2(a)*pow2(par->mm_M);
  gg.KB                 = par->mm_KB;

  gg.num_grids = log(par->grid)/log(2.0) - 2;  // total number of grids to use.

  // we go until grid 4
  // this is the index of the array...
  gg.g = 0;     // Initial grid (0=finest)   <=============== Important!!!  The cycles are not V, but inverted.

  // Allocate grids:
  //================
  allocate_grids_mm(&gg, dens, phi, chi);

  // Do the iterations:
  //===================
  conv = cycle_v_mm(&gg, a);
  
  // Clean up:
  //==========
  deallocate_grids_mm(&gg, dens, phi, chi);

  return 0;
    
}

/**
 * Allocate additional coarse grids for the multigrid solver.
 * For the fine grid, we do not make a copy, but just keep a pointer 
 * to the variables in the main structure par.
 **/
void allocate_grids_mm(struct grids_gravity_mm *gg, double ***dens, double ***phi, double ***chi)
{
  // Allocates all the coarser grids.

  int i, j, k, l;
  long size;
  
  // Allocate grids:
  //================
  size = (gg->par)->grid;
  for(i=0; i<=gg->num_grids; i++)
    {
      if(i==0)
	{
	  // Finest grid:
	  // No symmetron
	  //gg->source[i] = dens;  // the thing that comes from outside.
	  // Symmetron
	  gg->rho[i] = dens;
	  // The source is not strictly needed for the fine level:
	  // We keep it equal to zero for the shake of simplicity in the Gauss-Seidel part of the code.
	  alloc_double_3(&gg->source_phi[i], size, size, size);
	  alloc_double_3(&gg->source_chi[i], size, size, size);
	  for(j=0; j<size; j++)
	    for(k=0; k<size; k++)
	      for(l=0; l<size; l++)
		{
		  gg->source_phi[i][j][k][l] = 0.0;
		  gg->source_chi[i][j][k][l] = 0.0;
		}
	  gg->phi[i] = phi;  // use variable that comes from the main code (which should have an initial guess).
	  gg->chi[i] = chi;
	}
      else
	{
	  // Coarse grids:
	  alloc_double_3(&gg->rho[i], size, size, size);
	  alloc_double_3(&gg->source_phi[i], size, size, size);
	  alloc_double_3(&gg->source_chi[i], size, size, size);
	  alloc_double_3(&gg->phi[i], size, size, size);
	  alloc_double_3(&gg->chi[i], size, size, size);
	}
      alloc_double_3(&gg->res_phi[i], size, size, size);
      alloc_double_3(&gg->res_chi[i], size, size, size);
      alloc_double_3(&gg->temp_phi[i], size, size, size);
      alloc_double_3(&gg->temp_chi[i], size, size, size);

      alloc_double_3(&gg->temp2_phi[i], size, size, size); // for trunc. err.
      alloc_double_3(&gg->temp2_chi[i], size, size, size); // for trunc. err.
      alloc_double_3(&gg->trunc_phi[i], size, size, size);
      alloc_double_3(&gg->trunc_chi[i], size, size, size);

      gg->n[i]     = size;
      gg->h[i]     = (gg->par)->box / ((double)size);
      gg->fourh[i] = 4.0*gg->h[i];

      size /= 2;
    }
  
}

/** 
 * Clean up the coarse grids.
 * 
 * Comments:
 *  - We do not clean up the fine grid because we only have a pointer to the 
 *    main variables in the main code (in the structure ``grid``).
 *  - We could call this after obtaining a solution.  Which means that we 
 *    have to allocate everything again in the following time step.
 *    I am not sure if that is the best thing to do, but in any case,
 *    the overhead is negligible compare with everything else.
 **/
void deallocate_grids_mm(struct grids_gravity_mm *gg, double ***dens, double ***phi, double ***chi)
{

  // Deallocates all the coarser grids.

  int i;
  long size;

  size = (gg->par)->grid;
  for(i=0; i<=gg->num_grids; i++)
    {
      free_double_3(gg->source_phi[i]);
      free_double_3(gg->source_chi[i]);
      if(i!=0)
	{
	  free_double_3(gg->rho[i]);
	  free_double_3(gg->phi[i]);
	  free_double_3(gg->chi[i]);
	}
      free_double_3(gg->res_phi[i]);
      free_double_3(gg->res_chi[i]);
      free_double_3(gg->temp_phi[i]);
      free_double_3(gg->temp_chi[i]);

      free_double_3(gg->temp2_phi[i]);
      free_double_3(gg->temp2_chi[i]);
      free_double_3(gg->trunc_phi[i]);
      free_double_3(gg->trunc_chi[i]);
      size /= 2;
    }

}

