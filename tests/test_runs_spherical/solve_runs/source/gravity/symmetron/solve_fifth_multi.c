#include<stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../libs/allocator.h"
#include "../../libs/releaser.h"

#include "../../solve.h"

#include "../commons_multigrid/commons_multigrid.h"
#include "../commons_multigrid/restrict.h"
#include "../commons_multigrid/prolonge.h"

#include "solve_fifth_multi.h"
#include "solve_fifth_multi_cycles.h"
#include "source_fifth.h"
#include "solve_fifth_multi_fine.h"

#include "init_phi_fifth.h"  // for teh non-convergent grids.

#include "dump_routines_fifth.h"

#include <omp.h>

// two colors
int sx_first[8] = {0, 1, 1, 0,  1, 0, 0, 1};
int sy_first[8] = {0, 1, 0, 1,  0, 1, 0, 1};
int sz_first[8] = {0, 0, 1, 1,  0, 0, 1, 1};


//============================
// solve_fifth_2d
//============================
int solve_fifth_multigrid(struct params *par, double ***dens, double ***phi, double box, double grid, double a)
{
  // Multigrid solver for the non-linear Bekenstein-Milgrom FIFTH equation 
  // in a 3D uniform grid given a source.
  // See Brent 1977, Trottenberg 2000 or Wesseling 1992.
  //
  // Comments:
  //     - No 4 pi G involved (should come with the source if there is any).
  //     - No units (again, they come with the density).
  //     - User should give and initial condition in phi for the iterations 
  //       (zero is wrong).
  //
      
  //void (*gs_func)(double);  // pointer to the gs function to use.

  struct grids_gravity gg;

  int conv;

  fprintf(stderr, "Entering solve_fifth_multigrid...\n");

  conv = 1;
  
  // Inicialize some constants:
  //===========================
  //gs_func = &gs;  // coloring scheme
  gg.num_grids = log(par->grid)/log(2.0) - 2;  // total number of grids to use.
  // we go until grid 4
  // this is the index of the array...
  gg.g = 0;     // Initial grid (0=finest).
  
  // Allocate grids:
  //================
  allocate_grids_fifth(par, &gg, dens, phi);

  // Do the iterations:
  //===================
  conv = cycle_v_fifth(par, &gg, a);
  
  // Clean up:
  //==========
  deallocate_grids_fifth(par, &gg, dens, phi);

  return 0;
    
}


//======================================
// allocate_grids
//======================================
void allocate_grids_fifth(struct params *par, struct grids_gravity *gg, double ***dens, double ***phi)
{
  // Allocates all the coarser grids.

  int i;
  long size;

  // Set up pointers to discretization of differential operator:
  //============================================================
  switch(par->symm_order)
    {
    case 0:
      // second order
      gg->l_fifth_operator = &l_fifth;
      gg->l_prime_fifth_operator = &l_prime_fifth;
      gg->L_operator = &L;
      break;
    case 1:
      // Fourth order
      gg->l_fifth_operator = &l_fifth_fourth;
      gg->l_prime_fifth_operator = &l_prime_fifth_fourth;
      gg->L_operator = &L_fourth;
      break;
    default:
      fprintf(stderr, "allocate_grids_fifth: wrong symm_order flag.  Aborting.\n");
      exit(0);
      break;
    }

  // Allocate grids:
  //================
  size = par->grid;
  for(i=0; i<=gg->num_grids; i++)
    {
      if(i==0)
	{
	  // Finest grid:
	  // No symmetron
	  //gg->source[i] = dens;  // the thing that comes from outside.
	  // Symmetron
	  gg->rho[i] = dens;  // do not allocate anyhting for source.
	                     // (it will be 0).
	  gg->fi[i] = phi;
	}
      else
	{
	  // Coarse grids:
	  alloc_double_3(&gg->rho[i], size, size, size);
	  alloc_double_3(&gg->source[i], size, size, size);
	  alloc_double_3(&gg->fi[i], size, size, size);
	}
      alloc_double_3(&gg->res[i], size, size, size);
      alloc_double_3(&gg->temp[i], size, size, size);

      alloc_double_3(&gg->temp2[i], size, size, size); // for trunc. err.
      alloc_double_3(&gg->trunc[i], size, size, size);

      gg->n[i] = size;
      gg->h[i] = par->box / ((double)size);
      //fprintf(stderr, "i = %d      size = %ld    h = %lf\n", i, size, gg->h[i]);
      size /= 2;
    }
  
  //fprintf(stderr, "gg->num_grids = %ld\n", gg->num_grids);
  //fprintf(stderr, "\n");
  

 
  
}

//======================================
// deallocate_grids
//======================================
void deallocate_grids_fifth(struct params *par, struct grids_gravity *gg, double ***dens, double ***phi)
{

  // Deallocates all the coarser grids.

  int i;
  long size;

  size = par->grid;
  for(i=0; i<=gg->num_grids; i++)
    {
      if(i!=0)
	{
	  free_double_3(gg->rho[i]);
	  free_double_3(gg->source[i]);
	  free_double_3(gg->fi[i]);
	}
      free_double_3(gg->res[i]);
      free_double_3(gg->temp[i]);

      free_double_3(gg->temp2[i]);
      free_double_3(gg->trunc[i]);
      size /= 2;
    }

}

//======================================
// go_down
//======================================
void go_down_fifth(struct grids_gravity *gg, double a)
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

  //fprintf(stderr, "Going down\n");

  // Restrict solution:
  //===================
  restrict_grid(gg->fi[gg->g-1], gg->temp[gg->g], gg->n[gg->g-1]);
  //write_pot_fifth(499, gg->g-1, a);
  //write_pot_fifth(499, gg->g, a);
  //write_temp(0, gg->g);
    
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
	  gg->temp[gg->g][i][j][k] = gg->fi[gg->g][i][j][k] - gg->temp[gg->g][i][j][k];
#ifdef OPENMP
  }
#endif
  //write_temp_fifth(1, gg->g);

  // Prolongue the fine temp:
  //=========================
  prolonge(gg->temp[gg->g], gg->temp[gg->g-1], gg->n[gg->g], gg->h[gg->g]);
  //write_temp_fifth(2, gg->g-1);

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
	    gg->fi[gg->g-1][i][j][k] += gg->temp[gg->g-1][i][j][k];
	    //if(gg->fi[gg->g-1][i][j][k]<0.0)
	    //  gg->fi[gg->g-1][i][j][k] = -gg->fi[gg->g-1][i][j][k];
	  }
#ifdef OPENMP
  }
#endif
	
 
}

