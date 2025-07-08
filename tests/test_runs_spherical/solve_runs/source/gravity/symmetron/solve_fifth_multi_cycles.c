#include<stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../libs/allocator.h"
#include "../../libs/releaser.h"

#include "../../solve.h"

#include "../commons_multigrid/commons_multigrid.h"
#include "dump_routines_fifth.h"

#include "solve_fifth_multi.h"
#include "solve_fifth_multi_fine.h"
#include "solve_fifth_multi_coarse.h"


void cycle_v_fifth_recursive(struct params *par, struct grids_gravity *gg, double a);

//====================
// cycle_v
//====================
int cycle_v_fifth(struct params *par, struct grids_gravity *gg, double a)
{
  // Cycle V routine.  Do iterations until convergence.  
  // See Brent 1977, Trottenberg 2000 or Wesseling 1992.
  // 
  // Arguments:
  //             a = expansion factor.
  //
  // Comments:
  //             - Assumes user did already a step in the finner grid.
  //

  int i;

  int iter;  // to test convergence
  double res_conv = 0.0;

  int do_multi = 1;

  int verbose = 0;

  // Show initial condition
  gg->cur_resid[gg->g] = residual_fifth_fine(par, gg, a);
  gg->cur_trunc[gg->g] = trunc_err_fine_fifth(gg, a);
  fprintf(stderr, "NGS = %d   verbose = %d\n", NGS, verbose);
  fprintf(stderr, "0: %le  %le  %le  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->max_trunc[gg->g], gg->g);
  if(verbose)
    write_pot_fifth(par, gg, 1, gg->g, a);

  // Check if we are already too close to the solution
  // This is totally stupid...
  if(0)
    {
      if(gg->cur_resid[gg->g] < 1.0e-9)
	{
	  fprintf(stderr, "solve_fifth_cycles:  Residual is too small.  Assume solution is converged before to iterate.\n");
	  return 1;
	}
    }
  
  // Do iterations
  iter = 1;
  while(1)
    {
      // Show residual
      if(verbose) 
	{
	  gg->cur_resid[gg->g] = residual_fifth_fine(par, gg, a);
	  gg->cur_trunc[gg->g] = trunc_err_fine_fifth(gg, a);
	  fprintf(stderr, "1: %le  %le  %le  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->max_trunc[gg->g], gg->g);
	  write_pot_fifth(par, gg, 1, gg->g, a);
	}
      
      // Pre-smoothing:
      for(i=0; i<NGS; i++)
      	gs_fine_fifth(par, gg, a);
      
      // Show residual (we need to calculate the residual to go up anyways)
      gg->cur_resid[gg->g] = residual_fifth_fine(par, gg, a);
      if(verbose) 
	{
	  gg->cur_trunc[gg->g] = trunc_err_fine_fifth(gg, a);
	  fprintf(stderr, "2: %le  %le  %le  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->max_trunc[gg->g], gg->g);
	  write_pot_fifth(par, gg, 2, gg->g, a);
	}
      if(iter==1)  // keep first residual for convergence test
	res_conv = gg->cur_resid[gg->g];
      
      // FAS
      if(do_multi)
	{
	  // Go up
	  go_up_fine_fifth(par, gg, a);
	  gg->g++;
	  
	  // Show residual
	  if(verbose)
	    {
	      gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a); // is this needed?
	      fprintf(stderr, "%le  %le  %ld   C\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->g);
	    }
	  
	  // Recursive call
	  //if(iter<3)
	  cycle_v_fifth_recursive(par, gg, a);
	  //else
	  //fprintf(stderr, "cycle_v_fifth:  WARNING: using only one or two grids.\n");
	  
	  // Show residual
	  if(verbose)
	    {
	      gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a); // is this needed?
	      fprintf(stderr, "%le  %le  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->g);
	    }
	  
	  // Do down
	  go_down_fifth(gg, a);
	  gg->g--;
	  
	}

      // Show residual
      if(verbose) 
	{
	  gg->cur_resid[gg->g] = residual_fifth_fine(par, gg, a);
	  gg->cur_trunc[gg->g] = trunc_err_fine_fifth(gg, a);
	  fprintf(stderr, "3: %le  %le  %le  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->max_trunc[gg->g], gg->g);
	  write_pot_fifth(par, gg, 3, gg->g, a);
	}
      
      // Post-smoothing:
      for(i=0; i<NGS; i++)
	gs_fine_fifth(par, gg, a);
      
      // Get residual and trunc.err. for convergence
      gg->cur_resid[gg->g] = residual_fifth_fine(par, gg, a);
      gg->cur_trunc[gg->g] = trunc_err_fine_fifth(gg, a);
      if(verbose)
	{
	  fprintf(stderr, "4: %le  %le  %le  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->max_trunc[gg->g], gg->g);
	  write_pot_fifth(par, gg, 4, gg->g, a);
	}
      fprintf(stderr, "Step %d    Residual = %le   Trunc = %le   Max_trunc = %le   Error = %le  level = %ld\n", 
	      iter, 
	      gg->cur_resid[gg->g], 
	      gg->cur_trunc[gg->g], 
	      gg->max_trunc[gg->g], 
	      gg->cur_resid[gg->g]/res_conv, 
	      gg->g);
      //fprintf(stderr, "\n");
      
      // Check for convergence
      //if( (gg->cur_resid[gg->g] < 1.0e-7 * res_conv) || 
      // 	  (gg->cur_resid[gg->g] < 1.0e-9) )
      if( (gg->cur_resid[gg->g] < 0.1*gg->cur_trunc[gg->g]) || 
	  (gg->cur_resid[gg->g] < 1.0e-6 * res_conv) || 
	  (gg->cur_resid[gg->g] < 1.0e-9) )
	{
	  fprintf(stderr, "cycle_v_fifth:  converged.  residual_fifth = %le   trunc = %le   max_trunc = %le   ", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->max_trunc[gg->g]);
	  if(gg->cur_resid[gg->g] < 1.0e-9)
	    fprintf(stderr, "(absolute criteria)\n");
	  else
	    fprintf(stderr, "\n");
	  //fprintf(stderr, "Rel. = %g   Abs = %g\n", 
	  //	  gg->cur_resid[gg->g] - 1.0e-7 * res_conv,
	  //	  gg->cur_resid[gg->g] - 1.0e-9);
	  return 1;
	}
      if(iter==MAX_ITERATION_SYMM)
	{
	  fprintf(stderr, "WARNING: solve_fifth_multigrid did not converged after %d iterations.\n", iter);
	  if(verbose)
	    write_pot_fifth(par, gg, 4, gg->g, a);
	  return 0;
	}
      //if(iter==20)
      //	do_multi = 0;
      
      iter++;
      
    }

}


//====================
// cycle__recursive
//====================
void cycle_v_fifth_recursive(struct params *par, struct grids_gravity *gg, double a)
{
  // Cycle V routine.  Do iterations until convergence.  
  // See Brent 1977, Trottenberg 2000 or Wesseling 1992.
  // 
  // Arguments:
  //             a = expansion factor.
  //
  // Comments:
  //             - Assumes user did already a step in the finner grid.
  //

  int i;
  
  int verbose = 0;

  if(verbose)
    fprintf(stderr, "\n");

  // Coarsets grid:
  if(gg->g == gg->num_grids)
    {
      if(verbose) 
	{
	  gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a);
	  fprintf(stderr, "%le  %lf  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->g);
	}
      for(i=0; i<3*NGS; i++)
	gs_coarse_fifth(par, gg, a);
      if(verbose) 
	{
	  gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a);
	  fprintf(stderr, "%le  %lf  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->g);
	}
      return;
    }
    
  // Show residual (we need to calculate residual to go up anyways)
  if(verbose) 
    {
      gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a);
      fprintf(stderr, "%le  %lf  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->g);
      write_pot_fifth(par, gg, 1, gg->g, a);
    }

  // Pre-smoothing:
  for(i=0; i<NGS; i++)
    gs_coarse_fifth(par, gg, a);
  
  // Show residual
  gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a);
  if(verbose) 
    {
      fprintf(stderr, "%le  %lf  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->g);
      write_pot_fifth(par, gg, 2, gg->g, a);
    }
  
  if(1)  // put 0 for two grid algorithm
    {
      // Go up
      go_up_coarse_fifth(par, gg, a);
      gg->g++;
      
      // Show residual (is this one needed?)
      if(verbose)
	{
	  gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a);
	  fprintf(stderr, "%le  %ld\n", gg->cur_resid[gg->g], gg->g);
	}
      
      // Recursive call
      cycle_v_fifth_recursive(par, gg, a);
      
      // Show residual (is this one needed?)
      if(verbose)
	{
	  gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a);
	  fprintf(stderr, "%le  %ld\n", gg->cur_resid[gg->g], gg->g);
	}
      
      // Do down
      go_down_fifth(gg, a);
      gg->g--;
      
    }
  else
    fprintf(stderr, "cycle_v_fifth_recursive:  WARNING: using only two grids.\n");

  // Show residual
  if(verbose) 
    {
      gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a);
      fprintf(stderr, "%le  %lf  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->g);
      write_pot_fifth(par, gg, 3, gg->g, a);
    }

  // Post-smoothing:
  for(i=0; i<NGS; i++)
    gs_coarse_fifth(par, gg, a);
  
  // Show residual
  if(verbose) 
    {
      gg->cur_resid[gg->g] = residual_fifth_coarse(par, gg, a);
      fprintf(stderr, "%le  %lf  %ld\n", gg->cur_resid[gg->g], gg->cur_trunc[gg->g], gg->g);
      write_pot_fifth(par, gg, 4, gg->g, a);
    }

}


