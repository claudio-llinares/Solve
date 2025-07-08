/**
 * Here is the V cycle of multigrid iterations.
 * 
 * Since our first grid is number zero, we go to coarse grids by going up,
 * so the cycle is inverted.
 * 
 * Converge criteria:
 *      - See box in page 157 of Trottenberg, Oosterlee, Schuller, "Multigrid" for details.
 *      - The convergence criteria are based on the mean value of the residual
 *      - What we mean by "mean value" is defined in ``grid/get_mean_max_two_fields``.
 *      - We have three criteria:
 *           - Absolute value of the residual (which has units, so it is not good).
 *           - Ratio of residual of latest step and first step.  This is the one implemented in Ramses.  It depends on how good or bad the initial guess is, so it makes no sence.
 *           - Comparison with truncation error.
 **/
#include<stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../libs/allocator.h"
#include "../../libs/releaser.h"

#include "../../solve.h"

#include "../../grid/get_mean_max_two_fields.h"

#include "solve_multigrid_mm.h"
#include "go_up_and_down_mm.h"
#include "gs_residual_truncation_error_mm.h"

void do_one_v_cycle(struct grids_gravity_mm *gg, double a);
void cycle_v_mm_recursive(struct grids_gravity_mm *gg, double a);
void verbose_fine(struct grids_gravity_mm *gg, double a, int write_pot);



int cycle_v_mm(struct grids_gravity_mm *gg, double a)
{
  
  int iter;                   // Iteration number.  If larger than MAX_ITERATION then stop.
  double res_conv;            // Residual to be used for one of the convergence criteria.
  double cur_res, cur_trunc;  // Current residual and truncation error.
  double maximum, max_res, max_trunc;  // Maximum residual and truncation error.
  
  int converged = 0;          // To deal with convergence criteria

  // Show initial conditions:
  //=========================
  //gg->cur_resid[gg->g] = residual_mm(gg, a);   // not sure we need these two lines.
  //gg->cur_trunc[gg->g] = trunc_err_mm(gg, a);
  fprintf(stderr, "NGS = %d   verbose = %d\n", NGS, (gg->par)->verbose_mm);
  if((gg->par)->verbose_mm) verbose_fine(gg, a, 1);

  // Keep first residual for convergence test:
  //==========================================
  residual_mm(gg, gg->g, gg->res_phi[gg->g], gg->res_chi[gg->g]     );
  get_mean_max_two_fields(gg->res_phi[gg->g], gg->res_chi[gg->g], gg->n[gg->g], &res_conv, &maximum);

  // Do iterations:
  //===============
  iter = 1;
  while(1)
    {
      
      // Do one V cycle:
      //================
      do_one_v_cycle(gg, a);

      // Get residual and trunc.err. to test for convergence:
      //=====================================================
      residual_mm(gg, gg->g, gg->res_phi[gg->g], gg->res_chi[gg->g]);
      get_mean_max_two_fields(gg->res_phi[gg->g], gg->res_chi[gg->g], gg->n[gg->g], &cur_res, &max_res);
      trunc_err_mm(gg, gg->g, gg->trunc_phi[gg->g], gg->trunc_chi[gg->g]);
      get_mean_max_two_fields(gg->trunc_phi[gg->g], gg->trunc_chi[gg->g], gg->n[gg->g], &cur_trunc, &max_trunc);
      
      // Show status:
      //=============
      fprintf(stderr, "Step %d    Res = %le   Trunc = %le   Max_trunc = %le   Error = %le  level = %d\n", 
      	      iter, cur_res, cur_trunc, max_trunc, cur_res/res_conv, gg->g);
            
      // Check convergence:
      //===================
      //if( (gg->cur_resid[gg->g] < 1.0e-7 * res_conv) || 
      // 	  (gg->cur_resid[gg->g] < 1.0e-9) )
      //(gg->cur_resid[gg->g] < 0.1*gg->cur_trunc[gg->g]) || 
      if( cur_res < 1.0e-6 * res_conv )     // Relative to initial residual (what Ramses has)
	converged = 1;
      else if(cur_res < 1.0e-9)             // Absolute criterion
	converged = 2;
      //else if( cur_res < 1.0e-3*cur_trunc)  // Relative to truncation error (what Brandt proposed in original paper).
      else if( cur_res < 1.0e-2*cur_trunc)  // Relative to truncation error (what Brandt proposed in original paper).
	converged = 3;
      if(converged)
	{
	  fprintf(stderr, 
		  "cycle_v_mm:  converged.  res = %le   trunc = %le   max_trunc = %le   id = %d\n", 
		  cur_res, cur_trunc, max_trunc, converged);

	  //fprintf(stderr, "Rel. = %g   Abs = %g\n", 
	  //	  gg->cur_resid[gg->g] - 1.0e-7 * res_conv,
	  //	  gg->cur_resid[gg->g] - 1.0e-9);
	  return 1;
	}
      
      if(iter==MAX_ITERATION)
	{
	  fprintf(stderr, "WARNING: solve_multigrid_mm did not converged after %d iterations.\n", iter);
	  if((gg->par)->verbose_mm) verbose_fine(gg, a, 0);
	  
	  return 0;
	}
            
      iter++;
      
    }

}

//=============================
// do_one_v_cycle
//=============================
void do_one_v_cycle(struct grids_gravity_mm *gg, double a)
{

  int i;

  // Coarsest grid:
  //===============
  if(gg->g == gg->num_grids)
    {
      if(gg->verbose) verbose_fine(gg, a, 0);
      for(i=0; i<3*NGS; i++)
	gs_mm(gg);
      if(gg->verbose) verbose_fine(gg, a, 0);
      return;
    }

  // Pre-smoothing:       <================ Equations evolve here!
  //===============
  if(gg->verbose) verbose_fine(gg, a, 1);
  for(i=0; i<NGS; i++)
    gs_mm(gg);
  
  // Calculate residual needed for going up (to check):
  //===================================================
  //gg->cur_resid[gg->g] = residual_mm(gg, a);
  if(gg->verbose) verbose_fine(gg, a, 2);
  
  // Do V cycle from this level:
  //============================
  if( (gg->g==0 && GO_UP_FIRST) || GO_UP_SECOND )
    {
      // Go up:
      //=======
      go_up_mm(gg->par, gg, a);
      gg->g++;
      if(gg->verbose) verbose_fine(gg, a, 0);
      
      // Recursive call to coarse levels:
      //=================================
      do_one_v_cycle(gg, a);
      if(gg->verbose) verbose_fine(gg, a, 0);
      
      // Do down:
      //=========
      go_down_mm(gg, a);
      gg->g--;
      if(gg->verbose) verbose_fine(gg, a, 3);
      
    }
  else
    fprintf(stderr, "do_one_v_cycle:  WARNING: not doing complete V cycle.\n");
  
  // Post-smoothing:       <================ Equations evolve here!
  //================
  for(i=0; i<NGS; i++)
    gs_mm(gg);
  if(gg->verbose) verbose_fine(gg, a, 4);
  

}



//================================
// verbose_fine
//================================
void verbose_fine(struct grids_gravity_mm *gg, double a, int write_pot)
{

  double mean, maximum;

  if (gg->verbose_coarse==0)
    return;
  else
    {
      residual_mm(gg, gg->g, gg->res_phi[gg->g], gg->res_chi[gg->g]     );
      get_mean_max_two_fields(gg->res_phi[gg->g], gg->res_chi[gg->g], gg->n[gg->g], &mean, &maximum);
      //gg->cur_resid[gg->g] = residual_mm(gg, a);
      //gg->cur_trunc[gg->g] = trunc_err_mm(gg, a);
      fprintf(stderr, "verbose_multi = %d  %le  %d\n", gg->g, mean, write_pot);  //, gg->cur_trunc[gg->g], gg->max_trunc[gg->g], gg->g);
      //if (write_pot>0)
      //  write_pot_mm(par, gg, write_pot, gg->g, a);
    }

  return;
}



