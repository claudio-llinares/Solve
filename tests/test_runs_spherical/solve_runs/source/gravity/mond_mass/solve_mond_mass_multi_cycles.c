#include<stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../libs/allocator.h"
#include "../../libs/releaser.h"

#include "../../solve.h"

#include "../../grid/get_mean_max_two_fields.h"

#include "solve_mond_mass_multi.h"

#include "go_up_and_down.h"
#include "solve_mond_mass_multi_fine.h"

void do_one_v_cycle(struct grids_gravity_mm *gg, double a);
void cycle_v_mm_recursive(struct grids_gravity_mm *gg, double a);
void verbose_fine(struct grids_gravity_mm *gg, double a, int write_pot);


//====================
// cycle_v
//====================
int cycle_v_mm(struct grids_gravity_mm *gg, double a)
{
  // Cycle V routine.  Do iterations until convergence.  
  // See Brent 1977, Trottenberg 2000 or Wesseling 1992.
  // 
  // Arguments:
  //             a = expansion factor.
  //
  // Comments:
  //

  int iter;  // to test convergence
  double res_conv, cur_res;   // residuals

  double mean, maximum;

  // Show initial conditions:
  //=========================
  //gg->cur_resid[gg->g] = residual_mm(gg, a);   // not sure we need these two lines.
  //gg->cur_trunc[gg->g] = trunc_err_mm(gg, a);
  fprintf(stderr, "NGS = %d   verbose = %d\n", NGS, (gg->par)->verbose_mm);
  if((gg->par)->verbose_mm) verbose_fine(gg, a, 1);

  // Keep first residual for convergence test:
  //==========================================
  residual_mm(gg, gg->g, gg->res_phi[gg->g], gg->res_chi[gg->g]     );
  get_mean_max_two_fields(gg->res_phi[gg->g], gg->res_chi[gg->g], gg->n[gg->g]-1, &mean, &maximum);
  res_conv = mean;

  // Do iterations:
  //===============
  iter = 1;
  while(1)
    {
      
      // Do one V cycle:
      //================
      do_one_v_cycle(gg, a);

      // Get residual and trunc.err. for convergence:
      //=============================================
      residual_mm(gg, gg->g, gg->res_phi[gg->g], gg->res_chi[gg->g]     );
      get_mean_max_two_fields(gg->res_phi[gg->g], gg->res_chi[gg->g], gg->n[gg->g], &cur_res, &maximum);
      //gg->cur_resid[gg->g] = residual_mm(gg, a);
      //gg->cur_trunc[gg->g] = trunc_err_mm(gg, a);
      
      // Show status:
      //=============
      //fprintf(stderr, "Step %d    Residual = %le   Trunc = %le   Max_trunc = %le   Error = %le  level = %ld\n", 
      //	      iter, 
      //      gg->cur_resid[gg->g], 
      //      gg->cur_trunc[gg->g], 
      //      gg->max_trunc[gg->g], 
      //      gg->cur_resid[gg->g]/res_conv, 
      //      gg->g);
            
      // Check convergence:
      //===================
      //if( (gg->cur_resid[gg->g] < 1.0e-7 * res_conv) || 
      // 	  (gg->cur_resid[gg->g] < 1.0e-9) )
        if( //(gg->cur_resid[gg->g] < 0.1*gg->cur_trunc[gg->g]) || 
	  (cur_res < 1.0e-6 * res_conv) || 
	  (cur_res < 1.0e-9) )
	{
	  //fprintf(stderr, "cycle_v_mm:  converged.  residual_mm = %le   trunc = %le   max_trunc = %le   ", cur_res);   //, gg->cur_trunc[gg->g], gg->max_trunc[gg->g]);
	  fprintf(stderr, "cycle_v_mm:  converged.  residual_mm = %le   ", cur_res);   //, gg->cur_trunc[gg->g], gg->max_trunc[gg->g]);
	  if(cur_res < 1.0e-9)
	    fprintf(stderr, "(absolute criteria)\n");
	  else
	    fprintf(stderr, "\n");
	  //fprintf(stderr, "Rel. = %g   Abs = %g\n", 
	  //	  gg->cur_resid[gg->g] - 1.0e-7 * res_conv,
	  //	  gg->cur_resid[gg->g] - 1.0e-9);
	  return 1;
	}
      
      if(iter==MAX_ITERATION)
	{
	  fprintf(stderr, "WARNING: solve_mm_multigrid did not converged after %d iterations.\n", iter);
	  if((gg->par)->verbose_mm) verbose_fine(gg, a, 0);

	  return 0;
	  }
      //if(iter==20)
      //	do_multi = 0;
      
      iter++;
      
    }

}

//=============================
// do_one_v_cycle
//=============================
void do_one_v_cycle(struct grids_gravity_mm *gg, double a)
{

  // For choosing if we want two grids, three grids of full multigrid:
  int go_up_first  = 0;  // default = 1
  int go_up_second = 0;  // detault = 1

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
  if( (gg->g==0 && go_up_first) || go_up_second )
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



