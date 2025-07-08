/**
 * Here is the main time loop.
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef OPENMP
#include <omp.h>
#endif
#include <sys/time.h>

#include "../solve.h"
#include "../pm/init_parts.h"
#include "../grid/init_grid.h"

#include "../libs/cosmology.h"
#include "../libs/utils.h"
#include "../libs/allocator.h"

#include "../pm/interp_force_cic.h"
#include "../pm/write_particles.h"

#include "../gravity/get_gravity.h"

#include "../grid/write_phi.h"


//=========================
// leap_frog
//=========================
void leap_frog_variable(struct params *par, struct grids *grid, struct particles *part, struct units *unit)
{

  double move_part;
  
  int i;
  
  double t, t_init, t_end, dt;
  double dt1; //, dt2;
  double a, da;
  double t_of_a_plus_da;
  int l;
  
  int n_parts_written = 0;

  struct timeval start, end;
  struct timeval start_step, end_step;
  double elapsed, elapsed_step;
  int hour, min, sec;

#ifdef OPENMP
#pragma omp parallel 
  {
    if(!omp_get_thread_num())
      fprintf(stderr, "leap_frog_variable: using %d threads\n", omp_get_num_threads());
  }
#endif

  fprintf(stderr, "Entering leap_frog...\n");

  // Initialize some stuff:
  //=======================
  t_init = integrate_time_of_a(par, 1.0/(1.0+par->z_init));
  t_end = integrate_time_of_a(par, 1.0/(1.0+par->z_end));
  t = t_init;
  a = 1.0/(1.0+par->z_init);
  da = (1.0/(1.0+par->z_end) - 1.0/(1.0+par->z_init))/par->nsteps;
  l = 0;
  
  fprintf(stderr, "da = %10.7lf\n", da);
      
  // Calculate potential and forces:
  //================================
  if(1)
    {
      get_gravity(par, grid, part, unit, t);
    }
  
  // Leap-frog:
  //===========
  fprintf(stderr, "Starting leap-frog:\n");
  fprintf(stderr, "===================\n");
  gettimeofday(&start, 0);
  gettimeofday(&start_step, 0);
  while(t<t_end)
    {
      // Meassure time:
      //===============
      gettimeofday(&end, 0);
      gettimeofday(&end_step, 0);
      // elapsed time since begining in sec
      elapsed = (double) ( ( (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec ) / 1000000.0 );
      hour = (int)floor(elapsed/3600.0);
      min = (int)floor( (elapsed - hour*3600.0)/60.0 );
      sec = (int)(elapsed - hour*3600.0 - min*60.0);
      elapsed_step = (double)( ( (end_step.tv_sec - start_step.tv_sec)*1000000 + end_step.tv_usec - start_step.tv_usec ) / 1000000.0 );
      gettimeofday(&start_step, 0);
      
      // Some verbose:
      //==============
      if(par->verbose)
	write_phi(par, grid);
      
      // Calculate dt:
      //==============
      t_of_a_plus_da = integrate_time_of_a(par, a+da);
      dt1 = t_of_a_plus_da - t;  // for uniform steps in a
      //dt2 = (t_end-t_init)/par->nsteps;
      dt = dt1;

      fprintf(stderr, "====================================================================================================\n");
      fprintf(stderr, "      n       t            a            z            dt           Dt(step)     Dt(total)\n");
      fprintf(stderr, "Step  %-5d   %-10.7lf   %-10.7lf   %-10.7lf   %-10.7lf   %-10lg   %02d:%02d:%02d\n", 
	      l, t, a_of_t(par,t), 1.0/a_of_t(par,t) - 1.0, dt, 
	      elapsed_step, 
	      hour, 
	      min, 
	      sec);

      // Advance positions (first):
      //===========================
      if(1)
	{
	  move_part =  (dt/2.0)/a_of_t(par, t);
#ifdef OPENMP
#pragma omp parallel private(i) shared(par, part, move_part)
	  {
#pragma omp for schedule(static)
#endif
	    for(i=0; i<par->nparts; i++)
	      {
		part->pos_1[i] = modulus( part->pos_1[i] + part->vel_1[i]*move_part, par->box);
		part->pos_2[i] = modulus( part->pos_2[i] + part->vel_2[i]*move_part, par->box);
		part->pos_3[i] = modulus( part->pos_3[i] + part->vel_3[i]*move_part, par->box);
	      }
#ifdef OPENMP
	  }
#endif
	}
      
      // Calculate delta, and forces:
      //=============================
      if(1)
	{
	  get_gravity(par, grid, part, unit, t+dt/2.0);
	}

      // Advance velocities:
      //====================
      if(1)
	{
	  move_part = a_of_t(par, t+dt/2.0) * (dt);
	  interp_force_cic(par, grid, part, grid->force_1, grid->force_2, grid->force_3, move_part);
	}

      // Advance positions (second):
      //============================
      if(1)
	{
	  move_part =  (dt/2.0)/a_of_t(par, t+dt);
#ifdef OPENMP
#pragma omp parallel private(i) shared(par, part, move_part)
	  {
#pragma omp for schedule(static)
#endif
	    for(i=0; i<par->nparts; i++)
	      {
		part->pos_1[i] = modulus( part->pos_1[i] + part->vel_1[i]*move_part, par->box);
		part->pos_2[i] = modulus( part->pos_2[i] + part->vel_2[i]*move_part, par->box);
		part->pos_3[i] = modulus( part->pos_3[i] + part->vel_3[i]*move_part, par->box);
	      }
#ifdef OPENMP
	  }
#endif
	}

      // Advance time:
      //==============
      t += dt;
      a += da;
      if(t>t_end)
	t = t_end;
      l++;
            
    }  // the time loop
  

  // Write final state:
  //===================
  write_phi(par, grid);
  write_particles(par, unit, 
		  part->pos_1, part->pos_2, part->pos_3, 
		  part->vel_1, part->vel_2, part->vel_3, -1);
  

}
