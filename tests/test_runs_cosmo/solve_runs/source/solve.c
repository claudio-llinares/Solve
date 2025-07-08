/****************************************************************************
 * Here is the main function.
 *
 * Units:
 *        input:
 *                - [x] = Mpc/h
 *                - [v] = km sec^-1
 *        internal:
 *                - [x] = Mpc
 *                - [v] = Mpc/Gyr
 *                - [source] = 1/Gyr^2
 *                - [phi] = Mpc^2/Gyr^2 (it takes units from distances)
 * Comments:
 *            - The program can NOT restart.  
 *****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "solve.h"

#include "libs/read_params.h"
#include "pm/init_parts.h"
#include "grid/init_grid.h"
#include "init_physics.h"
#include "grid/clean_up_grid.h"
#include "libs/cosmology.h"

#include "tests/static_test.h"

// Time evolution
#include "time_evolution/time_evolution_standard_variable.h"

/**
 Main routine.  It initializes everything (physics, fine grid and particles) and gives control to the routine that does the time evolution. */
int main(int argc, char *argv[])
{
  
  struct params par;
  struct grids grid;
  struct particles part;
  struct units unit;

  time_t start, end;
  int delta_t, hour, min, sec;

  // Initializes everything:
  //========================
  read_params(&par);
  init_grid(&par, &grid);
  init_physics(&par, &unit);
  init_parts(&par, &part, &unit);
  
  // Do static test and exit:
  //=========================
  if (par.static_test==1)
    static_test(&par, &grid, &part, &unit);

  // Do time evolution:
  //===================
  time(&start);
  leap_frog_variable(&par, &grid, &part, &unit);
  time(&end);

  // Outputs total wall time:
  //=========================
  delta_t = difftime (end, start);
  hour = (int)floor(delta_t/3600.0);
  min = (int)floor( (delta_t - hour*3600.0)/60.0 );
  sec = (int)(delta_t - hour*3600.0 - min*60.0);
  fprintf(stderr, "Total time:     %02d:%02d:%02d\n", 
	  hour, min, sec);
  
  // Clean up:
  //==========
  clean_up_grid(&par, &grid);

  fprintf(stderr, "Finished with everything ...\n"); 
 
  return 0;
  
}

