#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

#include "solve.h"

#include "libs/read_params.h"
#include "pm/init_parts.h"
#include "grid/init_grid.h"
#include "init_physics.h"
#include "grid/clean_up_grid.h"

#include "pm/write_particles.h"
#include "pm/zeldovich.h"

// Solves 2D growth equation.
//
// Output:
//    In output file:  x, y, delta
//
// Units:
//        input:
//                [x] = mpc/h
//                [rho] = g cm^-2
//        internal:
//                [x] = mpc/h
//                [source] = dimensionless
//                [phi] = mpc^2 h^-2 (takes units from distances)
// Comments:
//          - The derivative of density is NOT stored in the snapshot files.
//            The program can NOT restart.  The option for using an expansion
//            factor taken from the spanshot is for static calculations
//            of histograms and stuff in the first step.
//

int main(int argc, char *argv[])
{

  struct params par;
  struct grids grid;
  struct particles part;
  struct units unit;

  FILE *fp_in;

  int i, j, k;

  double in[10];
  double dummy_d;

  double sigma_delta, max_delta, max_delta_proc;

  // Initializes everything:
  //========================
      fprintf(stderr, "read_params...\n");
  read_params(&par);
      fprintf(stderr, "init_grid...\n");
  init_grid(&par, &grid);
      fprintf(stderr, "init_physics...\n");
  init_physics(&par, &unit);
      fprintf(stderr, "init_parts...\n");
  init_parts(&par, &part, &unit);
  
  // More initialization (we need a):
  //=================================
  par.a_init = 1.0/(1.0+par.z_init);
  par.a_end = 1.0/(1.0+par.z_end);
  
  // Read initial density (only if needed):
  //=======================================
  fprintf(stderr, "Reading density...\n");
  if( (fp_in = fopen(par.file_in, "r")) == NULL)
    {
      fprintf(stderr, "Problems opening the input file %s\n", par.file_in);
      exit(0);
    }
  // Read initial expansion factor:
  //fread(&dummy_d, sizeof(double), 1, fp_in);
  // Read the rest of the file:
  double mean_delta = 0.0;
  for(i=0; i<par.grid; i++)
    for(j=0; j<par.grid; j++)
      for(k=0; k<par.grid; k++)
	{
	  fread(in, sizeof(double), 1, fp_in);
	  grid.delta[i][j][k] = par.sign_delta * in[0];
	  mean_delta += grid.delta[i][j][k];
	}
  fclose(fp_in);
  fprintf(stderr, "Mean delta = %lg\n", mean_delta);
  
  // Zeldovich:
  //===========
  fprintf(stderr, "Calling zeldovich...\n");
  zeldovich(&par, &unit, &grid, &part, part.kernel);
 
  // Output the ics:
  //================
  write_particles(&par, &unit, part.pos_1, part.pos_2, part.pos_3, part.vel_1, part.vel_2, part.vel_3, 1);
  write_initial_potential(&par, &unit, &grid, 1);
  
  // Get maximun of overdensity:
  //============================
  sigma_delta = 0.0;
  max_delta = 0.0;
#ifdef OPENMP
#pragma omp parallel reduction(+:sigma_delta) private(i, j, k, max_delta_proc) shared(grid, max_delta, par)
  {
    max_delta_proc = 0.0;
#pragma omp for schedule(static)
#endif
    for(i=0; i<par.grid; i++)
      for(j=0; j<par.grid; j++)
	for(k=0; k<par.grid; k++)
	  {
	    sigma_delta += pow2(grid.delta[i][j][k]);
	    if(max_delta_proc < grid.delta[i][j][k])
	      max_delta_proc = grid.delta[i][j][k];
	  }
    if(max_delta_proc>max_delta)
      {
#pragma omp critical
	{
	  if( max_delta_proc > max_delta)
	    max_delta = max_delta_proc;
	}
      }
#ifdef OPENMP
  }
#endif
  fprintf(stderr, "delta_main: sigma = %le  max = %le  a = %le  z = %le\n", 
	  sqrt(sigma_delta/(pow3(par.grid)-1.0)), max_delta,
	  par.a_init, 1.0/par.a_init - 1.0);

  // Clean up:
  //==========
  clean_up_grid(&par, &grid);

  fprintf(stderr, "Finish with everything ...\n");
 
  return 0;
  
}
