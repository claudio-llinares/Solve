/**
 * This is a static test.  If the corresponding flag is on in the parameter file, then the code will read initial conditions, calculate all the relevant fields, write them in a file and stops execution.
 * Comments:
 *      - There is a comparison between solutions obtained here and analytic solutions in the directory tests.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "../solve.h"

#include "../pm/init_parts.h"
#include "../grid/init_grid.h"
#include "../init_physics.h"

#include "../libs/cosmology.h"
#include "../libs/utils.h"

#include "../pm/calc_delta_cic.h"
#include "../gravity/get_forces.h"
#include "../gravity/newton/calc_newt_source.h"
#include "../gravity/newton/fft_solver/solve_poisson_fftw.h"

#include "../gravity/mond_mass/initial_guess_mm.h"
#include "../gravity/mond_mass/solve_multigrid_mm.h"


//=====================================
// static_test
//=====================================
void static_test(struct params *par, struct grids *grid, struct particles *part, struct units *unit)
{

  int i, j, k;
  double h, x_node, y_node, z_node, box_over_two;
  double out[10];
  FILE *fp;   //, *fp_ascii;
  char file[1000];

  h            = par->box/par->grid;
  box_over_two = par->box/2.0;

  // Calculate overdensity:
  //=======================
  calc_delta_cic(par, unit, part->pos_1, part->pos_2, part->pos_3, grid->delta);

  // Calculate potential and forces:
  //================================
  calc_newt_source(par, grid, unit, 1.0/(1.0+par->z_init));  // puts source in grid->source
  solve_poisson_fftw(par, grid, unit);
  forces(par, grid->phi, grid->force_1, grid->force_2, grid->force_3);

  // Calculate potential and forces MOND with mass (mm):
  //====================================================
  initial_guess_mm(par, grid->chi_bar_mm, grid->phi_bar_mm);
  solve_multigrid_mm(par, grid->delta, grid->phi_bar_mm, grid->chi_bar_mm, par->box, par->grid, par->a);
  forces(par, grid->phi_bar_mm, grid->force_phi_bar_mm_1, grid->force_phi_bar_mm_2, grid->force_phi_bar_mm_3);
  forces(par, grid->chi_bar_mm, grid->force_chi_bar_mm_1, grid->force_chi_bar_mm_2, grid->force_chi_bar_mm_3);

  // Write potential and forces:
  //============================
  fprintf(stderr, "%s\n", par->file_out);
  sprintf(file, "%s_static_test.dat", par->file_out);
  if( (fp = fopen(file, "w")) == NULL)
    {
      fprintf(stderr, "Problems opening the output file %s\n", file);
      exit(0);
    }
  //===========
  //sprintf(file, "%s_static_test_ascii.dat", par->file_out);
  //if( (fp_ascii = fopen(file, "w")) == NULL)
  //  {
  //    fprintf(stderr, "Problems opening the output file %s\n", file);
  //    exit(0);
  //  }
  int grid_over_two = par->grid/2;
  for(i=0; i<par->grid; i++)  // x directon
    for(j=0; j<par->grid; j++)  // y direction
      for(k=0; k<par->grid; k++)  // z direction
	{
	  if( !(  (i==grid_over_two && j==grid_over_two) || \
		  (i==grid_over_two && k==grid_over_two) || \
		  (j==grid_over_two && k==grid_over_two)  )  )
	    continue;

	  x_node = i*h + h/2.0;
	  y_node = j*h + h/2.0;
	  z_node = k*h + h/2.0;

	  out[0] = dist(x_node-box_over_two, y_node-box_over_two, z_node-box_over_two);
	  out[1] = grid->delta[i][j][k];
	  out[2] = grid->phi[i][j][k];
	  out[3] = dist(grid->force_1[i][j][k], grid->force_2[i][j][k], grid->force_3[i][j][k]);
	  out[4] = grid->chi_bar_mm[i][j][k];
	  out[5] = dist(grid->force_chi_bar_mm_1[i][j][k], grid->force_chi_bar_mm_2[i][j][k], grid->force_chi_bar_mm_3[i][j][k]);
	  out[6] = grid->phi_bar_mm[i][j][k];
	  out[7] = dist(grid->force_phi_bar_mm_1[i][j][k], grid->force_phi_bar_mm_2[i][j][k], grid->force_phi_bar_mm_3[i][j][k]);
	  fwrite(out, sizeof(double), 8, fp);
	  //fpintf(fp_ascii, "%f  %f  %f  %f\n", out[0], out[1], out[2], out[3]);
	  //fprintf(stderr, "%f  %f\n", x_node, out[0]);
	}
  fclose(fp);
  //fclose(fp_ascii);

  fprintf(stderr, "Done with static test.  We stop execution here.\n");
  exit(0);

}

