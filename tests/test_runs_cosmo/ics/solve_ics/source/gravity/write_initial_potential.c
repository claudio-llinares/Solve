#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../solve.h"
#include "../grid/init_grid.h"

#include "../libs/D_lib.h"
#include "../libs/cosmology.h"

//==============================
// write_initial_potential
//==============================
void write_initial_potential(struct params *par, struct units *unit, struct grids *grid, int num)
{

  // Output potential at its time derivative.
  // The time derivative is not phi_dot, but
  // q = phi_dot*a**5

  char file[1000];
  FILE *fp_out;
  double out[10];

  double D_phi, D_q, norm_D, H;
  double D_prime;
  double move_phi, move_phi_dot;

  int i, j, k;

  double da_over_two;
  
  fprintf(stderr, "Writing fields...\n");

  da_over_two = (1.0-par->a_end)/par->nsteps_nbody/2.0;

  sprintf(file, "%s.%05d.potential", par->file_out, num);
  if( (fp_out = fopen(file, "w")) == NULL)
    {
      fprintf(stderr, "Problems opening the output file\n");
      exit(0);
    }
  D_phi = D_a(par->a_end, par->omegam, par->omegar, par->omegal, par->H0, &norm_D);
  D_q   = D_a(par->a_end+da_over_two, par->omegam, par->omegar, par->omegal, par->H0, &norm_D);
  H = hubble(par, unit, par->a_end+da_over_two);
  move_phi = D_phi/(par->a_end);

  D_prime = 1.0/(par->omegam + pow3(par->a_end)*par->omegal) * 
    ( -1.5*par->omegam/par->a_end*D_phi + 1.0/(pow2(par->H0)*norm_D) );
  move_phi_dot = 
    ( par->a_end*D_prime - D_phi ) * 
    par->H0 * 
    sqrt(par->omegam/pow3(par->a_end)+par->omegal);

  fprintf(stderr, "norm_D       = %lg\n", norm_D);
  fprintf(stderr, "move_phi     = %lg\n", move_phi);
  fprintf(stderr, "move_phi_dot = %lg\n", move_phi_dot);
  fprintf(stderr, "D_q          = %lg\n", D_q);

  for(i=0; i<par->grid; i++)  // x directon
    for(j=0; j<par->grid; j++)  // x directon
      for(k=0; k<par->grid; k++)  // x directon
	{
	  out[0] = move_phi     * grid->phi[i][j][k];
	  out[1] = move_phi_dot * grid->phi[i][j][k];
	  out[2] = D_phi * grid->delta[i][j][k];
	  fwrite(out, sizeof(double), 3, fp_out);
	  //fprintf(fp_out, "%lg  %lg  %lg  %lg  %lg  %lg\n",	
	  //	      out[0], out[1], out[2], out[3], out[4], out[5]);
	}
  fclose(fp_out);

}
