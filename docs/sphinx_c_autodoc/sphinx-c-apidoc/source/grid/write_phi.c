#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

#include "../solve.h"
#include "init_grid.h"

void write_phi(struct params *par, struct grids *grid)
{

  FILE *fp;
  char file[1000];

  int i, j, k;

  double out[10];

  static int n_file = 1;

  //double h = par->box/par->grid;
  double delta_max;
  int k_max;

  fprintf(stderr, "Entering write_newton...\n");

  // Find z of maximum delta:
  //=========================
  if(0)
    {
      delta_max = -10.0;
      k_max = 0;
      for(i=0; i<par->grid; i++)
	for(j=0; j<par->grid; j++)
	  for(k=0; k<par->grid; k++)  // z direction
	    {
	      if(grid->delta[i][j][k] > delta_max)
		k_max = k;
	    }
      k = k_max;
    }
  else
    k = par->grid/2;

  // Write stuff:
  //=============
  fprintf(stderr, "%s\n", par->file_out);
  sprintf(file, "%s.phi.%05d", par->file_out, n_file);
  n_file++;
  if( (fp = fopen(file, "w")) == NULL)
    {
      fprintf(stderr, "Problems opening the output file %s\n", file);
      exit(0);
    }
  for(i=0; i<par->grid; i++)  // x directon
    {
      for(j=0; j<par->grid; j++)  // y direction
	{
	  out[0] = grid->delta[i][j][k];
	  out[1] = grid->phi[i][j][k];
	  out[2] = sqrt( pow2(grid->force_1[i][j][k]) + pow2(grid->force_2[i][j][k]) + pow2(grid->force_3[i][j][k]) );
	  fwrite(out, sizeof(double), 3, fp);
	}
    }
  fclose(fp);

}
