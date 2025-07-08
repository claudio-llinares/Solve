/**
 *  Here we read initial conditions for particles.
 * Comments:
 *      - Units in the input file are expected to be Mpc/h for the positions and km/sec for the velocities.
 *      - Units are converted later to Mpc and Mpc/Gyr.
 **/

#include <stdio.h>
#include <stdlib.h>

#include "../solve.h"
#include "../grid/init_grid.h"
#include "../libs/utils.h"

//==============================
// read_particles
//==============================
void read_particles(struct params *par, struct units *unit, double *pos_1, double *pos_2, double *pos_3, double *vel_1, double *vel_2, double *vel_3, char *file_name)
{

  int i;
  double out[10];
  FILE *fp_in;
  int fre;

  fprintf(stderr, "Reading particles...\n");

  // Read initial conditions for the particles:
  //===========================================
  if( (fp_in = fopen(file_name, "r")) == NULL)
    {
      fprintf(stderr, "Problems opening the input file %s\n", par->file_in_parts);
      exit(0);
    }
  // The order is three loop for x, y, z directions
  // (that the original for 3D arrays).
  for(i=0; i<pow3(par->grid_p); i++)
    {  
	  fre = fread(out, sizeof(double), 6, fp_in);
	  if(feof(fp_in))
	    {
	      fprintf(stderr, "Wrong number of particles.  Aborting.\n");
	      exit(0);
	    }
	  pos_1[i] = ( out[0] / par->hsmall ); // / par->box;
	  pos_2[i] = ( out[1] / par->hsmall ); // / par->box;
	  pos_3[i] = ( out[2] / par->hsmall ); // / par->box;
	  vel_1[i] = ( out[3] * unit->km2mpc / unit->sec2gyr ); // / par->clight;// mpc/gyr
	  vel_2[i] = ( out[4] * unit->km2mpc / unit->sec2gyr ); // / par->clight;
	  vel_3[i] = ( out[5] * unit->km2mpc / unit->sec2gyr ); // / par->clight;
	  
	  pos_1[i] = modulus(pos_1[i], par->box);
	  pos_2[i] = modulus(pos_2[i], par->box);
	  pos_3[i] = modulus(pos_3[i], par->box);
	}
  fclose(fp_in);

  fprintf(stderr, "done\n");

}
