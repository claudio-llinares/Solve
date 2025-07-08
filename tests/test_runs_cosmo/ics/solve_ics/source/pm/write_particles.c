#include <stdlib.h>
#include <stdio.h>

#include "../solve.h"

//=========================
// write_particles
//=========================
void write_particles(struct params *par, struct units *unit, double *pos_1, double *pos_2, double *pos_3, double *vel_1, double *vel_2, double *vel_3, int num)
{

  char file[1000];
  FILE *fp_out;
  double out[10];

  int i;
  
  sprintf(file, "%s.%05d.part", par->file_out, num);
  if( (fp_out = fopen(file, "w")) == NULL)
    {
      fprintf(stderr, "Problems opening the output file\n");
      exit(0);
    }
  for(i=0; i<par->nparts; i++)  // x directon
    {
      out[0] = pos_1[i] * par->hsmall;
      out[1] = pos_2[i] * par->hsmall;
      out[2] = pos_3[i] * par->hsmall;
      out[3] = vel_1[i] * unit->mpc2km / unit->gyr2sec;
      out[4] = vel_2[i] * unit->mpc2km / unit->gyr2sec;
      out[5] = vel_3[i] * unit->mpc2km / unit->gyr2sec;
      fwrite(out, sizeof(double), 6, fp_out);
      //fprintf(fp_out, "%lg  %lg  %lg  %lg  %lg  %lg\n",	
      //	      out[0], out[1], out[2], out[3], out[4], out[5]);
    }
  fclose(fp_out);

}
