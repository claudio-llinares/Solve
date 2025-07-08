#include <stdlib.h>
#include <stdio.h>

#include "../solve.h"
#include "../pm/init_parts.h"
//=========================
// write_particles
//=========================
void write_particles(struct params *par, struct units *unit, double *pos_1, double *pos_2, double *pos_3, double *vel_1, double *vel_2, double *vel_3, int num)
{

  char file[1000];
  FILE *fp_out;
  double out[10];

  int i;
  
  //sprintf(file, "%s.%05d.part", par->file_out, num);
  sprintf(file, "conv.0.0");
  if( (fp_out = fopen(file, "w")) == NULL)
    {
      fprintf(stderr, "Problems opening the output file\n");
      exit(0);
    }
  //fp_out_ascii = fopen("conv_ascii.0.0", "w");

  for(i=0; i<par->nparts; i++)  // x directon
    {
      out[0] = pos_1[i] * par->hsmall;
      out[1] = pos_2[i] * par->hsmall;
      out[2] = pos_3[i] * par->hsmall;
      out[3] = vel_1[i] * unit->mpc2km / unit->gyr2sec;
      out[4] = vel_2[i] * unit->mpc2km / unit->gyr2sec;
      out[5] = vel_3[i] * unit->mpc2km / unit->gyr2sec;
      out[6] = (double)i;
      fwrite(out, sizeof(double), 7, fp_out);
      //fprintf(fp_out_ascii, "%20.15lf  %20.15lf  %20.15lf\n", 
      //      out[0], out[1], out[2]);
      //fprintf(fp_out, "%lg  %lg  %lg  %lg  %lg  %lg\n",	
      //	      out[0], out[1], out[2], out[3], out[4], out[5]);
    }
  fclose(fp_out);

}


