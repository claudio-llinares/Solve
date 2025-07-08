#include <stdio.h>
#include <stdlib.h>

#include "../solve.h"

// Read and write time information for restarting a simulation.

//=====================
// read_times
//=====================
void read_times(struct params *par, int *l, int *l_out, double *a_out, double delta_a_out)
{

  FILE *fp;
  char file[1000];

  double z_init, z_end;
  double delta_a_out_ff;

  double in[10];
  int in_int[10];

  int por;

  // ff means from file.

  fprintf(stderr, "Reading times from file...\n");

  // Read stuff:
  //============
  sprintf(file, "%s.times", par->file_in_parts);
  if( (fp = fopen(file, "r")) == NULL)
    {
      fprintf(stderr, "Problems opening the input file %s\n", file);
      exit(0);
    }
  por = fread(in, 5, sizeof(double), fp);
  par->a         = in[0];
  z_init         = in[1];
  z_end          = in[2];
  *a_out         = in[3];
  delta_a_out_ff = in[4];
  
  por = fread(in_int, 2, sizeof(int), fp);
  *l          = in_int[0];
  *l_out      = in_int[1];
  
  // Some tests:
  //============
  if(delta_a_out)
    {
      if(delta_a_out_ff != delta_a_out)
	fprintf(stderr, "     Incompatible number of outputs...\n");
    }
  else
    fprintf(stderr, "read_times: not checking number of outputs...\n");
  if( (z_init!=par->z_init) || (z_end!=par->z_end) )
    fprintf(stderr, "     Incompatible initial or final redshifts...\n");
  

  fprintf(stderr, "Parameters from %s.times:\n", par->file_in_parts);
  fprintf(stderr, "  a           = %le\n", par->a);
  fprintf(stderr, "  a_out       = %le\n", *a_out);
  fprintf(stderr, "  l           = %d\n", *l);
  fprintf(stderr, "  l_out       = %d\n", *l_out);
  fprintf(stderr, "\n");
  
  fclose(fp);

}


//=====================
// write_times
//=====================
void write_times(struct params *par, int num, int l, int l_out, double a_out, double delta_a_out)
{

  FILE *fp;
  char file[1000];

  double out[10];
  int out_int[10];

  // ff means from file.

  fprintf(stderr, "writing times to file...\n");

  // Read stuff:
  //============
  sprintf(file, "%s.%05d.part.times", par->file_out, num);
  if( (fp = fopen(file, "w")) == NULL)
    {
      fprintf(stderr, "Problems opening the output file %s\n", file);
      exit(0);
    }
  out[0] = par->a;
  out[1] = par->z_init;
  out[2] = par->z_end;
  out[3] = a_out;
  out[4] = delta_a_out;
  fwrite(out, 5, sizeof(double), fp);
  
  out_int[0] = l;
  out_int[1] = l_out;
  fwrite(out_int, 2, sizeof(int), fp);
  
  fclose(fp);

}
