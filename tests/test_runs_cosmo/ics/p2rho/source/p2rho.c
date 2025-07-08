#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include <fftw3.h>

#include "libs/allocator.h"
#include "libs/releaser.h"

#include "libs/ran3.h"

#include "libs/z2t_lib.h"

#include "libs/power_lib_fftw_corrected.h"  // for testing pourposes.
#include "libs/D_lib.h"  // for testing pourposes.
#include "libs/sigma8_lib.h"  // for testing pourposes.

#include "p2rho_lib.h"  // The core.

#define pow3(x) ((x)*(x)*(x))
#define pow2(x) ((x)*(x))

struct units
{
  double H0;  // units:  sec^-1 (real H0 = hsmall*100).
  double pi;
  double mpc2km;
  double mpch2int;
  double gyr2sec;
} unit;

int main(void)
{
  double ***gr, ***gi;  // realization of the spectrum and delta
  
  int i, j, k;

  double out[10];   // to output

  // fftw data (for testing)
  fftw_complex *data;
  fftw_plan p_forward;
  fftw_plan p_backward;
  
  FILE *fp_out;

  double pi;
  pi = 4.0*atan(1.0);

  struct param_p2rho par;

  // for testing power spectrum
  char file[1000];
  double *phase, *p;

  int por;
  
  // Read input parameters:
  //=======================
  por = scanf("%s", par.prefix);
  por = scanf("%lf", &par.box);
  por = scanf("%d", &par.grid);
  por = scanf("%lf", &par.omega_0);
  por = scanf("%ld", &par.Nseed);
  por = scanf("%lf", &par.hsmall);
  por = scanf("%d", &par.model);
  por = scanf("%d", &par.show_spectr);
  por = scanf("%d", &par.test);
  por = scanf("%d", &par.mean_delta);
  por = scanf("%d", &par.show_sigma_8);
  por = scanf("%d", &par.nthreads);

  fprintf(stderr, "par.prefix       = %s\n", par.prefix);
  fprintf(stderr, "par.box          = %lf\n", par.box);
  fprintf(stderr, "par.grid         = %d\n", par.grid);
  fprintf(stderr, "par.omega_0      = %lf\n", par.omega_0);
  fprintf(stderr, "par.Nseed        = %ld\n", par.Nseed);
  fprintf(stderr, "par.hsmall       = %lf\n", par.hsmall);
  fprintf(stderr, "par.model        = %d\n", par.model);
  fprintf(stderr, "par.show_spectr  = %d\n", par.show_spectr);
  fprintf(stderr, "par.test         = %d\n", par.test);
  fprintf(stderr, "par.mean_delta   = %d\n", par.mean_delta);
  fprintf(stderr, "par.show_sigma_8 = %d\n", par.show_sigma_8);
  fprintf(stderr, "par.nthreads = %d\n", par.nthreads);

  // Allocates stuff for test:
  //==========================
  p     = calloc(par.grid, sizeof(double));
  phase = calloc(par.grid, sizeof(double));

  // Some cosmology/units:
  //======================
  unit.mpc2km  = 3.08568025e19;
  unit.mpch2int = 1.0 / par.box;
  //par.box /= par.hsmall;   // scale it to real megaparsecs
  unit.H0 = par.hsmall * 100.0 / unit.mpc2km;  // sec^-1
  unit.gyr2sec = 1.0e+9 * 365.25 * 23.56 * 3600.0;
  unit.pi = 4.0*atan(1.0);

  par.box /= par.hsmall;   // put things in Mpc.

  par.test_corr = 0;
  par.omega_l = 1.0 - par.omega_0;
  par.z = 1000.0;
  par.a = 1.0 / (1.0 + par.z);
  par.t = z2t(par.omega_0, par.omega_l, 100.0*par.hsmall, par.z) * 
    unit.gyr2sec;    // sec
  par.Delta = par.box / par.grid;
  fprintf(stderr, "z = %lf\n", par.z);
  fprintf(stderr, "a = %lf\n", par.a);
  fprintf(stderr, "t = %lf sec/h = %lf Gyr\n", par.t, par.t/unit.gyr2sec);

#ifdef OPENMP
  if( !fftw_init_threads())
    {
      fprintf(stderr, "Problems with fftw_init_threads.  Aborting...\n");
      exit(0);
    }
  else
    fprintf(stderr, "fftw threads initialized.\n");
#endif

  data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * par.grid * par.grid * par.grid);

#ifdef OPENMP
  fftw_plan_with_nthreads(par.nthreads);
#endif
  p_backward = fftw_plan_dft_3d(par.grid, par.grid, par.grid, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
#ifdef OPENMP
  fftw_plan_with_nthreads(par.nthreads);
#endif
  p_forward = fftw_plan_dft_3d(par.grid, par.grid, par.grid, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
  alloc_double_3(&gr, par.grid, par.grid, par.grid);
  alloc_double_3(&gi, par.grid, par.grid, par.grid);

  // Generate the overdensity realization:
  //======================================
  p2rho_lib(gr, gi, par, data, p_forward, p_backward);

  // Output:
  //========
  fprintf(stderr, "Writing...\n");
  if( (fp_out=fopen(par.prefix, "w")) == NULL)
    {
      fprintf(stderr, "Problems opening output file %s\n", par.prefix);
      exit(0);
    }
  // Write header
  //out[0] = -1;
  //fwrite(out, sizeof(double), 1, fp_out);
  // Write the rest of the file.
  for(i=0; i<=par.grid-1; i++)
    {
      for(j=0; j<=par.grid-1; j++)
	{
	  for(k=0; k<=par.grid-1; k++)
	    {
	      if(0)
		{
		  out[0] = ((double)i*par.Delta + par.Delta/2.0) * par.hsmall;
		  out[1] = ((double)j*par.Delta + par.Delta/2.0) * par.hsmall;
		  out[2] = ((double)k*par.Delta + par.Delta/2.0) * par.hsmall;
		  out[3] = gr[i][j][k];
		  out[4] = gi[i][j][k];
		  fwrite(out, sizeof(double), 5, fp_out);
		}
	      else
		{
		  out[0] = gr[i][j][k];
		  fwrite(out, sizeof(double), 1, fp_out);
		}
	      //fprintf(fp_out, "%lg  %lg  %lg  %g  %g\n",
	      //	      	      i*par.Delta + par.Delta/2.0, 
	      //	      	      j*par.Delta + par.Delta/2.0, 
	      //	      	      k*par.Delta + par.Delta/2.0,
	      //	      	      gr[i][j][k], 
	      //	      	      gi[i][j][k]);
	    }
	}
    }
  fclose(fp_out);

  // Calculates and output P:
  //=========================
  if(par.test)
    {  

      if(par.mean_delta)
    {
      // Calculate mean of delta
      double mean, sigma2;
      double delta_min, delta_max;
      mean = 0.0;
      delta_min = 1.0e+50;
      delta_max = -1.0e+50;
      for(i=0; i<par.grid; i++)
	for(j=0; j<par.grid; j++)
	  for(k=0; k<par.grid; k++)
	    {
	      mean += gr[i][j][k];
	      // Calculate max and min density
	      if(gr[i][j][k]>delta_max)
		delta_max = gr[i][j][k];
	      if(gr[i][j][k]<delta_min)
		delta_min = gr[i][j][k];
	    }
      
      mean /= pow3(par.grid);

      // Calculate sigma:
      sigma2 = 0.0;
      for(i=0; i<par.grid; i++)
	for(j=0; j<par.grid; j++)
	  for(k=0; k<par.grid; k++)
	    sigma2 += pow2(gr[i][j][k] - mean);

      sigma2 /= pow3(par.grid)-1.0;
      
      // Output mean and sigma:
      fprintf(stderr, "Delta mean = %g\n", mean);
      fprintf(stderr, "Delta sigma = %g\n", sqrt(sigma2));
      fprintf(stderr, "Delta max = %g\n", delta_max);
      fprintf(stderr, "Delta min = %g\n", delta_min);
      fprintf(stderr, "min max mean sigma = %g  %g  %g  %g\n", delta_min, delta_max, mean, sqrt(sigma2));
    }
  
      // To test real densities:
      if(0)
	{
	  for(i=0; i<=par.grid-1; i++)
	    for(j=0; j<=par.grid-1; j++)
	      for(k=0; k<=par.grid-1; k++)
		gi[i][j][k] = 0.0;
	}
      
      // Calculates P
      // THIS CALL CHANGES gr AND gi!!!!
      power_corrected(gr, gi, par.box, par.grid, phase, p, 0, data, p_backward);
      
      // Output it
      sprintf(file, "%s.p_out", par.prefix);
      if( (fp_out=fopen(file, "w")) == NULL)
	{
	  fprintf(stderr, "Problems opening output file %s\n", par.prefix);
	  exit(0);
	}
      //for(k=1; k<=par.grid/2.0; k++)//(int)(sqrt(3.0)*par.grid+.5); k++)
      
      for(k=1; k<par.grid; k++)//(int)(sqrt(3.0)*par.grid+.5); k++)
	fprintf(fp_out, "%lg  %lg\n", 
		phase[k] / par.hsmall, 
		p[k] * pow3(par.hsmall));
      fclose(fp_out);

      free(phase);
      free(p);

    }

  fprintf(stderr, "Finish with everything.\n");

  return 0;
  
}

