#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include <fftw3.h>
#include "libs/ran3.h"

#include "p.h"

#include "libs/power_lib_fftw_corrected.h"  // for testing pourposes.
#include "libs/D_lib.h"  // for testing pourposes.
#include "libs/sigma8_lib.h"  // for testing pourposes.

#include "p2rho_lib.h"

extern float gasdev(long *idum);

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow6(x) ((x)*(x)*(x)*(x)*(x)*(x))

#define index(i,j,k) ((i) + par.grid * ((j) + par.grid * (k)))


//================================
// p2rho_lib
//================================
void p2rho_lib(double ***gr, double ***gi, struct param_p2rho par, fftw_complex *data, fftw_plan p_forward, fftw_plan p_backward)
{

  // Generate an overdensity realization of a power spectrum
  //
  // Arguments:
  //            - gr, gi = space for the grid (real and imaginary part).
  //                       (it should be allocated before).
  //            - par = structure with many parameters.
  //            - data, p_forward = space for fftw and plan.
  //                       (it should be allocated outside).
  // Returns:
  //            - gr, gi = overdenty (real and imaginary parts).
  //
  // Comments:
  //            - The power spectrum is taken from a file.  
  //              Parameters for that are given in par.
  //             - gr, gi must be allocate with alloc_double_3.

  // Get a realization of spectrum:
  //===============================
  fprintf(stderr, "Generating a realization of P...\n");
  spectr(gr, gi, par);
  
  // Gets delta:
  //============
  fprintf(stderr, "Calculating overdensity...\n");
  get_delta(gr, gi, data, p_forward, par);

}


//==========================
// spectr
//==========================
void spectr(double ***gr, double ***gi, struct param_p2rho par)
{
  // 
  // Makes a realization of the power spectrum given by the funciton P
  // 
  // Arguments:
  //            - gr = real part of the overdensity in fourier.
  //            - gi = real part of the overdensity in fourier.
  //             
  // Returns:
  //            - The function evaluated in points
  //                   f_{(n,m,o)} = 1.0/(par.grid*Delta) * (n,m,o)
  //                   n = -N/2, ... , N/2
  //                   m = -N/2, ... , N/2
  //                   o = -N/2, ... , N/2
  //
  // Comments:
  //            - You should allocate before with a grid given by the
  //              number of particles in each direction.
  //            - See why this:
  //                      if(i3 + j3 + k3 == 0)
  //                      {
  //                          grx[1][1][1] = 0.0;
  //                          gry[1][1][1] = 0.0;
  //                          grz[1][1][1] = 0.0;
  //                      }

  int m, n, o;
  int ix, iy, iz;
  
  double f2, f;
  
  double ar, ai;  // real and imaginary part of P(k)/k^2
  
  double wave;  // wave number for the P.

  int z_start, z_end; // begining of loop
  int go2, g;   // par.grid/2
  
  double norm;    // To go from P to a Fourier transform.
                  // The inverse of what we have in power_lib_fftw_corrected.c

  // To show the power specturm
  FILE *fp_show;
  int i, j, k;
  double Delta;
  char file_name[1000];
  double f_max;
  double d2;  // Newtonian growth factor at z=0 normalized to 1 at z=1000

  double kl;
  double p[1000];
  double phase[1000];

  double pi;

  double sqrt2;
  
  // Some optimization:
  double sqrt_P;

  pi = 4.0*atan(1.0);

  fprintf(stderr, "Generating realization of the power spectrum\n");

  // Real density:
  //==============
  // Get frecuencies with kz>0
  // looks like par.grid/2 has frecuency 0.
  norm = pow6(par.grid) / pow3(par.box);
  norm = 1.0;
  z_start = par.grid/2+1;
  z_end = par.grid-2;
  go2 = par.grid/2;
  g = par.grid;
  sqrt2 = sqrt(2.0);
  m = -par.grid/2;
  for(ix=0; ix<par.grid; ix++)
    {
      n = -par.grid/2;
      for(iy=0; iy<par.grid; iy++)
	{
	  o = -par.grid/2;
	  //for(iz=par.grid/2+1; iz<=par.grid-1; iz++)
	  for(iz=0; iz<par.grid; iz++)
	    {
	      f2 = (double)pow2(m) + (double)pow2(n) + (double)pow2(o);
	      f = sqrt(f2);
	      
	      if(f2==0)
		{
		  ar = 0.0;
		  ai = 0.0;
		}
	      else
		{
		  //WHY SQRT? (because P=<|delta|^2>)
		  wave = 2.0*pi*f / (par.grid*par.Delta);
		  sqrt_P = sqrt(P(wave, par.omega_0, par.hsmall, par.prefix, par.model, par.box, par.grid));
		  //ar = sqrt_P * gasdev(&par.Nseed);
		  //ai = sqrt_P * gasdev(&par.Nseed);
		  // The sqrt(2) distributes the power in real and imaginary.
		  //gr[ix][iy][iz] = norm * (double)ar / sqrt2;
		  //gi[ix][iy][iz] = norm * (double)ai / sqrt2;

		  // Optimized:
		  gr[ix][iy][iz] = norm * (double)(sqrt_P * gasdev(&par.Nseed)) / sqrt2;
		  gi[ix][iy][iz] = norm * (double)(sqrt_P * gasdev(&par.Nseed)) / sqrt2;

		  //fprintf(stderr, "A %d  %d  %d  %lg  %lg  %d  %d\n", 
		  //	  ix, iy, iz,
		  //	  gr[ix][iy][iz], gi[ix][iy][iz], m, n);
		}
	      ++o;
	    }
	  ++n;
	}
      ++m;
    }
  
 /*  // Enforce hermitian (Wagner) */
/*   int nby2,nby2p1,nby2p1xn; */
/*   int i,j,k; */
  
/*   nby2 = n/2; */
/*   nby2p1 = nby2 + 1; */
/*   nby2p1xn = nby2p1 * n; */
  
/*   // do k=0 and k=Nyquist planes */
/*   for (k=0;k<nby2+1;k+=nby2)  */
/*     {        */
/*       array[index(0,0,k)] = sqrt2*creal(array[index(0,0,k)]); */
/*       array[index(0,nby2,k)] = sqrt2*creal(array[index(0,nby2,k)]); */
/*       array[index(nby2,0,k)] = sqrt2*creal(array[index(nby2,0,k)]); */
/*       array[index(nby2,nby2,k)] = sqrt2*creal(array[index(nby2,nby2,k)]); */
      
/*       for (i=1;i<nby2;i++)  */
/* 	{ */
/* 	  array[index(0,n-i,k)] = conj(array[index(0,i,k)]); */
/* 	  array[index(nby2,n-i,k)] = conj(array[index(nby2,i,k)]); */
/* 	  array[index(n-i,0,k)] = conj(array[index(i,0,k)]); */
/* 	  array[index(n-i,nby2,k)] = conj(array[index(i,nby2,k)]); */
/* 	} */
      
/*       for (j=1;j<nby2;j++)  */
/* 	{ */
/* 	  for (i=1;i<nby2;i++) */
/* 	    array[index(n-j,n-i,k)] = conj(array[index(j,i,k)]); */
/* 	  for (i=nby2+1;i<n;i++) */
/* 	    array[index(n-j,n-i,k)] = conj(array[index(j,i,k)]); */
/* 	} */
/*     } */


  // The other quadrants:
  //m = -par.grid/2;
  for(ix=1; ix<=par.grid-1; ix++)
    {
      m = -par.grid/2 + ix;
      //n = -par.grid/2;
      for(iy=1; iy<=par.grid-1; iy++)
	{
	  n = -par.grid/2 + iy;
	  //o = -par.grid/2;
	  // There was iz=0 in the original.  It gives a memory leack.
	  for(iz=1; iz<=par.grid/2; iz++)
	    {
	      o = -par.grid/2 + iz;
	      f2 = (double)pow2(m) + (double)pow2(n) + (double)pow2(o);
	      f = sqrt(f2);
	      
	     /*  if(f2) */
/* 		{ */
/* 		  if((ix==0      && iy==0     ) || */
/* 		     (ix==0      && iy==par.grid/2-1) || */
/* 		     (ix==par.grid/2-1 && iy==0     ) || */
/* 		     (ix==par.grid/2-1 && iy==par.grid/2-1)) */
/* 		    { */
/* 		      gr[ix][iy][iz] = sqrt2 * gr[ix][iy][iz]; */
/* 		      gi[ix][iy][iz] = 0.0; */
/* 		    } */
/* 		  else */
/* 		    { */
/* 		      gr[ix][iy][iz] =  gr[g-ix][g-iy][g-iz]; */
/* 		      gi[ix][iy][iz] = -gi[g-ix][g-iy][g-iz]; */
/* 		    } */
/* 		} */

	      // Original
	      if(f2)
		{
		  if(ix != 0 &&
		     iy != 0 )
		    {
		      gr[ix][iy][iz] =  gr[g-ix][g-iy][g-iz];
		      gi[ix][iy][iz] = -gi[g-ix][g-iy][g-iz];
		    }
		  else
		    {
		      gr[ix][iy][iz] = sqrt2*gr[ix][iy][iz];
		      gi[ix][iy][iz] = 0.0;
		    }
		}
	      ++o;
	    }
	  ++n;
	}
      ++m;
    }

  fprintf(stderr, "Finish with realization of the power spectrum\n");
  
  // Shows the spectrum:
  //====================
  if(par.show_spectr)
    {
      // Get D
      d2 = pow2(D_1000(1.0, par.omega_0, 0.0, par.omega_l));
      fprintf(stderr, "D_1000(a=1) = %lg\n", sqrt(d2));

      // Output
      Delta = par.box/par.grid;
      sprintf(file_name, "%s.spect", par.prefix);
      fp_show=fopen(file_name, "w");
      f_max = 1.0e+50;
      m = -par.grid/2;
      for(k = 0; k<=par.grid-1; k++)
	{
	  n = -par.grid/2;
	  for(j = 0; j<=par.grid-1; j++)
	    {
	      o = -par.grid/2;
	      for(i = 0; i<=par.grid-1; i++)  
		{
		  f = sqrt((double)pow2(m) + (double)pow2(n) + (double)pow2(o));
		  wave = 2.0*pi*f / (par.grid*par.Delta);
		  //wave = f / (par.grid*par.Delta);
		  if(f<f_max)
		    {
		      fprintf(fp_show, "%d  %d  %d  %g  %g  %g  %g  %g  %g\n", 
			      m,
			      n, 
			      o,
			      wave / par.hsmall, 
			      P(wave, par.omega_0, par.hsmall, par.prefix, par.model, par.box, par.grid) * pow3(par.hsmall),
			      d2 * P(wave, par.omega_0, par.hsmall, par.prefix, par.model, par.box, par.grid) * pow3(par.hsmall),
			      pow2(gr[i][j][k])+pow2(gi[i][j][k]),
			      gr[i][j][k],
			      gi[i][j][k]);
		      f_max = f;
		    }
		  //2.0*unit.pi*f/(grid*Delta), // why THIS?
		  ++o;
		}
	      ++n;
	    }
	  ++m;
	}
      fclose(fp_show);
    }

  // Shows sigma_8:
  //===============
  if(par.show_sigma_8)
    {
      fprintf(stderr, "Writing sigma_8 stuff...\n");
      sprintf(file_name, "%s.spect.sigma8", par.prefix);
      fp_show=fopen(file_name, "w");
      i = 0;
      for(kl=-2; kl<=2; kl+=.01)
	{
	  phase[i] = pow(10, kl) / par.hsmall;
	  p[i] = P(pow(10,kl), par.omega_0, par.hsmall, par.prefix, par.model, par.box, par.grid) * pow3(par.hsmall);
	  fprintf(fp_show, "%lg  %lg\n", phase[i], p[i]);
	  i++;
	}
      fclose(fp_show);
      fprintf(stderr, "sigma_8(z=1000) = %g\n", sigma8(phase, p, 6.0/.01, 8) );
      fprintf(stderr, "sigma_8(z=0)    = %g\n", d2 * sigma8(phase, p, 6.0/.01, 8) );
    }
      
}

//==========================
// get_delta
//==========================
void get_delta(double ***gr, double ***gi, fftw_complex *data, fftw_plan p_forward, struct param_p2rho par)
{

  // Calculates nabla(phi) given the spectrum generated by spectr.
  //
  // Returns:
  //             - Derivatives of the potential where before
  //               was the fourier transform.
  // Comments:
  //             - Strictly speaking, we should multiply by (-i) but
  //               as the numbers are gaussian, will be the same, so
  //               we don't do it.
   
  int ix, iy, iz;

  double norm;

  double pi;

  int grid_over_2;

  // some optimization:
  int grid;
  
  pi = 4.0*atan(1.0);

  grid = par.grid;
  grid_over_2 = par.grid/2;

  // Fill data according to fftw:
  //=============================
  fprintf(stderr, "Transfering to fftw format...\n");
#ifdef OPENMP
#pragma omp parallel private(ix, iy, iz) shared(grid_over_2, par, data, gr, gi) //lap_phi)
  {
#pragma omp for schedule(static)
#endif
  for(ix=0; ix<grid_over_2; ix++)
    {
      for(iy=0; iy<grid_over_2; iy++)
	{
	  for(iz=0; iz<grid_over_2; iz++)
	    {
	      // --- frecuencies:
	      data[index(par.grid/2+ix,par.grid/2+iy,par.grid/2+iz)][0] = (float)gr[ix][iy][iz];
	      data[index(par.grid/2+ix,par.grid/2+iy,par.grid/2+iz)][1] = (float)gi[ix][iy][iz];
	      // ++- frecuencies:
	      data[index(ix,iy,par.grid/2+iz)][0] = (float)gr[ix+par.grid/2][iy+par.grid/2][iz];
	      data[index(ix,iy,par.grid/2+iz)][1] = (float)gi[ix+par.grid/2][iy+par.grid/2][iz];
	      // -+- frecuencies:
	      data[index(par.grid/2+ix,iy,par.grid/2+iz)][0] = (float)gr[ix][iy+par.grid/2][iz]; 
	      data[index(par.grid/2+ix,iy,par.grid/2+iz)][1] = (float)gi[ix][iy+par.grid/2][iz]; 
	      // +-- frecuencies:
	      data[index(ix,iy+par.grid/2,par.grid/2+iz)][0] = (float)gr[ix+par.grid/2][iy][iz]; 
	      data[index(ix,iy+par.grid/2,par.grid/2+iz)][1] = (float)gi[ix+par.grid/2][iy][iz]; 
	      
	      // --+ frecuencies:
	      data[index(par.grid/2+ix,par.grid/2+iy,iz)][0] = (float)gr[ix][iy][iz+par.grid/2]; 
	      data[index(par.grid/2+ix,par.grid/2+iy,iz)][1] = (float)gi[ix][iy][iz+par.grid/2];
	      // +++ frecuencies:
	      data[index(ix,iy,iz)][0] = (float)gr[ix+par.grid/2][iy+par.grid/2][iz+par.grid/2];
	      data[index(ix,iy,iz)][1] = (float)gi[ix+par.grid/2][iy+par.grid/2][iz+par.grid/2];
	      // -++ frecuencies:
	      data[index(par.grid/2+ix,iy,iz)][0] = (float)gr[ix][iy+par.grid/2][iz+par.grid/2];
	      data[index(par.grid/2+ix,iy,iz)][1] = (float)gi[ix][iy+par.grid/2][iz+par.grid/2];
	      // +-+ frecuencies:
	      data[index(ix,iy+par.grid/2,iz)][0] = (float)gr[ix+par.grid/2][iy][iz+par.grid/2];
	      data[index(ix,iy+par.grid/2,iz)][1] = (float)gi[ix+par.grid/2][iy][iz+par.grid/2];
	    }
	}
    }
#ifdef OPENMP
  }
#endif

  // Do Fourier transform:
  //======================
  fprintf(stderr, "Antitransforming FFT...\n");
  fftw_execute(p_forward);

  // Copy into output vector g (and normalizes):
  //============================================
  fprintf(stderr, "Normalizing...\n");
  // OLD ONE
  // This one is compatible with the normalization used in power for 
  // the test that where in the first version of the linear chapter.
  norm = 1.0/pow3(par.Delta)/pow3(par.grid) * pow3(sqrt(2.0*pi));  // OLD ONE
  // NEW ONE
  // This one comes from making things ompatible with the test made
  // with grafic1
  norm = 1.0/pow(par.Delta*par.grid, 3.0/2);
#ifdef OPENMP
#pragma omp parallel private(ix, iy, iz) shared(grid, data, gr, gi, norm) //lap_phi)
  {
#pragma omp for schedule(static)
#endif
  for(ix=0; ix<grid; ix++)
    {
      for(iy=0; iy<grid; iy++)      
	{
	  for(iz=0; iz<grid; iz++)      
	    {
	      gr[ix][iy][iz] = (double)data[index(ix,iy,iz)][0] * norm;
	      gi[ix][iy][iz] = (double)data[index(ix,iy,iz)][1] * norm;
	    }
	}
    }
#ifdef OPENMP
  }
#endif
  
}

