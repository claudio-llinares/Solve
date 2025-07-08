#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../solve.h"

#include "../commons_multigrid/commons_multigrid.h"  // for fourth order Laplacian

#include "source_mond_mass.h"


//==========================
// l_fith
//==========================
double l_mm(struct params *par, struct grids_gravity *gg, int grid, double source, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1)
{

  double val;
  
  double nb_sum, nabla_phi;

  double phi;

  phi = gg->chi_bar_mm[grid][i][j][k];

  nb_sum = 
    gg->chi_bar_mm[grid][i  ][j  ][km1] + gg->chi_bar_mm[grid][i  ][j  ][kp1] +
    gg->chi_bar_mm[grid][i  ][jm1][k  ] + gg->chi_bar_mm[grid][i  ][jp1][k  ] +
    gg->chi_bar_mm[grid][im1][j  ][k  ] + gg->chi_bar_mm[grid][ip1][j  ][k  ];

  nabla_phi = (nb_sum - 6.0*phi) / pow2(gg->h[grid]);

  // This one has better convergence in the small void problem:
  val = nabla_phi -
    0.5*pow2(a/par->symm_range)*phi * 
    (-1.0 + pow2(phi) +
     pow3(par->symm_assb/a)*(gg->rho[grid][i][j][k]+1.0)) - source;

  //fprintf(stderr, "%g  %g  %g  %g  %g  %g  %d\n", phi, nb_sum, par->symm_range, par->symm_assb, gg->rho[grid][i][j][k], val, grid);

  return val;

}

//==========================
// l_prime_mm
//==========================
double l_prime_mm(struct params *par, double phi, double rho, double source, double h, double a)
{

  double val;

  val = -6.0/pow2(h) -
    0.5*pow2(a/par->symm_range) *
    (-1.0 + 3.0*pow2(phi) + pow3(par->symm_assb/a)*(rho+1.0));

  // Some test:
  //===========
  if(val==0)
    {
      fprintf(stderr, "l_prime_mm:  the value is 0.  Aborting...");
      exit(0);
    }

  return val;

}

//==========================
// l_fith_fourth
//==========================
double l_mm_fourth(struct params *par, struct grids_gravity *gg, int grid, double source, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1)
{

  // This works only for the fine grid!!!!!

  double val;
  
  double nabla_phi;

  int im2, ip2, jm2, jp2, km2, kp2;
  int n = gg->n[grid];
  double nb_sum;

  double phi = gg->chi_bar_mm[grid][i][j][k];

  ip2 = (i<(n-2))? i+2 : (i<(n-1)? 0 : 1);
  jp2 = (j<(n-2))? j+2 : (j<(n-1)? 0 : 1);
  kp2 = (k<(n-2))? k+2 : (k<(n-1)? 0 : 1);
  
  im2 = (i>1)? i-2 : (i>0? (n-1) : (n-2));
  jm2 = (j>1)? j-2 : (j>0? (n-1) : (n-2));
  km2 = (k>1)? k-2 : (k>0? (n-1) : (n-2));
  
  nb_sum = 
    -1.0/12*gg->chi_bar_mm[grid][im2][j][k] - 1.0/12.0*gg->chi_bar_mm[grid][ip2][j][k] + 
    4.0/3.0*gg->chi_bar_mm[grid][im1][j][k] +  4.0/3.0*gg->chi_bar_mm[grid][ip1][j][k];
  nb_sum += 
    -1.0/12*gg->chi_bar_mm[grid][i][jm2][k] - 1.0/12.0*gg->chi_bar_mm[grid][i][jp2][k] + 
    4.0/3.0*gg->chi_bar_mm[grid][i][jm1][k] +  4.0/3.0*gg->chi_bar_mm[grid][i][jp1][k];
  nb_sum += 
    -1.0/12*gg->chi_bar_mm[grid][i][j][km2] - 1.0/12.0*gg->chi_bar_mm[grid][i][j][kp2] + 
    4.0/3.0*gg->chi_bar_mm[grid][i][j][km1] +  4.0/3.0*gg->chi_bar_mm[grid][i][j][kp1];
  
  nabla_phi = (nb_sum - 3.0*5.0/2.0*gg->chi_bar_mm[grid][i][j][k] ) / 
    pow2(gg->h[grid]);

  // Wrong version, but convergent
  val = nabla_phi -
    0.5*pow2(a/par->symm_range)*phi * 
    (-1.0 + pow2(phi) +
     pow3(par->symm_assb/a)*(gg->rho[grid][i][j][k]+1.0)) - source;
  
  return val;

}

//==========================
// l_prime_mm_fourth
//==========================
double l_prime_mm_fourth(struct params *par, double phi, double rho, double source, double h, double a)
{

  double val;

  // Symmetron:
  //===========
  // Original:
  val = -3.0*5.0/2.0/pow2(h) -
    0.5*pow2(a/par->symm_range) *
    (-1.0 + 3.0*pow2(phi) + pow3(par->symm_assb/a)*(rho+1.0));
 
  // Some test:
  //===========
  if(val==0)
    {
      fprintf(stderr, "l_prime_mm:  the value is 0.  Aborting...");
      exit(0);
    }

  return val;

}


//======================================
// L_mm
//======================================
double L_mm(int i, int j, int k, double ***phi, double h, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int n)
{

  //
  // Differential operator (second order)
  //

  double nb_sum;
  double L_mm;

  nb_sum = 
    phi[i  ][j  ][km1] + phi[i  ][j  ][kp1] +
    phi[i  ][jm1][k  ] + phi[i  ][jp1][k  ] +
    phi[im1][j  ][k  ] + phi[ip1][j  ][k  ];

  L_mm = (nb_sum - 6.0*phi[i][j][k])/pow2(h);

  return L_mm;
  
}

//======================================
// L_mm_fourth
//======================================
double L_mm_fourth(int i, int j, int k, double ***phi, double h, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int n)
{

  //
  // Differential operator (fourth order)
  //

  double nb_sum;
  double L;

  int im2, ip2, jm2, jp2, km2, kp2;

  ip2 = (i<(n-2))? i+2 : (i<(n-1)? 0 : 1);
  jp2 = (j<(n-2))? j+2 : (j<(n-1)? 0 : 1);
  kp2 = (k<(n-2))? k+2 : (k<(n-1)? 0 : 1);
  
  im2 = (i>1)? i-2 : (i>0? (n-1) : (n-2));
  jm2 = (j>1)? j-2 : (j>0? (n-1) : (n-2));
  km2 = (k>1)? k-2 : (k>0? (n-1) : (n-2));

  nb_sum = 
    -1.0/12*phi[im2][j][k] - 1.0/12.0*phi[ip2][j][k] + 
    4.0/3.0*phi[im1][j][k] +  4.0/3.0*phi[ip1][j][k];
  nb_sum += 
    -1.0/12*phi[i][jm2][k] - 1.0/12.0*phi[i][jp2][k] + 
    4.0/3.0*phi[i][jm1][k] +  4.0/3.0*phi[i][jp1][k];
  nb_sum += 
    -1.0/12*phi[i][j][km2] - 1.0/12.0*phi[i][j][kp2] + 
    4.0/3.0*phi[i][j][km1] +  4.0/3.0*phi[i][j][kp1];
  
  // Nabla phi:
  L = (nb_sum - 3.0*5.0/2.0*phi[i][j][k] ) / 
    pow2(h);

  return L;
  
}


