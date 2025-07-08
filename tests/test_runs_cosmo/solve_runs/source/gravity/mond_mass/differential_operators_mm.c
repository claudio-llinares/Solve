/**
 * Definition of the differential operator for both equations (Poisson like and MOND like).
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"
#include "../../libs/utils.h"

#include "solve_multigrid_mm.h"


double mu(double x);


/** Evaluates differential operators on a node **/
void L_mm(double ***phi, double ***chi, double aa0, double a2M2, double h, double fourh, int i, int j, int k, int ip1, int im1, int jp1, int jm1, int kp1, int km1, double *operator_phi, double *operator_chi)
{

  //
  // Differential operator of both equations (these are Laplacian and MOND operator)
  //
  // The source is not taken into account (not even the mass term).
  //
  // To be used for calculating residual and truncation error.
  //

  double nb_sum;
  double ev_mu[3][2];

  double sum;
  
  // Evaluate the Laplacian:
  //========================
  nb_sum = 
    phi[i  ][j  ][km1] + phi[i  ][j  ][kp1] +
    phi[i  ][jm1][k  ] + phi[i  ][jp1][k  ] +
    phi[im1][j  ][k  ] + phi[ip1][j  ][k  ];

  *operator_phi = (nb_sum - 6.0*phi[i][j][k])/pow2(h) + a2M2*(phi[i][j][k]+0.0*chi[i][j][k]);

  // Evaluate the MOND operator:
  //============================
  ev_mu[X][MINUS] = mu(dist( ( chi[i  ][j  ][k  ] - chi[im1][j  ][k  ] ) / h,
			     ( chi[i  ][jp1][k  ] + chi[im1][jp1][k  ] - chi[i  ][jm1][k  ] - chi[im1][jm1][k  ] ) / fourh,
			     ( chi[i  ][j  ][kp1] + chi[im1][j  ][kp1] - chi[i  ][j  ][km1] - chi[im1][j  ][km1] ) / fourh )/aa0 );
	      
  ev_mu[X][PLUS]  = mu(dist( ( chi[ip1][j  ][k  ] - chi[i  ][j  ][k  ] ) / h,
			     ( chi[ip1][jp1][k  ] + chi[i  ][jp1][k  ] - chi[ip1][jm1][k  ] - chi[i  ][jm1][k  ] ) / fourh,
			     ( chi[ip1][j  ][kp1] + chi[i  ][j  ][kp1] - chi[ip1][j  ][km1] - chi[i  ][j  ][km1] ) / fourh )/aa0 );
	      
  ev_mu[Y][MINUS] = mu(dist( ( chi[ip1][j  ][k  ] + chi[ip1][jm1][k  ] - chi[im1][j  ][k  ] - chi[im1][jm1][k  ] ) / fourh,
			     ( chi[i  ][j  ][k  ] - chi[i  ][jm1][k  ] ) / h,
			     ( chi[i  ][j  ][kp1] + chi[i  ][jm1][kp1] - chi[i  ][j  ][km1] - chi[i  ][jm1][km1] ) / fourh )/aa0 );
	      
  ev_mu[Y][PLUS]  = mu(dist( ( chi[ip1][jp1][k  ] + chi[ip1][j  ][k  ] - chi[im1][jp1][k  ] - chi[im1][j  ][k  ] ) / fourh,
			     ( chi[i  ][jp1][k  ] - chi[i  ][j  ][k  ] ) / h,
			     ( chi[i  ][jp1][kp1] + chi[i  ][j  ][kp1] - chi[i  ][jp1][km1] - chi[i  ][j  ][km1] ) / fourh )/aa0 );
	      
  ev_mu[Z][MINUS] = mu(dist( ( chi[ip1][j  ][k  ] + chi[ip1][j  ][km1] - chi[im1][j  ][k  ] - chi[im1][j  ][km1] ) / fourh,
			     ( chi[i  ][jp1][k  ] + chi[i  ][jp1][km1] - chi[i  ][jm1][k  ] - chi[i  ][jm1][km1] ) / fourh,
			     ( chi[i  ][j  ][k  ] - chi[i  ][j  ][km1] ) / h )/aa0 );
	      
  ev_mu[Z][PLUS]  = mu(dist( ( chi[ip1][j  ][kp1] + chi[ip1][j  ][k  ] - chi[im1][j  ][kp1] - chi[im1][j  ][k  ] ) / fourh,
			     ( chi[i  ][jp1][kp1] + chi[i  ][jp1][k  ] - chi[i  ][jm1][kp1] - chi[i  ][jm1][k  ] ) / fourh,
			     ( chi[i  ][j  ][kp1] - chi[i  ][j  ][k  ] ) / h )/aa0 );
	      
  sum = 
    ev_mu[X][MINUS] * (chi[im1][j  ][k  ] - chi[i  ][j  ][k  ]) +
    ev_mu[X][PLUS]  * (chi[ip1][j  ][k  ] - chi[i  ][j  ][k  ]) +
    ev_mu[Y][MINUS] * (chi[i  ][jm1][k  ] - chi[i  ][j  ][k  ]) +
    ev_mu[Y][PLUS]  * (chi[i  ][jp1][k  ] - chi[i  ][j  ][k  ]) +
    ev_mu[Z][MINUS] * (chi[i  ][j  ][km1] - chi[i  ][j  ][k  ]) +
    ev_mu[Z][PLUS]  * (chi[i  ][j  ][kp1] - chi[i  ][j  ][k  ]);

  *operator_chi = sum/pow2(h) + a2M2*(0.0*phi[i][j][k]+chi[i][j][k]);
	      
  
}


//====================================
// mu
//====================================
double mu(double x)
{
  
  // Standard mu function.

  return x/(1.0+x);

}
 

/** Evaluates differential operators on the grid given by the second argument.
 * Note that the grid does not have to necessarily be the grid gg->g.
 * we do not store the result in the correponding variables in gg (instead 
 * we store it in the arguments res_phi, res_chi).  This is in case you wat 
 * to store it somewhere else.  **/
void evaluate_differential_operator_mm(struct grids_gravity_mm *gg, int grid, double ***phi, double ***chi, double ***operator_phi, double ***operator_chi)
{

    int i, j, k;

  // Optimization:
  int ip1, im1, jp1, jm1, kp1, km1;
  int nm1;
  double aa0, a2M2;
  double h, fourh;

  aa0                = gg->aa0;
  a2M2               = gg->a2M2;
  nm1                = gg->n[grid]-1;
  h                  = gg->h[grid];
  fourh              = 4.0*h;

#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(gg, aa0, h, fourh, nm1, phi, chi, operator_phi, operator_chi)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<= nm1; i++)
      {
	ip1 = (i==nm1 ? 0   : (i+1));  // periodic stuff
	im1 = (i==0   ? nm1 : (i-1));
	for(j=0; j<= nm1; j++)
	  {
	    jp1 = (j==nm1 ? 0   : (j+1));  // periodic stuff
	    jm1 = (j==0   ? nm1 : (j-1));
	    for(k=0; k<= nm1; k++)
	      {
		kp1 = (k==nm1 ? 0   : (k+1));  // periodic stuff
		km1 = (k==0   ? nm1 : (k-1));
		
		// Evaluate differential operators:
		L_mm(phi, chi, 
		     aa0, a2M2, 
		     h, fourh, 
		     i, j, k, ip1, im1, jp1, jm1, kp1, km1, 
		     &operator_phi[i][j][k], &operator_chi[i][j][k]);
		
	      }
	  }
      }
#ifdef OPENMP
  }
#endif
  
}
