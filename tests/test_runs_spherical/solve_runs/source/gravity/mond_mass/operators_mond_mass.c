#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"
#include "../../libs/utils.h"

#include "solve_mond_mass_multi.h"


double mu(double x);


//======================================
// L_mm
//======================================
//void L_mm(double rho, double ***phi, double ***chi, double cte_poisson_over_a, double aa0, double h, double fourh, int i, int j, int k, int ip1, int im1, int jp1, int jm1, int kp1, int km1, double *operator_phi, double *operator_chi)
void L_mm(double ***phi, double ***chi, double aa0, double h, double fourh, int i, int j, int k, int ip1, int im1, int jp1, int jm1, int kp1, int km1, double *operator_phi, double *operator_chi)
{

  //
  // Differential operator including the source of the original equation.
  //
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

  *operator_phi = (nb_sum - 6.0*phi[i][j][k])/pow2(h); // - cte_poisson_over_a*rho;

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

  *operator_chi = sum/pow2(h);   // - cte_poisson_over_a*rho;
	      
  
}


//====================================
// mu
//====================================
double mu(double x)
{
  
  // There is an additional 1 because this is 1+mu

  return x/(1.0+x);

}
 
