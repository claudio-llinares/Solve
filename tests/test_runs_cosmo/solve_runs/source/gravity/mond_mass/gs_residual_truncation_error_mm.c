/**
 * Routines for doing one Gauss-Seidel sweep, calculating residual and 
 * truncation error (all in one grid).
 * 
 *  Note:
 *  This is basically the core of the solver where the equations are included.
 *  See also ``differential_operators.c``.
**/
#include<stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"

#include "solve_multigrid_mm.h"
#include "differential_operators_mm.h"

#include "../commons_multigrid/restrict.h"
#include "../commons_multigrid/prolonge.h"   // For truncation error

#include "../../libs/utils.h"                // For distance

#ifdef OPENMP
#include <omp.h>
#endif


// two colors
int sx_first_mm[8] = {0, 1, 1, 0,  1, 0, 0, 1};
int sy_first_mm[8] = {0, 1, 0, 1,  0, 1, 0, 1};
int sz_first_mm[8] = {0, 0, 1, 1,  0, 0, 1, 1};

void gs_node_mm(double rho, double source_phi, double source_chi, double ***phi, double ***chi, double cte_poisson_over_a, double aa0, double a2M2, double h, double fourh, int i, int j, int k, int ip1, int im1, int jp1, int jm1, int kp1, int km1);



//======================================
// gs_mm
//======================================
void gs_mm(struct grids_gravity_mm *gg)
{
  // Do Gauss-Seidel on the grid gg->g.
    
  int i, j, k;
  
  int color;

  // optimization
  int ggn;  
  int ip1, im1, jp1, jm1, kp1, km1;
  int nm1;
  int grid;

  // Optimization:
  double ***rho, ***source_phi, ***source_chi;
  double ***phi, ***chi;
  double cte_poisson_over_a, aa0, a2M2;
  double h, fourh;

  grid               = gg->g;
  rho                = gg->rho[grid];
  source_phi         = gg->source_phi[grid];
  source_chi         = gg->source_chi[grid];
  phi                = gg->phi[grid];
  chi                = gg->chi[grid];
  cte_poisson_over_a = gg->cte_poisson_over_a;
  aa0                = gg->aa0;
  a2M2               = gg->a2M2;
  ggn                = gg->n[gg->g];
  nm1                = gg->n[grid]-1;
  h                  = gg->h[grid];
  fourh              = 4.0*h;
  
  // Two colors:
  //============
  for(color=0; color<=7; color++)
    {
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(gg, sx_first_mm, sy_first_mm, sz_first_mm, color, ggn, nm1, cte_poisson_over_a, aa0, a2M2, h, fourh, rho, source_phi, source_chi, phi, chi)
      {
#pragma omp for schedule(static)
#endif
	for(i=sx_first_mm[color]; i< ggn; i+=2)
	  {
	    ip1 = (i==nm1 ? 0   : (i+1));  // periodic stuff
	    im1 = (i==0   ? nm1 : (i-1));
	    for(j=sy_first_mm[color]; j< ggn; j+=2)
	      {
		jp1 = (j==nm1 ? 0   : (j+1));  // periodic stuff		
		jm1 = (j==0   ? nm1 : (j-1));
		for(k=sz_first_mm[color]; k< ggn; k+=2)
		  {
		    kp1 = (k==nm1 ? 0   : (k+1));  // periodic stuff
		    km1 = (k==0   ? nm1 : (k-1));
		    
		    // There are a lot of arguments here, but hopefully they will help in gaining speed...
		    gs_node_mm(rho[i][j][k], source_phi[i][j][k], source_chi[i][j][k], 
			       phi, chi, 
			       cte_poisson_over_a, aa0, a2M2, 
			       h, fourh, 
			       i, j, k, ip1, im1, jp1, jm1, kp1, km1);
		  }
	      }
	  }
#ifdef OPENMP
      }
#endif
    }
   
}


/** 
 * This routine does Gauss-Seidel (GS) iteration in one node.
 * 
 * Comments:
 *      - There are a lot of arguments to prevent the routine from accessing the pointer to ``grids_gravity_mm``
 *        and make it more efficient.
 *      - Whatever you put here is crucial for the performance of the code.  In particular, "if" statements can kill the performance.
 **/
void gs_node_mm(double rho, double source_phi, double source_chi, double ***phi, double ***chi, double cte_poisson_over_a, double aa0, double a2M2, double h, double fourh, int i, int j, int k, int ip1, int im1, int jp1, int jm1, int kp1, int km1)
{
  //
  // Determination of the new potential on a given node for one
  // iteration step of the explicit solver.
  //

  double nb_sum;
  double ev_mu[3][2];
  double sum1, sum2;

  double phi_old = phi[i][j][k];   // this is just to see if we can improve convergence.

  // Evolve phi:
  //============
  nb_sum = 
    phi[i  ][j  ][km1] + phi[i  ][j  ][kp1] +
    phi[i  ][jm1][k  ] + phi[i  ][jp1][k  ] +
    phi[im1][j  ][k  ] + phi[ip1][j  ][k  ];
  
  phi[i][j][k] = ( nb_sum - pow2(h)*(-   0.0   *   a2M2*chi[i][j][k] + cte_poisson_over_a*rho + source_phi) ) / (6.0-a2M2*pow2(h));

  // Evolve chi:
  //============
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
  
  sum1 = 
    ev_mu[X][MINUS] * chi[im1][j  ][k  ]    +    ev_mu[X][PLUS] * chi[ip1][j  ][k  ] +
    ev_mu[Y][MINUS] * chi[i  ][jm1][k  ]    +    ev_mu[Y][PLUS] * chi[i  ][jp1][k  ] +
    ev_mu[Z][MINUS] * chi[i  ][j  ][km1]    +    ev_mu[Z][PLUS] * chi[i  ][j  ][kp1];
  
  sum2 = 
    ev_mu[X][MINUS]    +    ev_mu[X][PLUS] +
    ev_mu[Y][MINUS]    +    ev_mu[Y][PLUS] +
    ev_mu[Z][MINUS]    +    ev_mu[Z][PLUS] -
    a2M2*pow2(h);

  chi[i][j][k] = ( sum1 - pow2(h)*(-   0.0    *     a2M2*phi_old + cte_poisson_over_a*rho + source_chi) ) / sum2;

  
}

/**
 * Calculates residual on the grid "grid".
 * Comments:
 *     - Note that the grid does not have to necessarily be the grid gg->g.
 *     - We do not store the result by default in the correponding variables in gg (instead we store it in 
 *       the arguments res_phi, res_chi).  This is in case you wat to store it somewhere else.
 **/
void residual_mm(struct grids_gravity_mm *gg, int grid, double ***res_phi, double ***res_chi)
{

   double operator_phi, operator_chi;
  int i, j, k;

  // Optimization:
  int ip1, im1, jp1, jm1, kp1, km1;
  int nm1;
  double ***rho, ***source_phi, ***source_chi;
  double ***phi, ***chi;
  double cte_poisson_over_a, aa0, a2M2;
  double h, fourh;

  rho                = gg->rho[grid];
  source_phi         = gg->source_phi[grid];
  source_chi         = gg->source_chi[grid];
  phi                = gg->phi[grid];
  chi                = gg->chi[grid];
  cte_poisson_over_a = gg->cte_poisson_over_a;
  aa0                = gg->aa0;
  a2M2               = gg->a2M2;
  nm1                = gg->n[grid]-1;
  h                  = gg->h[grid];
  fourh              = 4.0*h;


#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1, operator_phi, operator_chi) shared(gg, cte_poisson_over_a, aa0, a2M2, h, fourh, nm1, rho, source_phi, source_chi, res_phi, res_chi)
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
		     &operator_phi, &operator_chi);
		
		// Residual of equation for phi and chi:
		res_phi[i][j][k] = operator_phi - ( cte_poisson_over_a*rho[i][j][k] + source_phi[i][j][k]);
		res_chi[i][j][k] = operator_chi - ( cte_poisson_over_a*rho[i][j][k] + source_chi[i][j][k]);
	      }
	  }
      }
#ifdef OPENMP
  }
#endif
  
}


/**
 * Estimates truncation error (i.e. error associated to the discretizatio of the differential equation).
 * 
 * Comments:
 *     - Note that the grid does not have to necessarily be the grid gg->g.
 *     - The estimator is: tau = P[ D(R(phi)) - R(D(phi)) ]
 **/
void trunc_err_mm(struct grids_gravity_mm *gg, int grid, double ***trunc_phi, double ***trunc_chi)
{
  
  int i, j, k;
  
  // optimization
  int nm1   = gg->n[grid]-1;
  int nm1m1 = gg->n[grid+1]-1;
  int ip1, im1, jp1, jm1, kp1, km1;

  double operator_phi, operator_chi;


  // Optimization:
  double ***phi, ***chi;
  double aa0, a2M2;
  double h, fourh;

  phi                = gg->phi[grid];
  chi                = gg->chi[grid];
  aa0                = gg->aa0;
  a2M2               = gg->a2M2;
  nm1                = gg->n[grid]-1;
  h                  = gg->h[grid];
  fourh              = 4.0*h;

  // Evaluate differential operator in this grid:
  //=============================================
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(gg, nm1, h, fourh, aa0, a2M2, phi, chi)
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
		   &gg->temp_phi[grid][i][j][k], &gg->temp_chi[grid][i][j][k]);
	    }
	}
    }
#ifdef OPENMP
  }
#endif

  // Restrict evaluation of differential operator:
  //==============================================
  restrict_grid(gg->temp_phi[grid], gg->temp_phi[grid+1], gg->n[grid]);
  restrict_grid(gg->temp_chi[grid], gg->temp_chi[grid+1], gg->n[grid]);

  // Restrict solution:
  //===================
  restrict_grid(gg->phi[grid], gg->temp2_phi[grid+1], gg->n[grid]);
  restrict_grid(gg->chi[grid], gg->temp2_chi[grid+1], gg->n[grid]);
  
  // Put everything together to get tau:
  //====================================
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1, operator_phi, operator_chi) shared(nm1m1, gg, aa0, grid)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<= nm1m1; i++)
      {
	ip1 = (i==nm1m1 ? 0     : (i+1));  // periodic stuff
	im1 = (i==0     ? nm1m1 : (i-1));
	for(j=0; j<= nm1m1; j++)
	  {
	    jp1 = (j==nm1m1 ? 0     : (j+1));  // periodic stuff
	    jm1 = (j==0     ? nm1m1 : (j-1));
	    for(k=0; k<= nm1m1; k++)
	      {
		kp1 = (k==nm1m1 ? 0     : (k+1));  // periodic stuff
		km1 = (k==0     ? nm1m1 : (k-1));

		// Evaluate operator on the restricted solution:
		// This routine returns the operators in the last two arguments:
		L_mm(gg->temp2_phi[grid+1], gg->temp2_chi[grid+1],  
		     aa0, a2M2, 
		     gg->h[grid+1], gg->fourh[grid+1], 
		     i, j, k, ip1, im1, jp1, jm1, kp1, km1, 
		     &operator_phi, &operator_chi);
		
		gg->trunc_phi[grid+1][i][j][k] = operator_phi - gg->temp_phi[grid+1][i][j][k];
		gg->trunc_chi[grid+1][i][j][k] = operator_chi - gg->temp_chi[grid+1][i][j][k];
	      }
	  }
      }
#ifdef OPENMP
  }
#endif
	  
  // Prolonge tau to fine grid:
  //===========================
  prolonge(gg->trunc_phi[grid+1], trunc_phi, gg->n[grid+1], gg->h[grid+1]);
  prolonge(gg->trunc_chi[grid+1], trunc_chi, gg->n[grid+1], gg->h[grid+1]);

        
}

