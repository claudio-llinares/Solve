#include<stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"
//#include "../../grid/get_mean_max_two_fields.h"

#include "solve_mond_mass_multi.h"
#include "operators_mond_mass.h"

#include "../commons_multigrid/restrict.h"
#include "../commons_multigrid/prolonge.h"   // For truncation error
//#include "../get_mond_operator.h"

#include "../../libs/utils.h"                // For distance

#ifdef OPENMP
#include <omp.h>
#endif


// two colors
int sx_first_mm[8] = {0, 1, 1, 0,  1, 0, 0, 1};
int sy_first_mm[8] = {0, 1, 0, 1,  0, 1, 0, 1};
int sz_first_mm[8] = {0, 0, 1, 1,  0, 0, 1, 1};

void gs_node_mm(double rho, double source_phi, double source_chi, double ***phi, double ***chi, double cte_poisson_over_a, double aa0, double h, double fourh, int i, int j, int k, int ip1, int im1, int jp1, int jm1, int kp1, int km1);



//======================================
// gs_mm
//======================================
void gs_mm(struct grids_gravity_mm *gg)
{
  // Do Gauss-Seidel in all the points.
    
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
  //double ***res_phi, ***res_chi;   // Optimization
  double cte_poisson_over_a, aa0, h, fourh;

  grid       = gg->g;
  rho        = gg->rho[grid];
  source_phi = gg->source_phi[grid];
  source_chi = gg->source_chi[grid];
  phi        = gg->phi[grid];
  chi        = gg->chi[grid];
  //res_phi    = gg->res_phi[grid];
  //res_chi    = gg->res_chi[grid];
  cte_poisson_over_a = gg->cte_poisson_over_a;
  aa0        = gg->aa0;
  nm1        = gg->n[grid]-1;
  h          = gg->h[grid];
  fourh      = 4.0*h;

  ggn = gg->n[gg->g];
  nm1 = ggn - 1;
  grid = gg->g;

  
  // Two colors:
  //============
  for(color=0; color<=7; color++)
    {
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(gg, sx_first_mm, sy_first_mm, sz_first_mm, color, ggn, nm1, grid, cte_poisson_over_a, aa0, h, fourh)
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
			       cte_poisson_over_a, aa0, h, fourh, 
			       i, j, k, ip1, im1, jp1, jm1, kp1, km1);
		  }
	      }
	  }
#ifdef OPENMP
      }
#endif
    }
   
}


//======================================
// gs_node_mm
//======================================
void gs_node_mm(double rho, double source_phi, double source_chi, double ***phi, double ***chi, double cte_poisson_over_a, double aa0, double h, double fourh, int i, int j, int k, int ip1, int im1, int jp1, int jm1, int kp1, int km1)
{
  //
  // Determination of the new potential on a given node for one
  // iteration step of the explicit solver.
  //

  double nb_sum;
  double ev_mu[3][2];
  double sum1, sum2;

  // Evolve phi:
  //============
  nb_sum = 
    phi[i  ][j  ][km1] + phi[i  ][j  ][kp1] +
    phi[i  ][jm1][k  ] + phi[i  ][jp1][k  ] +
    phi[im1][j  ][k  ] + phi[ip1][j  ][k  ];
  
  phi[i][j][k] = ( nb_sum - pow2(h)*(cte_poisson_over_a*rho+source_phi) ) / 6.0;

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
  
  sum1 = ev_mu[X][MINUS] * chi[im1][j  ][k  ] +
    ev_mu[X][PLUS]  * chi[ip1][j  ][k  ] +
    ev_mu[Y][MINUS] * chi[i  ][jm1][k  ] +
    ev_mu[Y][PLUS]  * chi[i  ][jp1][k  ] +
    ev_mu[Z][MINUS] * chi[i  ][j  ][km1] +
    ev_mu[Z][PLUS]  * chi[i  ][j  ][kp1];
  
  sum2 = ev_mu[X][MINUS] +
    ev_mu[X][PLUS] +
    ev_mu[Y][MINUS] +
    ev_mu[Y][PLUS] +
    ev_mu[Z][MINUS] +
    ev_mu[Z][PLUS];

  chi[i][j][k] = ( sum1 - pow2(h)*(cte_poisson_over_a*rho+source_chi) ) / sum2;
  
}

//======================================
// residual_mm
//======================================
void residual_mm(struct grids_gravity_mm *gg, int grid, double ***res_phi, double ***res_chi)
{

  // Calculates residual on the grid grid.
  // Note that the grid does not have to necessarily be the grid gg->g.

  double operator_phi, operator_chi;
  int i, j, k;

  // Optimization:
  int ip1, im1, jp1, jm1, kp1, km1;
  int nm1;
  double ***rho, ***source_phi, ***source_chi;
  double ***phi, ***chi;
  //double ***res_phi, ***res_chi;   // Optimization
  double cte_poisson_over_a, aa0, h, fourh;

  rho        = gg->rho[grid];
  source_phi = gg->source_phi[grid];
  source_chi = gg->source_chi[grid];
  phi = gg->phi[grid];
  chi = gg->chi[grid];
  //res_phi    = gg->res_phi[grid];
  //res_chi    = gg->res_chi[grid];
  cte_poisson_over_a = gg->cte_poisson_over_a;
  aa0        = gg->aa0;
  nm1        = gg->n[grid]-1;
  h          = gg->h[grid];
  fourh      = 4.0*h;


#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1, operator_phi, operator_chi) shared(gg, cte_poisson_over_a, aa0, h, fourh, nm1, source_phi, source_chi, res_phi, res_chi)
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
		     aa0, h, fourh, 
		     i, j, k, ip1, im1, jp1, jm1, kp1, km1, 
		     &operator_phi, &operator_chi);
		
		// Residual of equation for phi and chi:
		res_phi[i][j][k] = operator_phi - cte_poisson_over_a*rho[i][j][k] - source_phi[i][j][k];
		res_chi[i][j][k] = operator_chi - cte_poisson_over_a*rho[i][j][k] - source_chi[i][j][k];
	      }
	  }
      }
#ifdef OPENMP
  }
#endif
  
}


//======================================
//  trunc_err
//======================================
void trunc_err_mm(struct grids_gravity_mm *gg, double a)
{
  
  // Estimates truncation error.

  //double max_tau;
  //double tau_phi, tau_chi;
  /*int i, j, k;
  
  // optimization
  int nm1   = gg->n[gg->g]-1;
  int nm1m1 = gg->n[gg->g+1]-1;
  int ip1, im1, jp1, jm1, kp1, km1;

  double aa0   = a*(gg->par)->mm_a0;
  double fourh = 4.0*gg->h[gg->g];
  double fourh_g_plus_one = 4.0*gg->h[gg->g+1];
  // tau = P[ L(R(fi)) - R(L(fi)) ]
  double eval_L_phi, eval_L_chi;

  // L(fi):
  //=======
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(a, gg, nm1, fourh, aa0)
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
	      // This routine returns the operators in the last two arguments:
	      L_mm(i, j, k, ip1, im1, jp1, jm1, kp1, km1, gg->phi[gg->g], gg->chi[gg->g], gg->h[gg->g], fourh, a, aa0, &(gg->temp_phi[gg->g][i][j][k]), &(gg->temp_chi[gg->g][i][j][k]));
	    }
	}
    }
#ifdef OPENMP
  }
#endif

  // R(L(fi)):
  //==========
  restrict_grid(gg->temp_phi[gg->g], gg->temp_phi[gg->g+1], gg->n[gg->g]);
  restrict_grid(gg->temp_chi[gg->g], gg->temp_chi[gg->g+1], gg->n[gg->g]);
  // R(fi):
  //=======
  restrict_grid(gg->phi[gg->g], gg->temp2_phi[gg->g+1], gg->n[gg->g]);
  restrict_grid(gg->chi[gg->g], gg->temp2_chi[gg->g+1], gg->n[gg->g]);
  
  // Tau in coarse grid:
  //====================
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1, eval_L_phi, eval_L_chi) shared(nm1m1, gg, a, four_g_plus_one, aa0)
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

		// This routine returns the operators in the last two arguments:
		L_mm(gg->rho[grid][i][j][k], gg->source_phi[grid][i][j][k], gg->source_chi[grid][i][j][k], gg->phi[grid], gg->chi[grid], 
		     gg->cte_poisson_over_a, gg->h[grid], gg->fourh[grid], gg->aa0, 
		     i, j, k, ip1, im1, jp1, jm1, kp1, km1, 
		     &operator_phi, &operator_chi);
		//L_mm(i, j, k, ip1, im1, jp1, jm1, kp1, km1, gg->temp2_phi[gg->g+1], gg->temp2_chi[gg->g+1], gg->h[gg->g+1], fourh_g_plus_one, a, aa0, &eval_L_phi, &eval_L_chi);

		gg->trunc_phi[gg->g+1][i][j][k] = eval_L_phi - gg->temp_phi[gg->g+1][i][j][k];
		gg->trunc_chi[gg->g+1][i][j][k] = eval_L_chi - gg->temp_chi[gg->g+1][i][j][k];
	      }
	  }
      }
#ifdef OPENMP
  }
#endif
	  
  // Extrapolate tau to fine grid:
  //==============================
  prolonge(gg->trunc_phi[gg->g+1], gg->trunc_phi[gg->g], gg->n[gg->g+1], gg->h[gg->g+1]);
  prolonge(gg->trunc_chi[gg->g+1], gg->trunc_chi[gg->g], gg->n[gg->g+1], gg->h[gg->g+1]);

  // Calculate norm:
  //================
  tau_phi = 0.0;
  tau_chi = 0.0;
  max_tau = 0.0;
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
	      
	      tau_phi += pow2(gg->trunc_phi[gg->g][i][j][k]);
	      tau_chi += pow2(gg->trunc_chi[gg->g][i][j][k]);
	      
	      if(   fabs(gg->trunc_phi[gg->g][i][j][k])>max_tau || fabs(gg->trunc_chi[gg->g][i][j][k])>max_tau   )
#ifdef OPENMP
#pragma omp critical
#endif
		{
		  if(   fabs(gg->trunc_phi[gg->g][i][j][k])>max_tau || fabs(gg->trunc_chi[gg->g][i][j][k])>max_tau   )
		    max_tau = new_max(  fabs(gg->trunc_phi[gg->g][i][j][k]), fabs(gg->trunc_chi[gg->g][i][j][k])  );
		}
	    }
	}
    }
  
  gg->max_trunc[gg->g] = max_tau;

  return sqrt(tau_phi) + sqrt(tau_chi);
		       */  
}
