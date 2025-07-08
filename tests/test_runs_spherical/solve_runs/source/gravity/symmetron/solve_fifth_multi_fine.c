#include<stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../solve.h"

#include "../commons_multigrid/commons_multigrid.h"
#include "../commons_multigrid/restrict.h"
#include "solve_fifth_multi.h"
#include "source_fifth.h"
#include "../commons_multigrid/prolonge.h"   // for truncation error

#ifdef OPENMP
#include <omp.h>
#endif
#include <gsl/gsl_poly.h>  // for cubic equation for explicit method

// two colors
int sx_first_fine[8] = {0, 1, 1, 0,  1, 0, 0, 1};
int sy_first_fine[8] = {0, 1, 0, 1,  0, 1, 0, 1};
int sz_first_fine[8] = {0, 0, 1, 1,  0, 0, 1, 1};

// lexicographic
//int sx_first_fine[8] = {0, 1, 0, 1, 0, 1, 0, 1};
//int sy_first_fine[8] = {0, 0, 1, 1, 0, 0, 1, 1};
//int sz_first_fine[8] = {0, 0, 0, 0, 1, 1, 1, 1};

void gs_node_fifth_fine(struct params *par, struct grids_gravity *gg, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1);
double residual_node_fifth_fine(struct params *par, struct grids_gravity *gg, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int g);

//======================================
// gs_fine_fifth
//======================================
void gs_fine_fifth(struct params *par, struct grids_gravity *gg, double a)
{
  // Actualices all interiors points.
    
  int i, j, k;
  
  int color;

  // optimization
  int ggn;  
  int ip1, im1, jp1, jm1, kp1, km1;
  int nm1;

  ggn = gg->n[gg->g];
  nm1 = ggn - 1;

  //FILE *fp;          // for testing colors.
  //char file[500];
  
  //fprintf(stderr, "Iteration on grid %ld (%ld)\n", gg->g, gg->n[gg->g]);
  // Two colors:
  //============
  //for(n=0; n<=2; n++)
    {
      for(color=0; color<=7; color++)
	{
	  //sprintf(file, "color.%d",color);
	  //fp = fopen(file, "w");
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(a, gg, sx_first, sy_first, sz_first,color, ggn, nm1)
	  {
#pragma omp for schedule(static)
#endif
	  for(i=sx_first[color]; i< ggn; i+=2)
	    {
	      ip1 = (i==nm1 ? 0   : (i+1));  // periodic stuff
	      im1 = (i==0   ? nm1 : (i-1));
	      for(j=sy_first[color]; j< ggn; j+=2)
		{
		  jp1 = (j==nm1 ? 0   : (j+1));  // periodic stuff
		  jm1 = (j==0   ? nm1 : (j-1));
		  for(k=sz_first[color]; k< ggn; k+=2)
		    {
		      kp1 = (k==nm1 ? 0   : (k+1));  // periodic stuff
		      km1 = (k==0   ? nm1 : (k-1));

		      gs_node_fifth_fine(par, gg, i, j, k, a, ip1, im1, jp1, jm1, kp1, km1);
		    }
		}
	    }
#ifdef OPENMP
	  }
#endif
	}
    }
   
}

//======================================
// gS_node_fifth_fine
//======================================
void gs_node_fifth_fine(struct params *par, struct grids_gravity *gg, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1)
{
  //
  // Determination of the new potential on a given node for one
  // iteration step of the explicit solver.
  //

  // optimization
  double h = gg->h[gg->g];
  
  double l, l_prime;
    
  l = (*(gg->l_fifth_operator))(par, gg, gg->g, 0.0, i, j, k, a, ip1, im1, jp1, jm1, kp1, km1);
  
  l_prime = (*(gg->l_prime_fifth_operator))(par, gg->fi[gg->g][i][j][k],
					    gg->rho[gg->g][i][j][k],
					    0.0,
					    h, a);

  gg->fi[gg->g][i][j][k] = gg->fi[gg->g][i][j][k] - l/l_prime;

  // For solution to be positive:
  //=============================
  if(1)
    {
      if(gg->fi[gg->g][i][j][k] < 0.0)
	gg->fi[gg->g][i][j][k] = -gg->fi[gg->g][i][j][k];
    }

  return;
  
}

//======================================
// residual_fifth_fine
//======================================
double residual_fifth_fine(struct params *par, struct grids_gravity *gg, double a)
{

  // Calculates residual.

  double res;
  double r;
  int i, j, k;

  // optimization
  int nm1 = gg->n[gg->g]-1;
  int ip1, im1, jp1, jm1, kp1, km1;

  res = 0.0;
  gg->max_resid[gg->g] = 0.0;
#ifdef OPENMP
#pragma omp parallel reduction(+:res) private(i, j, k, ip1, im1, jp1, jm1, kp1, km1, r) shared(nm1, gg, a, par)
  {
    res = 0.0;
    gg->max_resid[gg->g] = 0.0;
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
		
		// Calculate residual in node
		r = residual_node_fifth_fine(par, gg, i, j, k, a, ip1, im1, jp1, jm1, kp1, km1, gg->g);
		gg->res[gg->g][i][j][k] = r;

		// Calculate max residual
		if(fabs(r) > fabs(gg->max_resid[gg->g]))
#ifdef OPENMO
#pragma omp critical
#endif
		  {
		    if(fabs(r) > fabs(gg->max_resid[gg->g]))
		      gg->max_resid[gg->g] = fabs(r);
		  }
		
		// Calculate residual as norm2
		res += pow2(r);
	      }
	  }
      }
#ifdef OPENMP
  }
#endif
  
  return sqrt(res);

}


//======================================
// residual_node_fifth_fine
//======================================
double residual_node_fifth_fine(struct params *par, struct grids_gravity *gg, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int g)
{
  //
  // Calculate residual defined by l_fifth_operator
  // (second or fourth order)
  //

  double l;

  l = (*(gg->l_fifth_operator))(par, gg, g, 0.0, i, j, k, a, ip1, im1, jp1, jm1, kp1, km1);
  
  return l;

}

//======================================
// go_up_fine_fifth
//======================================
void go_up_fine_fifth(struct params *par, struct grids_gravity *gg, double a)
{
  //
  // Moves everything to the next coarser grid.
  //
    
  int i, j, k;

  // optimization:
  int nm1_p1 = gg->n[gg->g+1]-1;
  int ip1, im1, jp1, jm1, kp1, km1;

  // Restrict residual, density and solution:
  //=========================================
  restrict_grid(gg->res[gg->g], gg->res[gg->g+1], gg->n[gg->g]);
  restrict_grid(gg->rho[gg->g], gg->rho[gg->g+1], gg->n[gg->g]);
  restrict_grid(gg->fi[gg->g], gg->fi[gg->g+1], gg->n[gg->g]);
  //write_pot(-10, gg->g, a);
  //write_pot(-10, gg->g+1, a);

  // Get new source:
  //================
  // This will be restored after the loop
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(gg, a, nm1_p1)
  {
#pragma omp for schedule(static)
#endif
  for(i=0; i<= nm1_p1; i++)   // Fills all the values in the grid.
    {
      ip1 = (i==nm1_p1 ? 0    : (i+1));  // periodic stuff
      im1 = (i==0    ? nm1_p1 : (i-1));
      for(j=0; j<= nm1_p1; j++)
	{
	  jp1 = (j==nm1_p1 ? 0    : (j+1));  // periodic stuff
	  jm1 = (j==0    ? nm1_p1 : (j-1));
	  for(k=0; k<= nm1_p1; k++)
	    {
	      kp1 = (k==nm1_p1 ? 0    : (k+1));  // periodic stuff
	      km1 = (k==0    ? nm1_p1 : (k-1));

	      // Fourth order:
	      //==============
	      gg->source[gg->g+1][i][j][k] = -gg->res[gg->g+1][i][j][k] + 
		(*(gg->l_fifth_operator))(par, gg, gg->g+1, 0.0, i, j, k, a, ip1, im1, jp1, jm1, kp1, km1);
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
double trunc_err_fine_fifth(struct grids_gravity *gg, double a)
{
  
  // Estimates truncation error.

  double tau, max_tau;
  int i, j, k;
  
  // optimization
  int nm1   = gg->n[gg->g]-1;
  int nm1m1 = gg->n[gg->g+1]-1;
  int ip1, im1, jp1, jm1, kp1, km1;

  // tau = P[ L(R(fi)) - R(L(fi)) ]

  // L(fi):
  //=======
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(a, gg, nm1)
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
	      gg->temp[gg->g][i][j][k] =
		(*(gg->L_operator))(i, j, k, gg->fi[gg->g], gg->h[gg->g], a, ip1, im1, jp1, jm1, kp1, km1, gg->n[gg->g]);
	    }
	}
    }
#ifdef OPENMP
  }
#endif

  // R(L(fi)):
  //==========
  restrict_grid(gg->temp[gg->g], gg->temp[gg->g+1], gg->n[gg->g]);

  // R(fi):
  //=======
  restrict_grid(gg->fi[gg->g], gg->temp2[gg->g+1], gg->n[gg->g]);
  
  // Tau in coarse grid:
  //====================
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1) shared(nm1m1, gg, a)
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
		
		gg->trunc[gg->g+1][i][j][k] = 
		  (*(gg->L_operator))(i, j, k, gg->temp2[gg->g+1], gg->h[gg->g+1], a, ip1, im1, jp1, jm1, kp1, km1, gg->n[gg->g+1]) - 
		  gg->temp[gg->g+1][i][j][k];
	      }
	  }
      }
#ifdef OPENMP
  }
#endif
	  
  // Extrapolate tau to fine grid:
  //==============================
  prolonge(gg->trunc[gg->g+1], gg->trunc[gg->g], gg->n[gg->g+1], gg->h[gg->g+1]);
  // Calculate norm:
  //================
  tau = 0.0;
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
	      
	      tau += pow2(gg->trunc[gg->g][i][j][k]);
	      
	      if(fabs(gg->trunc[gg->g][i][j][k])>max_tau)
#ifdef OPENMP
#pragma omp critical
#endif
		{
		  if(fabs(gg->trunc[gg->g][i][j][k])>max_tau)
		    max_tau = fabs(gg->trunc[gg->g][i][j][k]);
		}
	    }
	}
    }
  
  gg->max_trunc[gg->g] = max_tau;

  return sqrt(tau);
    
}
