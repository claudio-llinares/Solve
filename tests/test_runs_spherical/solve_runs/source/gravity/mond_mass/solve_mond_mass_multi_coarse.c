#include<stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../../solve.h"

#include "solve_mond_mass_multi.h"
#include "source_mond_mass.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_poly.h>  // for cubic equation for explicit method

// two colors
int sx_first_coarse_mm[8] = {0, 1, 1, 0,  1, 0, 0, 1};
int sy_first_coarse_mm[8] = {0, 1, 0, 1,  0, 1, 0, 1};
int sz_first_coarse_mm[8] = {0, 0, 1, 1,  0, 0, 1, 1};

// lexicographic
//int sx_first_coarse[8] = {0, 1, 0, 1, 0, 1, 0, 1};
//int sy_first_coarse[8] = {0, 0, 1, 1, 0, 0, 1, 1};
//int sz_first_coarse[8] = {0, 0, 0, 0, 1, 1, 1, 1};

double residual_node_mm_coarse(struct params *par, struct grids_gravity_mm *gg, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int g);


//======================================
// residual_mm_coarse
//======================================
double residual_mm_coarse(struct params *par, struct grids_gravity_mm *gg, double a)
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
#pragma omp parallel reduction(+:res) private(i, j, k, ip1, im1, jp1, jm1, kp1, km1, r) shared(nm1, gg, a)
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
		r = residual_node_mm_coarse(par, gg, i, j, k, a, ip1, im1, jp1, jm1, kp1, km1, gg->g);
		gg->res[gg->g][i][j][k] = r;

		// Calculate max residual
		if(fabs(r) > fabs(gg->max_resid[gg->g]))
#ifdef OPENMP
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
// gS_node_newton
//======================================
double residual_node_mm_coarse(struct params *par, struct grids_gravity_mm *gg, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int g)
{
  //
  // Calculate residual defined by l_mm_operator
  // (second or fourth order)
  //

  double l;

  // Fourth order:
  //==============
  l = (*(gg->l_mm_p))(par, gg, g, gg->chi_bar_mm[gg->g], gg->source[g][i][j][k], i, j, k, a, ip1, im1, jp1, jm1, kp1, km1);
  
  return l;

}


