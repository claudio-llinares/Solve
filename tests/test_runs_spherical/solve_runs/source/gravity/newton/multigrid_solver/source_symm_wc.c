#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../../solve.h"
#include "../commons_multigrid/commons_multigrid.h"
#include "solve_symm_wc_multi.h"  // for X, Y, Z, MINUS, PLUS
#include "source_symm_wc.h"

//===================================
// diff_operator_symm_wc_optimized
//===================================
void diff_operator_symm_wc(struct params *par, double ***phi, double h, int i, int j, int k, int ip1, int im1, int jp1, int jm1, int kp1, int km1, long int g, double *operator, double *doperator_dw)
{

  // Calculates the differential operator (left hand side of original equation)
  // in a given node.
  //
  // Arguments:
  //          - g = grid in which we are sitting.
  // Returns:
  //          - operator.
  //          - derivative with respect to the field assumiong the coefficients 
  //            are constant.
  // Comments:
  //          - The first order formulas have one half that comes from
  //            the definition of f_{i+1/2).
  //
  
  double f[3][2];   // function f in between the nodes.
                    // f[X][MINUS] = f(\omega_{i-1/2})
  double sum1, sum2;
  
  double change;
  //double h = gg->h[g];
  
  change = par->symm_wc_change;

  // First order
  f[X][MINUS] = 0.5*(dchange_domega_symm_wc(change, phi[im1][j  ][k  ]) +
  		     dchange_domega_symm_wc(change, phi[i  ][j  ][k  ]));
  f[X][PLUS]  = 0.5*(dchange_domega_symm_wc(change, phi[i  ][j  ][k  ]) +
  		     dchange_domega_symm_wc(change, phi[ip1][j  ][k  ]));
  
  f[Y][MINUS] = 0.5*(dchange_domega_symm_wc(change, phi[i  ][jm1][k  ]) +
  		     dchange_domega_symm_wc(change, phi[i  ][j  ][k  ]));
  f[Y][PLUS]  = 0.5*(dchange_domega_symm_wc(change, phi[i  ][j  ][k  ]) +
  		     dchange_domega_symm_wc(change, phi[i  ][jp1][k  ]));
  
  f[Z][MINUS] = 0.5*(dchange_domega_symm_wc(change, phi[i  ][j  ][km1]) +
  		     dchange_domega_symm_wc(change, phi[i  ][j  ][k  ]));
  f[Z][PLUS]  = 0.5*(dchange_domega_symm_wc(change, phi[i  ][j  ][k  ]) +
  		     dchange_domega_symm_wc(change, phi[i  ][j  ][kp1]));
  
  sum1 =
    f[X][MINUS] * phi[im1][j  ][k  ] +
    f[X][PLUS]  * phi[ip1][j  ][k  ] +
    f[Y][MINUS] * phi[i  ][jm1][k  ] +
    f[Y][PLUS]  * phi[i  ][jp1][k  ] +
    f[Z][MINUS] * phi[i  ][j  ][km1] +
    f[Z][PLUS]  * phi[i  ][j  ][kp1];
  
  sum2 =
    f[X][MINUS] + f[X][PLUS] +
    f[Y][MINUS] + f[Y][PLUS] +
    f[Z][MINUS] + f[Z][PLUS];
  
  *operator = ( sum1 - sum2*phi[i][j][k] ) / pow2(h);
  
  *doperator_dw = -sum2 / pow2(h);

}
//==========================
// l_symm_wc
//==========================
double l_symm_wc(struct params *par, double operator, double phi, double rho, double source, double a)
{

  double val;
   
  fprintf(stderr, "l_symm_wc: check what is rho in this routine!!!  Aboerting.\n");
  exit(0);

  if(par->solver==3)
    {
      // General change of variables (symm_wceleon):
      val = operator - 
	par->symm_wc_cte_l * ( rho - 
			    pow(change_symm_wc(par->symm_wc_change, phi), 
				-(1.0+par->symm_wc_power))
			    ) - 
	source;
    }
  else
    {
      // General change of variables (symmetron):
      val = operator - 
	0.5*pow2(a/par->symm_range)*
	(-change_symm_wc(par->symm_wc_change, phi) + pow3(change_symm_wc(par->symm_wc_change, phi)) +
	 pow3(par->symm_assb/a)*rho*change_symm_wc(par->symm_wc_change, phi)) - source;
    }
  
  return val;

}

//==========================
// l_prime_symm_wc
//==========================
double l_prime_symm_wc(struct params *par, double doperator_dw, double phi, double rho, double source, double a)
{

  double val;

  // Symm_Wcaleon:
  //===========
  // Optimized version of a power law change of variables:
  // NOT WORKING
  //val = doperator_dw -  
  //  par->symm_wc_cte_l * ( par->symm_wc_change*(1.0+par->symm_wc_power) * 
  //			pow(phi, -(1.0+par->symm_wc_power)*(par->symm_wc_change-2)-1) );

  if(par->solver==3)
    {
      // General change of variables (symm_wceleon):
      val = doperator_dw -  
	(1.0/pow2(par->symm_wc_range)) * ( pow(change_symm_wc(par->symm_wc_change, phi), 
					    -(par->symm_wc_power+2)) * 
					dchange_domega_symm_wc(par->symm_wc_change, phi));
    }
  else
    {
      // General change of variables (symmetron):
      val = doperator_dw -  
	0.5*pow2(a/par->symm_range) *
	(-1.0 + 
	 3.0*pow2(change_symm_wc(par->symm_wc_change, phi)) +      
	 pow3(par->symm_assb/a)*rho) * dchange_domega_symm_wc(par->symm_wc_change, phi);
    }
  
  // Some test:
  //===========
  //if(val==0)
  //  {
  //    fprintf(stderr, "l_prime_symm_wc:  the value is 0.  Aborting...");
  //    exit(0);
  //  }

  return val;

}

//==========================
// dchange_domega_symm_wc
//==========================
double dchange_domega_symm_wc(double change, double omega)
{
  
  // Free function in the operator:
  //     DIV( f(w) GRAD(w) )
  //
  // This comes fro mthe fact that:
  //     chi = g(omega)
  //     f = dg\d omega
  // implies
  //     LAPLACIAN chi = DIV( f(w) GRAD(w) )
  //
  
  //return log(10.0)*pow(10,omega);
  //return exp(omega);
  return (change) * pow(omega, (change)-1);

}

//==========================
// change_symm_wc
//==========================
double change_symm_wc(double change, double omega)
{
  
  // Change of variable that put the symm_wcaleon operator as:
  //
  //     DIV( f(w) GRAD(w) )
  //
  // This comes from the fact that:
  //     chi = g(omega)
  //     f = dg/d omega
  // implies
  //     LAPLACIAN chi = DIV( f(w) GRAD(w) )
  //
  // This function returns g.
  // This routine is needed for the source of the equation.
  // (the power law change is optimited in the routines l_symm_wc and l_prime_symm_wc,
  // so check them if you want to make changes here).
  // This is also used to recover the original potential chi.
  //
  
  //return pow(10,omega);
  //return exp(omega);
  return pow(omega, change);

}

//==========================
// change_symm_wc
//==========================
double inverse_change_symm_wc(double change, double omega)
{
  
  // Change of variable that put the symm_wcaleon operator as:
  //
  //     DIV( f(w) GRAD(w) )
  //
  // This comes from the fact that:
  //     chi = g(omega)
  //     f = dg/d omega
  // implies
  //     LAPLACIAN chi = DIV( f(w) GRAD(w) )
  //
  // This function returns inverse of g.
  // This routine is not needed to solve the equation, but only to get
  // the initial guess.
  //
  
  //return pow(10,omega);
  //return log(omega);
  return pow(omega, 1.0/change);

}
