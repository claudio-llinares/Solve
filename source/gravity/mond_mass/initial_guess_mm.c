/**
 * Initial guess for multigrid solver.
 * 
 * Comments:
 *  - Use this if you do not have a better solution (such as the solution 
 *    is previous time step).
 *  - This must be called before you call the multigrid solver.
 *  - For the MOND case, zero is not good enough (that's why we put a 
 *    random number).
 *  - A better initial guesses can be otained by what we have in ``solve_spherical_mond``.
 * 
 **/
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"

#include "solve_multigrid_mm.h"

#include "../../libs/ran3.h"


//=======================
// initial_guess_mm
//=======================
void initial_guess_mm(struct params *par, double ***chi_bar, double ***phi_bar)
{

  int i, j, k;

  long int iseed = -123456;

  // Initialize potential to random values:
  //=======================================
#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(par, chi_bar, phi_bar)
  {
#pragma omp for schedule(static)
#endif
  for(i=0; i<par->grid; i++)  // x directon
    for(j=0; j<par->grid; j++)  // y direction
      for(k=0; k<par->grid; k++)  // z direction
	{
	  //phi_bar[i][j][k] = 0.0;
#pragma omp critical               //  not including this line end up giving zeros in some places of chi_bar.
	  // A small perturbation around zero.
	  // Using simply ran3 shifts the boundary condition in infinity for the Newtonian limit.
	  chi_bar[i][j][k] = 0.01*(ran3(&iseed)-0.5);
	  phi_bar[i][j][k] = chi_bar[i][j][k];
	}
  
#ifdef OPENMP
  }
#endif

}
