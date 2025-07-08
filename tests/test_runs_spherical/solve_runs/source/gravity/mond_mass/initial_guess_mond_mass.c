#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"

#include "solve_mond_mass_multi.h"

#include "../../libs/ran3.h"


//=======================
// initial_guess_mm
//=======================
void initial_guess_mm(struct params *par, double ***chi_bar, double ***phi_bar)
{

  // Initial guess for multigrid solver.
  // Use this if you do not have a better solution
  // (such as the solution is previous time step).

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
	  phi_bar[i][j][k] = 0.0;
	  chi_bar[i][j][k] = ran3(&iseed);
	}
  
#ifdef OPENMP
  }
#endif

}
