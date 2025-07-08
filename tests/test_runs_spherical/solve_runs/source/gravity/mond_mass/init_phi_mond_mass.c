#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

#include <fftw3.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"

#include "../commons_multigrid/commons_multigrid.h"  // for fourth order Laplacian
#include "source_mond_mass.h"

#include "../../libs/ran3.h"


//=======================
// init_phi_mond_mass
//=======================
void init_phi_mm(struct params *par, double ***chi_bar, double ***phi_bar)
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
	  chi_bar[i][j][k] = ran3(&iseed);
	  phi_bar[i][j][k] = 0.0;
	}
  
#ifdef OPENMP
  }
#endif

}
