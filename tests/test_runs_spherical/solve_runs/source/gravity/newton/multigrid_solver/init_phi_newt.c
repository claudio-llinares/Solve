#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

#include <fftw3.h>

#include "../../../solve.h"
#include "../../../grid/init_grid.h"


#include "../../commons_multigrid/commons_multigrid.h"  // for fourth order Laplacian
#include "source_newt.h"

//=======================
// init_phi_newt
//=======================
void init_phi_newt(struct params *par, double ***phi, double ***rho)
{

  int i, j, k;

  // Initialize potential to random values:
  //=======================================
#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(par, phi, rho)
  {
#pragma omp for schedule(static)
#endif
  for(i=0; i<par->grid; i++)  // x directon
    for(j=0; j<par->grid; j++)  // y direction
      for(k=0; k<par->grid; k++)  // z direction
	phi[i][j][k] = 0.0;
#ifdef OPENMP
  }
#endif

}
