#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

#include <fftw3.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"


#include "../commons_multigrid/commons_multigrid.h"  // for fourth order Laplacian
#include "source_fifth.h"

//=======================
// init_phi_fifth
//=======================
void init_phi_fifth(struct params *par, double ***phi, double ***rho)
{

  int i, j, k;

  fprintf(stderr, "init_phi_fifth: phib = %g\n", phi_b(par, par->a));
 
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
	{ 
	  
	  if(par->symm_assb>=par->a)
	    phi[i][j][k] = 1.0e-13;
	  else
	    phi[i][j][k] = 0.5;

	}
#ifdef OPENMP
  }
#endif

}
