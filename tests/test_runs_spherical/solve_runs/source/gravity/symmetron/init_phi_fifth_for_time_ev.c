#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

#include <fftw3.h>

#include "../../solve.h"
#include "../../grid/init_grid.h"

#include "../commons_multigrid/commons_multigrid.h"  // for fourth order Laplacian
#include "source_fifth.h"

#include "../../libs/ran3.h"

//==============================
// init_phi_fifth_for_time_ev
//==============================
void init_phi_fifth_for_time_ev(struct params *par, struct grids *grid)
{

  int i, j, k;

  long seed;
  
  seed = -1;

  fprintf(stderr, "init_phi_fifth_for_time_ev: phib = %g\n", phi_b(par, 1.0));
 
  // Initialize potential to random values:
  //=======================================
  // Not possible to paraleize thans to ran3
  //#ifdef OPENMP
  //#pragma omp parallel private(i, j, k) shared(par, grid)
  //  {
  //#pragma omp for schedule(static)
  //#endif
  for(i=0; i<par->grid; i++)  // x directon
    for(j=0; j<par->grid; j++)  // y direction
      for(k=0; k<par->grid; k++)  // z direction
	{
	  //grid->phi_fifth[i][j][k] = phi_b(par, 1.0) + 1.0e-3*(double)ran3(&seed);  // symmetron
	  //grid->phi_fifth[i][j][k] = 1.0e-8*(double)ran3(&seed);  // symmetron
	  grid->q[i][j][k] = 1.0e-3*(double)ran3(&seed);
	}
  //#ifdef OPENMP
  //  }
  //#endif

}
