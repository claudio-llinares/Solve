#include "stdio.h"
#include "stdlib.h"

#include <fftw3.h>

#include "../solve.h"
#include "init_grid.h"

//=======================
// init_fftw
//=======================
void init_fftw(struct params *par, struct grids *grid)
{

  // Initializes fftw stuff.

#ifdef OPENMP
  if( !fftw_init_threads())
    {
      fprintf(stderr, "Problems with fftw_init_threads.  Aborting...\n");
      exit(0);
    }
#endif
  (grid->data) = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * par->grid * par->grid * par->grid);
  fprintf(stderr, "Planing... ");
#ifdef OPENMP
  fftw_plan_with_nthreads(par->nthreads);
#endif
  fprintf(stderr, "forward... ");
  grid->p_forward = 
    fftw_plan_dft_3d(par->grid, 
		     par->grid, 
		     par->grid, 
		     (grid->data), (grid->data), 
		     FFTW_FORWARD, FFTW_ESTIMATE);
//		     FFTW_FORWARD, FFTW_MEASURE);

#ifdef OPENMP
  fftw_plan_with_nthreads(par->nthreads);
#endif
  fprintf(stderr, "backward... ");
  grid->p_backward = 
    fftw_plan_dft_3d(par->grid, 
		     par->grid, 
		     par->grid, 
		     (grid->data), (grid->data), 
		     FFTW_BACKWARD, FFTW_ESTIMATE);
//		     FFTW_BACKWARD, FFTW_MEASURE);
  
}
