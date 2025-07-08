#include <stdio.h>
#include <stdlib.h>

#include <fftw3.h>

#include "../libs/allocator.h"
#include "../libs/releaser.h"

#include "../solve.h"
#include "init_grid.h"
#include "init_fftw.h"

//====================
// init_grid
//====================
void clean_up_grid(struct params *par, struct grids *grid)
{

  // Releases all the stuff that was allocated by init_grid.

  // Allocation:
  //============
  free_double_3(grid->delta);
  free_double_3(grid->source);
  free_double_3(grid->phi);
  free_double_3(grid->force_1);
  free_double_3(grid->force_2);
  free_double_3(grid->force_3);

  // FFTW stuff:
  //============
  fftw_destroy_plan(grid->p_forward);
  fftw_destroy_plan(grid->p_backward);
  fftw_free(grid->data);
  fftw_cleanup();

#ifdef OPENMP
  fftw_cleanup_threads();
#endif
  

}
