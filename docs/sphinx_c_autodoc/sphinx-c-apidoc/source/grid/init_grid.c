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
void init_grid(struct params *par, struct grids *grid)
{

  // Allocation:
  //============
  // Density:
  //=========
  alloc_double_3(&grid->delta,     par->grid, par->grid, par->grid);

  // Newtonian fields:
  //==================
  alloc_double_3(&grid->source,    par->grid, par->grid, par->grid);
  alloc_double_3(&grid->phi,       par->grid, par->grid, par->grid);
  alloc_double_3(&grid->force_1,   par->grid, par->grid, par->grid);
  alloc_double_3(&grid->force_2,   par->grid, par->grid, par->grid);
  alloc_double_3(&grid->force_3,   par->grid, par->grid, par->grid);

  // Symmetron fields:
  //==================
  //alloc_double_3(&grid->phi_fifth, par->grid, par->grid, par->grid);
  //alloc_double_3(&grid->force_fifth_1,   par->grid, par->grid, par->grid);
  //alloc_double_3(&grid->force_fifth_2,   par->grid, par->grid, par->grid);
  //alloc_double_3(&grid->force_fifth_3,   par->grid, par->grid, par->grid);

  //alloc_double_3(&grid->phi_f,       par->grid, par->grid, par->grid);
  //alloc_double_3(&grid->phi_b,       par->grid, par->grid, par->grid);

  // MOND with mass fields:
  //=======================
  alloc_double_3(&grid->chi_bar_mm, par->grid, par->grid, par->grid);
  alloc_double_3(&grid->phi_bar_mm, par->grid, par->grid, par->grid);
  alloc_double_3(&grid->force_chi_bar_mm_1,   par->grid, par->grid, par->grid);
  alloc_double_3(&grid->force_chi_bar_mm_2,   par->grid, par->grid, par->grid);
  alloc_double_3(&grid->force_chi_bar_mm_3,   par->grid, par->grid, par->grid);
  alloc_double_3(&grid->force_phi_bar_mm_1,   par->grid, par->grid, par->grid);
  alloc_double_3(&grid->force_phi_bar_mm_2,   par->grid, par->grid, par->grid);
  alloc_double_3(&grid->force_phi_bar_mm_3,   par->grid, par->grid, par->grid);

  // FFTW stuff:
  //============
  //if(par->solver == 1)
  init_fftw(par, grid);
  fprintf(stderr, "done\n");

}
