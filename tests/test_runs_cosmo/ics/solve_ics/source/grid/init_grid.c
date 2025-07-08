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
  fprintf(stderr, "Allocating...\n");
  fprintf(stderr, "1/8\n");
  alloc_double_3(&grid->delta,     par->grid, par->grid, par->grid);
  fprintf(stderr, "2/8\n");
  alloc_double_3(&grid->source,    par->grid, par->grid, par->grid);
  fprintf(stderr, "3/8\n");
  alloc_double_3(&grid->phi,       par->grid, par->grid, par->grid);
  fprintf(stderr, "4/8\n");
  alloc_double_3(&grid->lap_phi,   par->grid, par->grid, par->grid);
  fprintf(stderr, "5/8\n");
  alloc_double_3(&grid->force_1,   par->grid, par->grid, par->grid);
  fprintf(stderr, "6/8\n");
  alloc_double_3(&grid->force_2,   par->grid, par->grid, par->grid);
  fprintf(stderr, "7/8\n");
  alloc_double_3(&grid->force_3,   par->grid, par->grid, par->grid);
  fprintf(stderr, "8/8\n");
  alloc_double_3(&grid->q,         par->grid, par->grid, par->grid);
  
  // FFTW stuff:
  //============
  fprintf(stderr, "Calling init_fftw...\n");
  init_fftw(par, grid);
  fprintf(stderr, "done\n");

}
