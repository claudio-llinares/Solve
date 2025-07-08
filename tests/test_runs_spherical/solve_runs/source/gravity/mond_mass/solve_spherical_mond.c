/***
 * Routines for solving for the spherical mond potential given by
 * the standard poisson's equation.
 * The calculation of the source and everything is made here.
 * - See appendix in Llinares' thesis 2010.
 ***/

#include<stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../libs/allocator.h"
#include "../libs/releaser.h"

#include "../solve_poisson_fftw.h"
#include "../solve_growth_pnb.h"
#include "../get_lap.h"

#include "../get_forces.h"

// faster way to calculate power of 2, 3 and 4:
//=============================================
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define pow5(x) ((x)*(x)*(x)*(x)*(x))
#define pow6(x) ((x)*(x)*(x)*(x)*(x)*(x))


/***
 * This routine comes from the thesis.  It has to be merged into the new 
 * code.  It is not included in the Makefile.
 * 
 * The routine solves an equation for the spherical mond potential which can 
 * be obtained by calculating the derivative of the integral (yes, that's 
 * what I wrote) of the MOND equation.
 * The method is:
 *                 - get newtonian forces.
 *                 - get spherical mond forces.
 *                 - take divergence.
 *                 - solve poisson's equation with that source.
 * 
 * Arguments:
 *           - phi = newtonian potential.
 *           - a = expansion factor.
 * Returns:
 *           - phi = spherical mondian potential.
 * Comments:
 *           - No units (they come with whatever you put as input).
 *           - THIS ROUTINE REWITES PHI.
 **/
int solve_spherical_mond(double ***phi, double a, fftw_complex *data, fftw_plan p_forward, fftw_plan p_backward, double power_a0)
{


  double ***dphi_dx, ***dphi_dy, ***dphi_dz;
  int i, j, k;

  int ip1, im1, jp1, jm1, kp1, km1;

  double ***div_force;

  double h;


  fprintf(stderr, "solve_spherical_mond: calculating spherical MONDian potential...\n");

  // Allocate grids:
  //================
  alloc_double_3(&dphi_dx, par.grid, par.grid, par.grid);
  alloc_double_3(&dphi_dy, par.grid, par.grid, par.grid);
  alloc_double_3(&dphi_dz, par.grid, par.grid, par.grid);
  alloc_double_3(&div_force, par.grid, par.grid, par.grid);
  h = par.box / par.grid;

  // Get laplacian of spherical potential:
  //======================================
  calc_lap_mond(phi, div_force, a, par.power_a0);

  // Solve poisson's equation with the divergence as source:
  //========================================================
  solve_poisson_fftw(div_force, phi, par.box, par.grid, data, p_forward, p_backward);

  // Deallocate grids:
  //==================
  free_double_3(dphi_dx, par.grid, par.grid, par.grid);
  free_double_3(dphi_dy, par.grid, par.grid, par.grid);
  free_double_3(dphi_dz, par.grid, par.grid, par.grid);
  free_double_3(div_force, par.grid, par.grid, par.grid);

  return 0;
    
}
