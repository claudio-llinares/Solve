#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../solve.h"

//=========================
// interp_force
//=========================
void interp_force(struct params *par, double ***force_1, double ***force_2, double ***force_3, double pos_1, double pos_2, double pos_3, double *force_i1, double *force_i2, double *force_i3, double (*kernel)(double, double, double, double, double, double, struct params *))
{
  // Interpolate forces in force_ from the grid to the position
  // given by pos.
  // 
  // Arguments:
  //            - kernel = pointer to the definition of the kernel.
  // Comments:
  //            - The indentification of the node that corresponds
  //              to pos is made here.
  //            - distance x, pos1 should be less that half of the grid size.

  int i, j, k;  // index of teh node
  double x, y, z;  // position of the node

  double h = par->box/par->grid;
  double ho2 = h/2.0;
  int gm1 = par->grid - 1;

  int ip1, im1, jp1, jm1, kp1, km1;  // i plus 1, etc.

  double grid_over_box = par->grid/par->box;

  //uble xph, xmh, yph, ymh, zph, zmh;

  double dist_xm, dist_x0, dist_xp;
  double dist_ym, dist_y0, dist_yp;
  double dist_zm, dist_z0, dist_zp;
  double wxm, wx0, wxp;
  double wym, wy0, wyp;
  double wzm, wz0, wzp;

  //static int first_problem = 1;

  // Identify node in which the particle is:
  //========================================
  i = (int)floor(pos_1*grid_over_box);  // index of teh node (the rho grid)
  j = (int)floor(pos_2*grid_over_box);
  k = (int)floor(pos_3*grid_over_box);
  x = ((double)i/grid_over_box + ho2);  // position of the node 
  y = ((double)j/grid_over_box + ho2);
  z = ((double)k/grid_over_box + ho2);

  //exit(0);
  /* if( ( x-pos_1 >= ho2 || */
  /* 	y-pos_2 >= ho2 || */
  /* 	z-pos_3 >= ho2 || */
  /* 	i<0 || j<0 || k<0 || */
  /* 	i>gm1 || j>gm1 || k>gm1 ) &&  */
  /*     first_problem ) */
  /*   { */
  /*     fprintf(stderr, "interp_force:  there is a problem...\n%d  %d  %d\n%lg  %lg  %lg\n%lg  %lg  %lg\n%20.17lg  %20.17lg  %20.17lg\n%20.17lg\n%lg  %lg  %lg\n", */
  /* 	      i, j, k, */
  /* 	      pos_1, pos_2, pos_3, */
  /* 	      x, y, z, */
  /* 	      x-pos_1, y-pos_2, z-pos_3, */
  /* 	      ho2, */
  /* 	      pos_1/par->box * par->grid, */
  /* 	      pos_2/par->box * par->grid, */
  /* 	      pos_3/par->box * par->grid); */
  /*     fprintf(stderr, "box = %lf\n", par->box); */
  /*     fprintf(stderr, "h   = %lf\n", h); */
  /*     fprintf(stderr, "ho2 = %lf\n", ho2); */
  /*     fprintf(stderr, "No more messages like this one will be printed.\n"); */
  /*     first_problem = 0; */
  /*     //exit(0); */
  /*   }      */

  // Calculate weigths:
  //===================
  dist_xm = fabs(pos_1-(x-h));
  dist_x0 = fabs(pos_1-(x  ));
  dist_xp = fabs(pos_1-(x+h));
  wxm = (dist_xm<=h ? (1.0 - dist_xm/h) : 0.0);
  wx0 = (dist_x0<=h ? (1.0 - dist_x0/h) : 0.0);
  wxp = (dist_xp<=h ? (1.0 - dist_xp/h) : 0.0);

  dist_ym = fabs(pos_2-(y-h));
  dist_y0 = fabs(pos_2-(y  ));
  dist_yp = fabs(pos_2-(y+h));
  wym = (dist_ym<=h ? (1.0 - dist_ym/h) : 0.0);
  wy0 = (dist_y0<=h ? (1.0 - dist_y0/h) : 0.0);
  wyp = (dist_yp<=h ? (1.0 - dist_yp/h) : 0.0);

  dist_zm = fabs(pos_3-(z-h));
  dist_z0 = fabs(pos_3-(z  ));
  dist_zp = fabs(pos_3-(z+h));
  wzm = (dist_zm<=h ? (1.0 - dist_zm/h) : 0.0);
  wz0 = (dist_z0<=h ? (1.0 - dist_z0/h) : 0.0);
  wzp = (dist_zp<=h ? (1.0 - dist_zp/h) : 0.0);

  // Calculate neibourgh nodes taking in account periodicity:
  //=========================================================
  ip1 = (i==(gm1) ? 0     : (i+1));
  im1 = (i==(0)   ? (gm1) : (i-1));
  jp1 = (j==(gm1) ? 0     : (j+1));
  jm1 = (j==(0)   ? (gm1) : (j-1));
  kp1 = (k==(gm1) ? 0     : (k+1));
  km1 = (k==(0)   ? (gm1) : (k-1));

  // Initialize force:
  //==================
  *force_i1 = 0.0;
  *force_i2 = 0.0;
  *force_i3 = 0.0;

  // Do the interpolation:
  //======================
  // x componet:
  *force_i1 += force_1[  i][  j][  k] * wx0 * wy0 * wz0;
  *force_i1 += force_1[  i][  j][km1] * wx0 * wy0 * wzm;
  *force_i1 += force_1[  i][  j][kp1] * wx0 * wy0 * wzm;
  *force_i1 += force_1[  i][jm1][  k] * wx0 * wym * wz0;
  *force_i1 += force_1[  i][jm1][km1] * wx0 * wym * wzm;
  *force_i1 += force_1[  i][jm1][kp1] * wx0 * wym * wzp;
  *force_i1 += force_1[  i][jp1][  k] * wx0 * wyp * wz0;
  *force_i1 += force_1[  i][jp1][km1] * wx0 * wyp * wzm;
  *force_i1 += force_1[  i][jp1][kp1] * wx0 * wyp * wzp;
  
  *force_i1 += force_1[im1][  j][  k] * wxm * wy0 * wz0;
  *force_i1 += force_1[im1][  j][km1] * wxm * wy0 * wzm;
  *force_i1 += force_1[im1][  j][kp1] * wxm * wy0 * wzp;
  *force_i1 += force_1[im1][jm1][  k] * wxm * wym * wz0;
  *force_i1 += force_1[im1][jm1][km1] * wxm * wym * wzm;
  *force_i1 += force_1[im1][jm1][kp1] * wxm * wym * wzp;
  *force_i1 += force_1[im1][jp1][  k] * wxm * wyp * wz0;
  *force_i1 += force_1[im1][jp1][km1] * wxm * wyp * wzm;
  *force_i1 += force_1[im1][jp1][kp1] * wxm * wyp * wzp;

  *force_i1 += force_1[ip1][  j][  k] * wxp * wy0 * wz0;
  *force_i1 += force_1[ip1][  j][km1] * wxp * wy0 * wzm;
  *force_i1 += force_1[ip1][  j][kp1] * wxp * wy0 * wzp;
  *force_i1 += force_1[ip1][jm1][  k] * wxp * wym * wz0;
  *force_i1 += force_1[ip1][jm1][km1] * wxp * wym * wzm;
  *force_i1 += force_1[ip1][jm1][kp1] * wxp * wym * wzp;
  *force_i1 += force_1[ip1][jp1][  k] * wxp * wyp * wz0;
  *force_i1 += force_1[ip1][jp1][km1] * wxp * wyp * wzm;
  *force_i1 += force_1[ip1][jp1][kp1] * wxp * wyp * wzp;
  // y componet:
  *force_i2 += force_2[  i][  j][  k] * wx0 * wy0 * wz0;
  *force_i2 += force_2[  i][  j][km1] * wx0 * wy0 * wzm;
  *force_i2 += force_2[  i][  j][kp1] * wx0 * wy0 * wzp;
  *force_i2 += force_2[  i][jm1][  k] * wx0 * wym * wz0;
  *force_i2 += force_2[  i][jm1][km1] * wx0 * wym * wzm;
  *force_i2 += force_2[  i][jm1][kp1] * wx0 * wym * wzp;
  *force_i2 += force_2[  i][jp1][  k] * wx0 * wyp * wz0;
  *force_i2 += force_2[  i][jp1][km1] * wx0 * wyp * wzm;
  *force_i2 += force_2[  i][jp1][kp1] * wx0 * wyp * wzp;
  
  *force_i2 += force_2[im1][  j][  k] * wxm * wy0 * wz0;
  *force_i2 += force_2[im1][  j][km1] * wxm * wy0 * wzm;
  *force_i2 += force_2[im1][  j][kp1] * wxm * wy0 * wzp;
  *force_i2 += force_2[im1][jm1][  k] * wxm * wym * wz0;
  *force_i2 += force_2[im1][jm1][km1] * wxm * wym * wzm;
  *force_i2 += force_2[im1][jm1][kp1] * wxm * wym * wzp;
  *force_i2 += force_2[im1][jp1][  k] * wxm * wyp * wz0;
  *force_i2 += force_2[im1][jp1][km1] * wxm * wyp * wzm;
  *force_i2 += force_2[im1][jp1][kp1] * wxm * wyp * wzp;

  *force_i2 += force_2[ip1][  j][  k] * wxp * wy0 * wz0;
  *force_i2 += force_2[ip1][  j][km1] * wxp * wy0 * wzm;
  *force_i2 += force_2[ip1][  j][kp1] * wxp * wy0 * wzp;
  *force_i2 += force_2[ip1][jm1][  k] * wxp * wym * wz0;
  *force_i2 += force_2[ip1][jm1][km1] * wxp * wym * wzm;
  *force_i2 += force_2[ip1][jm1][kp1] * wxp * wym * wzp;
  *force_i2 += force_2[ip1][jp1][  k] * wxp * wyp * wz0;
  *force_i2 += force_2[ip1][jp1][km1] * wxp * wyp * wzm;
  *force_i2 += force_2[ip1][jp1][kp1] * wxp * wyp * wzp;
  // z componet:
  *force_i3 += force_3[  i][  j][  k] * wx0 * wy0 * wz0;
  *force_i3 += force_3[  i][  j][km1] * wx0 * wy0 * wzm;
  *force_i3 += force_3[  i][  j][kp1] * wx0 * wy0 * wzp;
  *force_i3 += force_3[  i][jm1][  k] * wx0 * wym * wz0;
  *force_i3 += force_3[  i][jm1][km1] * wx0 * wym * wzm;
  *force_i3 += force_3[  i][jm1][kp1] * wx0 * wym * wzp;
  *force_i3 += force_3[  i][jp1][  k] * wx0 * wyp * wz0;
  *force_i3 += force_3[  i][jp1][km1] * wx0 * wyp * wzm;
  *force_i3 += force_3[  i][jp1][kp1] * wx0 * wyp * wzp;
  
  *force_i3 += force_3[im1][  j][  k] * wxm * wy0 * wz0;
  *force_i3 += force_3[im1][  j][km1] * wxm * wy0 * wzm;
  *force_i3 += force_3[im1][  j][kp1] * wxm * wy0 * wzp;
  *force_i3 += force_3[im1][jm1][  k] * wxm * wym * wz0;
  *force_i3 += force_3[im1][jm1][km1] * wxm * wym * wzm;
  *force_i3 += force_3[im1][jm1][kp1] * wxm * wym * wzp;
  *force_i3 += force_3[im1][jp1][  k] * wxm * wyp * wz0;
  *force_i3 += force_3[im1][jp1][km1] * wxm * wyp * wzm;
  *force_i3 += force_3[im1][jp1][kp1] * wxm * wyp * wzp;

  *force_i3 += force_3[ip1][  j][  k] * wxp * wy0 * wz0;
  *force_i3 += force_3[ip1][  j][km1] * wxp * wy0 * wzm;
  *force_i3 += force_3[ip1][  j][kp1] * wxp * wy0 * wzp;
  *force_i3 += force_3[ip1][jm1][  k] * wxp * wym * wz0;
  *force_i3 += force_3[ip1][jm1][km1] * wxp * wym * wzm;
  *force_i3 += force_3[ip1][jm1][kp1] * wxp * wym * wzp;
  *force_i3 += force_3[ip1][jp1][  k] * wxp * wyp * wz0;
  *force_i3 += force_3[ip1][jp1][km1] * wxp * wyp * wzm;
  *force_i3 += force_3[ip1][jp1][kp1] * wxp * wyp * wzp;

}
