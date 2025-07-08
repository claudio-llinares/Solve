/***
 * Interpolates a force vector from the grid to the particles usign CIC scheme.
 * Comments:
 *      - Details in the CIC scheme can be found in 
 *        Hockney R. W., Eastwood J. W., 1988, Computer simulation using particles.
**/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../solve.h"
#include "../grid/init_grid.h"
#include "../pm/init_parts.h"

//==============================
// interp_force_cic
//==============================
void interp_force_cic(struct params *par, struct grids *grid, struct particles *part, double ***force_1, double ***force_2, double ***force_3, double move_part)
{
  // Interpolate forces in force_ from the grid to the position
  // given by pos.
  //
  // Comments:
  //            - The indentification of the node that corresponds
  //              to pos is made here.
  //            - distance x, pos1 should be less that half of the grid size.

  int i, j, k;  // index of teh node

  double h = par->box/par->grid;
  
  int ip1, jp1, kp1;  // i plus 1, etc.
  
  double wx0, wxp;
  double wy0, wyp;
  double wz0, wzp;

  double p1, p2, p3;
  int i_part;

  double w[3][3][3];

  double force_i1, force_i2, force_i3;

  int sg_1, sg_2, sg_3;


  //static int first_problem = 1;
#ifdef OPENMP
#pragma omp parallel private(i_part, i, j, k, p1, p2, p3, wx0, wxp, wy0, wyp, wz0, wzp, ip1, jp1, kp1, w, force_i1, force_i2, force_i3) shared(par, part, grid, move_part)
  {
#pragma omp for schedule(static)
#endif
    for(i_part=0; i_part<par->nparts; i_part++)
      {
	i = (int)(part->pos_1[i_part]/h);  // index of node that contains the particle
	j = (int)(part->pos_2[i_part]/h);
	k = (int)(part->pos_3[i_part]/h);
	p1 = (part->pos_1[i_part]/h - i);  // fractional part of (x_particle/h)
	p2 = (part->pos_2[i_part]/h - j);
	p3 = (part->pos_3[i_part]/h - k);

	sg_1 = p1>0.5? 1 : -1;
	sg_2 = p2>0.5? 1 : -1;
	sg_3 = p3>0.5? 1 : -1;

	//===========================
	// Calculate weigths:
	//===================
	wxp = sg_1*(p1-0.5);
	wx0  = 1.0 - wxp;
	wyp = sg_2*(p2-0.5);	
	wy0  = 1.0 - wyp;
	wzp = sg_3*(p3-0.5);
	wz0  = 1.0 - wzp;

	// Calculate neibourgh nodes taking in account periodicity:
	//=========================================================
	ip1 = i+sg_1;  // this can be the node to the right or left according to sign.
	jp1 = j+sg_2;
	kp1 = k+sg_3;
	if(ip1<0)
	  ip1 = ip1 + par->grid;
	else
	  if(ip1 >= par->grid)
	    ip1 = ip1 - par->grid;
	if(jp1<0)
	  jp1 = jp1 + par->grid;
	else
	  if(jp1 >= par->grid)
	    jp1 = jp1 - par->grid;
	if(kp1<0)
	  kp1 = kp1 + par->grid;
	else
	  if(kp1 >= par->grid)
	    kp1 = kp1 - par->grid;
	
	// Calculate complete weights:
	//============================
	w[0][0][0] = wx0 * wy0 * wz0;
	w[0][0][1] = wx0 * wy0 * wzp;
	w[0][1][0] = wx0 * wyp * wz0;
	w[0][1][1] = wx0 * wyp * wzp;
	w[1][0][0] = wxp * wy0 * wz0;
	w[1][0][1] = wxp * wy0 * wzp;
	w[1][1][0] = wxp * wyp * wz0;
	w[1][1][1] = wxp * wyp * wzp;
	
	// Do the interpolation:
	//======================
	// x componet:
	force_i1 = 
	  force_1[  i][  j][  k] * w[0][0][0] +
	  force_1[  i][  j][kp1] * w[0][0][1] +
	  force_1[  i][jp1][  k] * w[0][1][0] +
	  force_1[  i][jp1][kp1] * w[0][1][1] +
	  force_1[ip1][  j][  k] * w[1][0][0] +
	  force_1[ip1][  j][kp1] * w[1][0][1] +
	  force_1[ip1][jp1][  k] * w[1][1][0] +
	  force_1[ip1][jp1][kp1] * w[1][1][1];

	force_i2 = 
	  force_2[  i][  j][  k] * w[0][0][0] +
	  force_2[  i][  j][kp1] * w[0][0][1] +
	  force_2[  i][jp1][  k] * w[0][1][0] +
	  force_2[  i][jp1][kp1] * w[0][1][1] +
	  force_2[ip1][  j][  k] * w[1][0][0] +
	  force_2[ip1][  j][kp1] * w[1][0][1] +
	  force_2[ip1][jp1][  k] * w[1][1][0] +
	  force_2[ip1][jp1][kp1] * w[1][1][1];

	force_i3 = 
	  force_3[  i][  j][  k] * w[0][0][0] +
	  force_3[  i][  j][kp1] * w[0][0][1] +
	  force_3[  i][jp1][  k] * w[0][1][0] +
	  force_3[  i][jp1][kp1] * w[0][1][1] +
	  force_3[ip1][  j][  k] * w[1][0][0] +
	  force_3[ip1][  j][kp1] * w[1][0][1] +
	  force_3[ip1][jp1][  k] * w[1][1][0] +
	  force_3[ip1][jp1][kp1] * w[1][1][1];

	part->vel_1[i_part] -= move_part*force_i1;
	part->vel_2[i_part] -= move_part*force_i2;
	part->vel_3[i_part] -= move_part*force_i3;

      }
#ifdef OPENMP
  }
#endif


  }
