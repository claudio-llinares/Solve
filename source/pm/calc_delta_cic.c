/*************************************************************************
* Calculates overdensity on the grid from the particle's positions using cloud 
* in cell (CIC) smoothing.
* Comments:
*      - Details in the CIC scheme can be found in 
*        Hockney R. W., Eastwood J. W., 1988, Computer simulation using particles.
****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../solve.h"
#include "../grid/init_grid.h"

#include "../libs/utils.h"


/**
 * Calculates overdensity on the grid from the particle's positions 
 * using cloud in cell (CIC) smoothing..
 * 
 * .. note:: This is a note for testign Sphinx.
 *
 * Parameters:
 *     par: par structure.
 *     unit: units structure.
 *     pos_1: X component of position vector.
 *     pos_2: X component of position vector.
 *     pos_3: X component of position vector.
 *     delta: to return overdensity.
 *
 * Returns:
 *    Overdensity in the parameter delta, which should be allocated.
 */
void calc_delta_cic(struct params *par, struct units *unit, double *pos_1, double *pos_2, double *pos_3, double ***delta)
{
    
  int i_part;
    
  int ngrid = par->grid;
  double h = par->box/par->grid;
  double factor_rho = pow3(ngrid)/par->nparts;

  int i, j, k;
  int ip1, jp1, kp1;  // i plus 1, etc.
  double p1, p2, p3;  // optimization
  int sg_1, sg_2, sg_3;

  double wx0, wxp;
  double wy0, wyp;
  double wz0, wzp;
  
  // Initialize overdensity:
  //========================
#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(delta, par)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<par->grid; i++)
      for(j=0; j<par->grid; j++)
	for(k=0; k<par->grid; k++)
	  delta[i][j][k] = 0.0;
#ifdef OPENMP
  }
#endif
  

  // Loop over the particles:
  //=========================
#ifdef OPENMP
#pragma omp parallel private(i_part, i, j, k, p1, p2, p3, sg_1, sg_2, sg_3, wx0, wxp, wy0, wyp, wz0, wzp, ip1, jp1, kp1) shared(pos_1, pos_2, pos_3, delta, h, ngrid)
  {
#pragma omp for schedule(static)
#endif
    for(i_part=0; i_part<par->nparts; i_part++)
      {
	i = (int)(pos_1[i_part]/h);  // index of node that contains the particle
	j = (int)(pos_2[i_part]/h);
	k = (int)(pos_3[i_part]/h);
	p1 = (pos_1[i_part]/h - i);  // fractional part of (x_particle/h)
	p2 = (pos_2[i_part]/h - j);
	p3 = (pos_3[i_part]/h - k);

	sg_1 = p1>0.5? 1 : -1;
	sg_2 = p2>0.5? 1 : -1;
	sg_3 = p3>0.5? 1 : -1;

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
	ip1 = ( i + sg_1 + ngrid ) % ngrid;
	jp1 = ( j + sg_2 + ngrid ) % ngrid;
	kp1 = ( k + sg_3 + ngrid ) % ngrid;
	
	// Add the particle to all the neibourgh nodes:
	//=============================================
#ifdef OPENMP
#pragma omp atomic
#endif
	delta[  i][  j][  k] += wx0 * wy0 * wz0;
#ifdef OPENMP
#pragma omp atomic
#endif
	delta[  i][  j][kp1] += wx0 * wy0 * wzp;
#ifdef OPENMP
#pragma omp atomic
#endif
	delta[  i][jp1][  k] += wx0 * wyp * wz0;
#ifdef OPENMP
#pragma	omp atomic
#endif
	delta[  i][jp1][kp1] += wx0 * wyp * wzp;

#ifdef OPENMP
#pragma omp atomic
#endif
	delta[ip1][  j][  k] += wxp * wy0 * wz0;
#ifdef OPENMP
#pragma omp atomic
#endif
	delta[ip1][  j][kp1] += wxp * wy0 * wzp;
#ifdef OPENMP
#pragma omp atomic
#endif
	delta[ip1][jp1][  k] += wxp * wyp * wz0;
#ifdef OPENMP
#pragma omp atomic
#endif
	delta[ip1][jp1][kp1] += wxp * wyp * wzp;
	
      }
#ifdef OPENMP
  }
#endif	


#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(delta, par, factor_rho)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<par->grid; i++)
      for(j=0; j<par->grid; j++)
	for(k=0; k<par->grid; k++)
	  delta[i][j][k] = delta[i][j][k]*factor_rho - 1.0;
#ifdef OPENMP
  }
#endif
  
}


