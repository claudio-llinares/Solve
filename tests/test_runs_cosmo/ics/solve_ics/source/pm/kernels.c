#include <math.h>

#include "../solve.h"

//=====================
// tsc_kernel
//=====================
double tsc_kernel(double xg, double yg, double zg, double xp, double yp, double zp, struct params *par)
{
  // Calculates TSC kernel (3D).
  // Arguments:
  //       yg, yg, zg = position of the node that has that point.
  //       xp, yp, zp = positions of the point where you want the interpolation
  //
    
  double dist_x, dist_y, dist_z;

  double wx, wy, wz;

  double h = par->box / par->grid;
  
  // Get distance betwen particle and center of the node.
  dist_x = fabs(xp - xg);
  dist_y = fabs(yp - yg);
  dist_z = fabs(zp - zg);
  // Calculate weigth x.
  if( dist_x <= .5*h )
    wx = .75 - pow(dist_x / h, 2);
  else
    if(dist_x <= 1.5*h )
      wx = .5 * pow(1.5 - dist_x / h, 2);
    else
      wx = 0.0;
	  
  // Calculate weigth y.
  if( dist_y <= .5*h )
    wy = .75 - pow(dist_y / h, 2);
  else
    if(dist_y <= 1.5*h )
      wy = .5 * pow(1.5 - dist_y / h, 2);
    else
      wy = 0.0;

  // Calculate weigth z.
  if( dist_z <= .5*h )
    wz = .75 - pow(dist_z / h, 2);
  else
    if(dist_z <= 1.5*h )
      wz = .5 * pow(1.5 - dist_z / h, 2);
    else
      wz = 0.0;
  
  return wx * wy * wz;

}

//=====================
// tsc_kernel
//=====================
double cic_kernel(double xg, double yg, double zg, double xp, double yp, double zp, struct params *par)
{
  // Calculates CIC kernel (3D).
  // Arguments:
  //       yg, yg, zg = position of the node that has that point.
  //       xp, yp, zp = positions of the point where you want the interpolation
  //
    
  double dist_x, dist_y, dist_z;

  double wx, wy, wz;

  double h = par->box / par->grid;
  
  // Get distance betwen particle and center of the node.
  dist_x = fabs(xp - xg);
  dist_y = fabs(yp - yg);
  dist_z = fabs(zp - zg);
  // Calculate weigth x.
  if( dist_x <= h )
    wx = 1.0 - dist_x / h;
  else
    wx = 0.0;
	  
  // Calculate weigth y.
  if( dist_y <= h )
    wy = 1.0 - dist_y / h;
  else
    wy = 0.0;

  // Calculate weigth z.
  if( dist_z <= h )
    wz = 1.0 - dist_z / h;
  else
    wz = 0.0;
  
  return wx * wy * wz;

}


