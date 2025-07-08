/** 
 * Prologation routines which interpolate a field from a coarse grid back to a fine grid.
 * 
 * Comments:
 *     - See Trottenberg, Oosterlee, Schuller, "Multigrid" for details.
 *     - These routines are needed by the multigrid solver to move between grids.
 *     - There are several routines that use different interpolation order.  Higher order is not necessarily better.
 *     - The only public routine is ``prolonge``.  The other routines are accessed through it.  Uncomment the call you want to use.
 * 
 **/
#include <math.h>

#define pow2(x) ((x)*(x))

extern void prolonge_first_order(double ***old, double ***new, int n_old, int n_new, double h_old, double h_new);
extern void prolonge_second_order(double ***old, double ***new, int n_old, int n_new, double h_old, double h_new);
extern void prolonge_trilinear(double ***old, double ***new, int n_old, int n_new, double h_old, double h_new);


void prolonge(double ***old, double ***new, int n_old, double h_old)
{

  int n_new;
  double h_new;

  n_new = n_old*2;
  h_new = h_old/2.0;

  //prolonge_first_order(old, new, n_old, n_new, h_old, h_new);
  prolonge_second_order(old, new, n_old, n_new, h_old, h_new);
  //prolonge_trilinear(old, new, n_new, n_old, h_old, h_new);


}

//======================================
// prolonge_first_order
//======================================
void prolonge_first_order(double ***old, double ***new, int n_old, int n_new, double h_old, double h_new)
{
  //
  // Prolonge values in array old and put them in new.
  // Size arrays should be ok.
  // g_old is the number of the old grid (coarse).
  //
  // Arguments:
  //           - n_new = size of new grid.
  //           - h_new = spacing of new grid.
  // Comments:
  //           - We do not touch the boundary.
  //
  
  
  int i, j, k;
  long i_old, j_old, k_old;
  double x_old, y_old, z_old, x_new, y_new, z_new;
  double dold_dx, dold_dy, dold_dz;

  int i_old_p1, i_old_m1, j_old_p1, j_old_m1, k_old_p1, k_old_m1;

  int n_oldm1;  // optimization

  //n = gg.n[g_old-1];
  n_oldm1 = n_old-1; //gg.n[g_old] - 1;

  // Limits are in a way that we only touch the inner-inner nodes
  // of the lower grid.
#ifdef OPENMP
#pragma omp parallel private(i, j, k, i_old, j_old, k_old, x_old, y_old, z_old, x_new, y_new, z_new, i_old_p1, i_old_m1, j_old_p1, j_old_m1, k_old_p1, k_old_m1, dold_dx, dold_dy, dold_dz) shared(n_new, n_oldm1, old, new, h_old, h_new)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<n_new; i+=1)   // Fills all the values in the grid.
      {
	for(j=0; j<n_new; j+=1)
	  {
	    for(k=0; k<n_new; k+=1)
	      {
		// indexes/positions in both grids:
		i_old = (int)floor(i/2.0);
		j_old = (int)floor(j/2.0);	    
		k_old = (int)floor(k/2.0);	    
		x_old = (double)i_old * h_old + h_old/2.0;
		y_old = (double)j_old * h_old + h_old/2.0;
		z_old = (double)k_old * h_old + h_old/2.0;
		x_new = (double)i * h_new + h_new/2.0;
		y_new = (double)j * h_new + h_new/2.0;
		z_new = (double)k * h_new + h_new/2.0;
		
		// periodic stuff:
		i_old_p1 = (i_old==(n_oldm1) ? 0         : (i_old+1));
		i_old_m1 = (i_old==(0)       ? (n_oldm1) : (i_old-1));
		j_old_p1 = (j_old==(n_oldm1) ? 0         : (j_old+1));
		j_old_m1 = (j_old==(0)       ? (n_oldm1) : (j_old-1));
		k_old_p1 = (k_old==(n_oldm1) ? 0         : (k_old+1));
		k_old_m1 = (k_old==(0)       ? (n_oldm1) : (k_old-1));
		
		// First derivatives:
		dold_dx = ( old[i_old_p1][j_old][k_old] - 
			    old[i_old_m1][j_old][k_old] ) / (2.0 * h_old);
		dold_dy = ( old[i_old][j_old_p1][k_old] - 
			    old[i_old][j_old_m1][k_old] ) / (2.0 * h_old);
		dold_dz = ( old[i_old][j_old][k_old_p1] - 
			    old[i_old][j_old][k_old_m1] ) / (2.0 * h_old);
				
		// First order:
		new[i][j][k] = old[i_old][j_old][k_old] + 
		  dold_dx * (x_new - x_old) +
		  dold_dy * (y_new - y_old) +
		  dold_dz * (z_new - z_old);
	      }
	  }
      }
#ifdef OPENMP
  }
#endif

}

//======================================
// prolonge_second_order
//======================================
void prolonge_second_order(double ***old, double ***new, int n_old, int n_new, double h_old, double h_new)
{
  //
  // Prolonge values in array old and put them in new.
  // Size arrays should be ok.
  // g_old is the number of the old grid (coarse).
  //
  // Arguments:
  //           - n_new = size of new grid.
  //           - h_new = spacing of new grid.
  // Comments:
  //           - We do not touch the boundary.
  //
  
  
  int i, j, k;
  long i_old, j_old, k_old;
  double x_old, y_old, z_old, x_new, y_new, z_new;
  double dold_dx, dold_dy, dold_dz;
  double ddold_dxdx, ddold_dydy, ddold_dzdz, ddold_dxdy, ddold_dxdz, ddold_dydz;

  int i_old_p1, i_old_m1, j_old_p1, j_old_m1, k_old_p1, k_old_m1;

  int n_oldm1;  // optimization

  //n = gg.n[g_old-1];
  n_oldm1 = n_old-1; //gg.n[g_old] - 1;

  // Limits are in a way that we only touch the inner-inner nodes
  // of the lower grid.
#ifdef OPENMP
#pragma omp parallel private(i, j, k, i_old, j_old, k_old, x_old, y_old, z_old, x_new, y_new, z_new, i_old_p1, i_old_m1, j_old_p1, j_old_m1, k_old_p1, k_old_m1, dold_dx, dold_dy, dold_dz, ddold_dxdx, ddold_dydy, ddold_dzdz, ddold_dxdy, ddold_dxdz, ddold_dydz) shared(n_new, n_oldm1, old, new, h_old, h_new)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<n_new; i+=1)   // Fills all the values in the grid.
      {
	for(j=0; j<n_new; j+=1)
	  {
	    for(k=0; k<n_new; k+=1)
	      {
		// indexes/positions in both grids:
		i_old = (int)floor(i/2.0);
		j_old = (int)floor(j/2.0);	    
		k_old = (int)floor(k/2.0);	    
		x_old = (double)i_old * h_old + h_old/2.0;
		y_old = (double)j_old * h_old + h_old/2.0;
		z_old = (double)k_old * h_old + h_old/2.0;
		x_new = (double)i * h_new + h_new/2.0;
		y_new = (double)j * h_new + h_new/2.0;
		z_new = (double)k * h_new + h_new/2.0;
		
		// periodic stuff:
		i_old_p1 = (i_old==(n_oldm1) ? 0         : (i_old+1));
		i_old_m1 = (i_old==(0)       ? (n_oldm1) : (i_old-1));
		j_old_p1 = (j_old==(n_oldm1) ? 0         : (j_old+1));
		j_old_m1 = (j_old==(0)       ? (n_oldm1) : (j_old-1));
		k_old_p1 = (k_old==(n_oldm1) ? 0         : (k_old+1));
		k_old_m1 = (k_old==(0)       ? (n_oldm1) : (k_old-1));
		
		// First derivatives:
		dold_dx = ( old[i_old_p1][j_old][k_old] - 
			    old[i_old_m1][j_old][k_old] ) / (2.0 * h_old);
		dold_dy = ( old[i_old][j_old_p1][k_old] - 
			    old[i_old][j_old_m1][k_old] ) / (2.0 * h_old);
		dold_dz = ( old[i_old][j_old][k_old_p1] - 
			    old[i_old][j_old][k_old_m1] ) / (2.0 * h_old);
		
		// Second derivatives:
		ddold_dxdx = ( old[i_old_p1][j_old][k_old] +
			       old[i_old_m1][j_old][k_old] -
			       2.0*old[i_old][j_old][k_old] ) / 
		  pow(h_old,2);
		
		ddold_dydy = ( old[i_old][j_old_p1][k_old] +
			       old[i_old][j_old_m1][k_old] -
			       2.0*old[i_old][j_old][k_old] ) / 
		  pow(h_old,2);
		
		ddold_dzdz = ( old[i_old][j_old][k_old_p1] +
			       old[i_old][j_old][k_old_m1] -
			       2.0*old[i_old][j_old][k_old] ) / 
		  pow(h_old,2);
		
		ddold_dxdy = ( old[i_old_p1][j_old_p1][k_old] - 
			       old[i_old_p1][j_old_m1][k_old] -
			       old[i_old_m1][j_old_p1][k_old] + 
			       old[i_old_m1][j_old_m1][k_old] ) / 
		  (4.0*pow2(h_old));
		
		ddold_dxdz = ( old[i_old_p1][j_old][k_old_p1] - 
			       old[i_old_p1][j_old][k_old_m1] -
			       old[i_old_m1][j_old][k_old_p1] + 
			       old[i_old_m1][j_old][k_old_m1] ) / 
		  (4.0*pow2(h_old));
		
		ddold_dydz = ( old[i_old][j_old_p1][k_old_p1] - 
			       old[i_old][j_old_m1][k_old_p1] -
			       old[i_old][j_old_p1][k_old_m1] + 
			       old[i_old][j_old_m1][k_old_m1] ) / 
		  (4.0*pow2(h_old));
		
		// First order:
		new[i][j][k] = old[i_old][j_old][k_old] + 
		  dold_dx * (x_new - x_old) +
		  dold_dy * (y_new - y_old) +
		  dold_dz * (z_new - z_old);
		
		// Add second order:
		if(1)
		  {
		    new[i][j][k] += ddold_dxdx * pow2(x_new - x_old) +
		      ddold_dydy * pow2(y_new - y_old) +
		      ddold_dzdz * pow2(z_new - z_old) +
		      2.0*ddold_dxdy * (x_new - x_old)*(y_new - y_old) +
		      2.0*ddold_dxdz * (x_new - x_old)*(z_new - z_old) +
		      2.0*ddold_dydz * (y_new - y_old)*(z_new - z_old);
		  }
	      }
	  }
      }
#ifdef OPENMP
  }
#endif

}
//======================================
// prolonge_trilinear
//======================================
void prolonge_trilinear(double ***old, double ***new, int n_old, int n_new, double h_old, double h_new)
{
  //
  // Prolonge values in array old and put them in new.
  // Size arrays should be ok.
  // g_old is the number of the old grid (coarse).
  //
  // Arguments:
  //           - n_new = size of new grid.
  //           - h_new = spacing of new grid.
  // Comments:
  //            - We do not touch the boundary.
  //

  int i, j, k;
  long i_old, j_old, k_old;
   
  int i_old_p1, j_old_p1, k_old_p1;

  int n_oldm1;  // optimization

  //n = gg.n[g_old-1];
  n_oldm1 = n_old-1; //gg.n[g_old] - 1;

  // Limits are in a way that we only touch the inner-inner nodes
  // of the lower grid.
  for(i=0; i<n_new; i+=1)   // Fills all the values in the grid.
    {
      for(j=0; j<n_new; j+=1)
	{
	  for(k=0; k<n_new; k+=1)
	    {
	      // indexes/positions in both grids:
	      i_old = (int)floor(i/2.0);
	      j_old = (int)floor(j/2.0);	    
	      k_old = (int)floor(k/2.0);

	      // periodic stuff:
	      i_old_p1 = (i_old==(n_oldm1) ? 0         : (i_old+1));
	      j_old_p1 = (j_old==(n_oldm1) ? 0         : (j_old+1));
	      k_old_p1 = (k_old==(n_oldm1) ? 0         : (k_old+1));
	      
	      if(i_old==i && j_old==j && k_old==k)
		new[i][j][k] = old[i_old][j_old][k_old];
	      else
		if(k_old==k && i_old==i)
		  new[i][j][k] = (old[i_old][j_old   ][k_old] + 
				  old[i_old][j_old_p1][k_old] ) / 2.0;
	      else
		if(k_old==k && j_old==j)
		  new[i][j][k] = (old[i_old   ][j_old][k_old] + 
				  old[i_old_p1][j_old][k_old] ) / 2.0;
	      else
		if(i_old==i && j_old==j)
		  new[i][j][k] = (old[i_old][j_old][k_old   ] + 
				  old[i_old][j_old][k_old_p1] ) / 2.0;
	      else
		if(k_old==k)
		  new[i][j][k] = (old[i_old   ][j_old   ][k_old] + 
				  old[i_old_p1][j_old   ][k_old] +
				  old[i_old   ][j_old_p1][k_old] + 
				  old[i_old_p1][j_old_p1][k_old] ) / 4.0;
	      else
		if(i_old==i)
		  new[i][j][k] = (old[i_old][j_old   ][k_old   ] + 
				  old[i_old][j_old   ][k_old_p1] +
				  old[i_old][j_old_p1][k_old   ] + 
				  old[i_old][j_old_p1][k_old_p1] ) / 4.0;
	      else
		if(j_old==j)
		  new[i][j][k] = (old[i_old][j_old][k_old] + 
				  old[i_old][j_old][k_old_p1] +
				  old[i_old_p1][j_old][k_old] + 
				  old[i_old_p1][j_old][k_old_p1] ) / 4.0;
		else
		  new[i][j][k] = (old[i_old   ][j_old   ][k_old   ] + 
				  old[i_old_p1][j_old   ][k_old   ] +
				  old[i_old   ][j_old_p1][k_old   ] + 
				  old[i_old   ][j_old   ][k_old_p1] +
				  old[i_old_p1][j_old_p1][k_old   ] + 
				  old[i_old_p1][j_old   ][k_old_p1] +
				  old[i_old   ][j_old_p1][k_old_p1] + 
				  old[i_old_p1][j_old_p1][k_old_p1]) / 8.0;
	    }
	}
    }

}
