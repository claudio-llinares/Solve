/** 
 * Restriction routines which interpolate a field from a fina grid to a coarse grid (with half the number of points per dimension).
 * 
 * Comments:
 *     - See Trottenberg, Oosterlee, Schuller, "Multigrid" for details.
 *     - These routines are needed by the multigrid solver to move between grids.
 *     - There are several routines that use different interpolation order.  Higher order is not necessarily better.
 *     - The only public routine is ``restrict_grid``.  The other routines are accessed through it.  Uncomment the call you want to use.
 * 
 **/

#include <math.h>

void restrict_grid_cospatial(double ***old, double ***new, long n_old);
void restrict_grid_original(double ***old, double ***new, long n_old);

//======================================
// restrict
//======================================
void restrict_grid(double ***old, double ***new, long n_old)
{
  
  restrict_grid_original(old, new, n_old);
  //restrict_grid_cospatial(old, new, n_old);

}

//=========================
// restrict_grid_cospatial
//=========================
void restrict_grid_cospatial(double ***old, double ***new, long n_old)
{
  // Restricts values in array old and put them in new.
  // Size of the arrays should be ok.
  // n_old = size of old grid (fine).
  //
  // old = fine
  // new = coarse
  
  int i, j, k;
  int i_new, j_new, k_new;

  int ip1, im1, jp1, jm1, kp1, km1;

  int n_oldm1;  // optimization
  
  //n_oldm1 = gg.n[g_old]-1;
  n_oldm1 = n_old - 1;

  // Limits are in a way that we only touch the inner
  // nodes of the upper grid.
  // -4 because we advance +2
  // Loops are for the fine grid.
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, im1, jp1, jm1, kp1, km1, i_new, j_new, k_new) shared(n_oldm1, new, old)
  {
#pragma omp for schedule(static)
#endif
  for(i=0; i<=n_oldm1; i+=2)   // Fills all the values in the grid.
    {
      for(j=0; j<=n_oldm1; j+=2)
	{
	  for(k=0; k<=n_oldm1; k+=2)
	    {
	      // periodic stuff:
	      ip1 = (i==(n_oldm1) ? 0       : (i+1));
	      im1 = (i==(0)       ? n_oldm1 : (i-1));
	      jp1 = (j==(n_oldm1) ? 0       : (j+1));
	      jm1 = (j==(0)       ? n_oldm1 : (j-1));
	      kp1 = (k==(n_oldm1) ? 0       : (k+1));
	      km1 = (k==(0)       ? n_oldm1 : (k-1));
	      

	      i_new = (int)floor(i/2.0);
	      j_new = (int)floor(j/2.0);	    
	      k_new = (int)floor(k/2.0);	    

	      new[i_new][j_new][k_new] = 
		8.0*old[i  ][j  ][k  ];

	      new[i_new][j_new][k_new] += 
		4.0 * ( old[ip1][j  ][k  ] +
		        old[i  ][jp1][k  ] +
		        old[i  ][j  ][kp1] + 
		        old[im1][j  ][k  ] +
		        old[i  ][jm1][k  ] +
		        old[i  ][j  ][km1] );

	      new[i_new][j_new][k_new] += 
		2.0 * ( old[ip1][jp1][k  ] + 
			old[ip1][jm1][k  ] +
			old[im1][jp1][k  ] +
			old[im1][jm1][k  ] +
			old[ip1][jp1][kp1] + 
			old[ip1][jm1][kp1] +
			old[im1][jp1][kp1] +
			old[im1][jm1][kp1] +
			old[ip1][jp1][km1] + 
			old[ip1][jm1][km1] +
			old[im1][jp1][km1] +
			old[im1][jm1][km1] );

	      new[i_new][j_new][k_new] += 
		old[ip1][jp1][kp1] + 
		old[ip1][jp1][km1] + 
		old[ip1][jm1][kp1] + 
		old[ip1][jm1][km1] + 
		old[im1][jp1][kp1] + 
		old[im1][jp1][km1] + 
		old[im1][jm1][kp1] + 
		old[im1][jm1][km1];

	      new[i_new][j_new][k_new] /= 64.0;

	      //fprintf(stderr, "%d  %d       ", i_new, j_new);
	    }
	}
      //fprintf(stderr, "\n");
    }
#ifdef OPENMP
  }
#endif
  
}

//=========================
// restrict original
//=========================
void restrict_grid_original(double ***old, double ***new, long n_old)
{

  // Restricts values in array old and put them in new.
  // Size of the arrays should be ok.
  // n_old = size of old grid (fine).
  //
  // old = fine
  // new = coarse
  
  int i, j, k;
  int i_new, j_new, k_new;

  int ip1, jp1, kp1;

  int n_oldm1;  // optimization
  
  //n_oldm1 = gg.n[g_old]-1;
  n_oldm1 = n_old - 1;

  // Limits are in a way that we only touch the inner
  // nodes of the upper grid.
  // -4 because we advance +2
  // Loops are for the fine grid.
#ifdef OPENMP
#pragma omp parallel private(i, j, k, ip1, jp1, kp1, i_new, j_new, k_new) shared(n_oldm1, new, old)
  {
#pragma omp for schedule(static)
#endif
  for(i=0; i<=n_oldm1; i+=2)   // Fills all the values in the grid.
    {
      for(j=0; j<=n_oldm1; j+=2)
	{
	  for(k=0; k<=n_oldm1; k+=2)
	    {
	      // periodic stuff:
	      ip1 = (i==(n_oldm1) ? 0       : (i+1));
	      jp1 = (j==(n_oldm1) ? 0       : (j+1));
	      kp1 = (k==(n_oldm1) ? 0       : (k+1));
	      
	      i_new = (int)floor(i/2.0);
	      j_new = (int)floor(j/2.0);	    
	      k_new = (int)floor(k/2.0);	    
	      new[i_new][j_new][k_new] = ( old[i  ][j  ][k  ] + 
					   old[ip1][j  ][k  ] +
					   old[i  ][jp1][k  ] +
					   old[i  ][j  ][kp1] +
					   old[ip1][jp1][k  ] + 
					   old[ip1][j  ][kp1] +
					   old[i  ][jp1][kp1] +
					   old[ip1][jp1][kp1] ) / 8.0;
	      //fprintf(stderr, "%d  %d       ", i_new, j_new);
	    }
	}
      //fprintf(stderr, "\n");
    }
#ifdef OPENMP
  }
#endif
  
}
