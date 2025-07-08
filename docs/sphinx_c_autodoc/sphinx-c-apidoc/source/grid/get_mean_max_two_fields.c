#include <math.h>
#include <omp.h>

#include "../solve.h"
//===================================
// get_mean_max_two_fields
//===================================
void get_mean_max_two_fields(double ***field1, double ***field2, int n, double *mean, double *maximum)
{

  int i, j, k;
  int nm1 = n-1;
  double m1, m2, maximum_local;   // temporary variables


#define new_max(x,y) ((x) >= (y) ? (x) : (y))

  m1 = 0.0;
  m2 = 0.0;
  maximum_local = -1.0e+50;
#ifdef OPENMP
#pragma omp parallel reduction(+:m1) reduction(+:m2) reduction(max:maximum_local) private(i, j, k) shared(field1, field2)
  {
    m1 = 0.0;
    m2 = 0.0;
    maximum_local = -1.0e+50;
#pragma omp for schedule(static)
#endif
    for(i=0; i<=nm1; i++)
      for(j=0; j<=nm1; j++)
	for(k=0; k<=nm1; k++)
	  {
	    
	    // Calculate norm2
	    m1 += pow2(field1[i][j][k]);
	    m2 += pow2(field2[i][j][k]);
	    
	    // Calculate maximum:
	    if(   fabs(field1[i][j][k])>maximum_local || fabs(field2[i][j][k])>maximum_local   )
#ifdef OPENMP
#pragma omp critical
#endif
	      {
		// I think we have the same "if" here to avoid synchronization problems when running with OpenMP
		if(   fabs(field1[i][j][k])>maximum_local || fabs(field2[i][j][k])>maximum_local   )
		  maximum_local = new_max( fabs(field1[i][j][k]), fabs(field2[i][j][k]) );
	      }
	  }
#ifdef OPENMP
  }
#endif
  
  *mean    = sqrt(m1) + sqrt(m2);
  *maximum = maximum_local;
  
}
