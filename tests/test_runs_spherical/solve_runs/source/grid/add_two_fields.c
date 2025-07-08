#include <math.h>
#include <omp.h>

#include "../solve.h"

//===================================
// add_two_fields
//===================================
void add_two_fields(double ***field1, double ***field2, int n, double ***result)
{
  
  int i, j, k;
  int nm1 = n-1;
  
#ifdef OPENMP
#pragma omp parallel private(i, j, k) shared(field1, field2, result)
  {
#pragma omp for schedule(static)
#endif
    for(i=0; i<=nm1; i++)
      for(j=0; j<=nm1; j++)
	for(k=0; k<=nm1; k++)
	  result[i][j][k] = field1[i][j][k] + field2[i][j][k];
#ifdef OPENMP
  }
#endif
  
  return ;
  
}
