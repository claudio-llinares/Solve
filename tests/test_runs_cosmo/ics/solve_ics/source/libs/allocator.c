#include <stdlib.h>
#include <stdio.h>


//*************************************
//  alloc_int_2
//*************************************
void alloc_int_2(int ***a, int d1, int d2)
{
  // allocs int[d1][d2]
  // d1 = rows
  // d2 = columns
  // returns pointer to the array

  int i;
  
  // Allocate memory:
  *a = (int **)malloc(d1 * sizeof(int*));
  if( ((*a)[0] = (int *) malloc(d1*d2 * sizeof(int))) == NULL) 
    {
      printf("alloc_int_2:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d1; i++)
    (*a)[i] = (*a)[i-1] + d2;
  
}

//*************************************
//  alloc_int_3
//*************************************
void alloc_int_3(int ****a, int d1, int d2, int d3)
{
  // allocs int[d1][d2][d3]
  // d1 = layers?
  // d2 = rows
  // d3 = columns
  // returns pointer to the array
  
  int i, j;
  
  // Allocate memory:
  *a = (int ***)malloc(d1 * sizeof(int**));
  (*a)[0] = (int **)malloc(d1*d2 * sizeof(int*));
  (*a)[0][0] = (int *)malloc(d1*d2*d3 * sizeof(int));
  if( (*a)[0][0] == NULL) 
    {
      printf("alloc_int_3:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d2; i++)
    (*a)[0][i] = (*a)[0][i-1] + d3;

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      for(j=1; j<d2; j++)
	(*a)[i][j] = (*a)[i][j-1] + d3;
    }
  
}

//*************************************
//  alloc_int_4
//*************************************
void alloc_int_4(int *****a, int d1, int d2, int d3, int d4)
{
  // allocs int[d1][d2][d3][d4]
  // returns pointer to the array

  int i, j, k;

  // Allocate memory:
  *a = (int ****)malloc(d1 * sizeof(int***));
  (*a)[0] = (int ***)malloc(d1*d2 * sizeof(int**));
  (*a)[0][0] = (int **)malloc(d1*d2*d3 * sizeof(int*));
  (*a)[0][0][0] = (int *)malloc(d1*d2*d3*d4 * sizeof(int));
  if( (*a)[0][0][0] == NULL) 
    {
      printf("alloc_int_4:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d3; i++)
    (*a)[0][0][i] = (*a)[0][0][i-1] + d4;

  for(i=1; i<d2; i++)
    {
      (*a)[0][i] = (*a)[0][i-1] + d3;
      (*a)[0][i][0] = (*a)[0][i-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[0][i][j] = (*a)[0][i-1][j-1] + d4;
    }

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      (*a)[i][0][0] = (*a)[i-1][d2-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[i][0][j] = (*a)[i][0][j-1] + d4;

      for(j=1; j<d2; j++)
	{
	  (*a)[i][j] = (*a)[i][j-1] + d3;
	  (*a)[i][j][0] = (*a)[i][j-1][d3-1] + d4;
	  for(k=1; k<d3; k++)
	    (*a)[i][j][k] = (*a)[i][j][k-1] + d4;
	}
    }

  
}

//*************************************
//  alloc_float_2
//*************************************
void alloc_float_2(float ***a, int d1, int d2)
{
  // allocs float[d1][d2]
  // d1 = rows
  // d2 = columns
  // returns pointer to the array
  
  int i;
 
  // Allocate memory:
  *a = (float **)malloc(d1 * sizeof(float*));
  if( ((*a)[0] = (float *) malloc(d1*d2 * sizeof(float))) == NULL) 
    {
      printf("alloc_float_2:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d1; i++)
    (*a)[i] = (*a)[i-1] + d2;
  
}

//*************************************
//  alloc_float_3
//*************************************
void alloc_float_3(float ****a, int d1, int d2, int d3)
{
  // allocs float[d1][d2][d3]
  // d1 = layers?
  // d2 = rows
  // d3 = columns
  // returns pointer to the array
  
  int i, j;
  
  // Allocate memory:
  *a = (float ***)malloc(d1 * sizeof(float**));
  (*a)[0] = (float **)malloc(d1*d2 * sizeof(float*));
  (*a)[0][0] = (float *)malloc(d1*d2*d3 * sizeof(float));
  if( (*a)[0][0] == NULL) 
    {
      printf("alloc_float_3:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d2; i++)
    (*a)[0][i] = (*a)[0][i-1] + d3;

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      for(j=1; j<d2; j++)
	(*a)[i][j] = (*a)[i][j-1] + d3;
    }

}


//*************************************
//  alloc_float_4
//*************************************
void alloc_float_4(float *****a, int d1, int d2, int d3, int d4)
{
  // allocs double[d1][d2][d3][d4]
  // returns pointer to the array
  
  int i, j, k;
  
  // Allocate memory:
  *a = (float ****)malloc(d1 * sizeof(float***));
  (*a)[0] = (float ***)malloc(d1*d2 * sizeof(float**));
  (*a)[0][0] = (float **)malloc(d1*d2*d3 * sizeof(float*));
  (*a)[0][0][0] = (float *)malloc(d1*d2*d3*d4 * sizeof(float));
  if( (*a)[0][0][0] == NULL) 
    {
      printf("alloc_float_4:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d3; i++)
    (*a)[0][0][i] = (*a)[0][0][i-1] + d4;

  for(i=1; i<d2; i++)
    {
      (*a)[0][i] = (*a)[0][i-1] + d3;
      (*a)[0][i][0] = (*a)[0][i-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[0][i][j] = (*a)[0][i-1][j-1] + d4;
    }

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      (*a)[i][0][0] = (*a)[i-1][d2-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[i][0][j] = (*a)[i][0][j-1] + d4;

      for(j=1; j<d2; j++)
	{
	  (*a)[i][j] = (*a)[i][j-1] + d3;
	  (*a)[i][j][0] = (*a)[i][j-1][d3-1] + d4;
	  for(k=1; k<d3; k++)
	    (*a)[i][j][k] = (*a)[i][j][k-1] + d4;
	}
    }

}

//*************************************
//  alloc_double_2
//*************************************
void alloc_double_2(double ***a, int d1, int d2)
{
  // allocs double[d1][d2]
  // d1 = rows
  // d2 = columns
  // returns pointer to the array
  
  int i;
  int j;
  
  // Allocate memory:
  *a = (double **)malloc(d1 * sizeof(double*));
  if( ((*a)[0] = (double *)malloc(d1*d2 * sizeof(double))) == NULL) 
    {
      printf("alloc_double_2:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d1; i++)
    (*a)[i] = (*a)[i-1] + d2;
    
  // Initialize:
  for(i=0; i<d1; i++)
    for(j=0; j<d2; j++)
      (*a)[i][j] = 0.0;
  
}

//*************************************
//  alloc_double_3
//*************************************
void alloc_double_3(double ****a, int d1, int d2, int d3)
{
  // allocs double[d1][d2][d3]
  // d1 = layers?
  // d2 = rows
  // d3 = columns
  // returns pointer to the array
  
  int i, j;
  int k;
  
  // Allocate memory:
  *a = (double ***)malloc(d1 * sizeof(double**));
  (*a)[0] = (double **)malloc(d1*d2 * sizeof(double*));
  (*a)[0][0] = (double *)malloc(d1*d2*d3 * sizeof(double));
  if( (*a)[0][0] == NULL) 
    {
      printf("alloc_double_3:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d2; i++)
    (*a)[0][i] = (*a)[0][i-1] + d3;

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      for(j=1; j<d2; j++)
	(*a)[i][j] = (*a)[i][j-1] + d3;
    }
  
  // Initialize:
  for(i=0; i<d1; i++)
    for(j=0; j<d2; j++)
      for(k=0; k<d3; k++)
      (*a)[i][j][k] = 0.0;
}

//*************************************
//  alloc_double_4
//*************************************
void alloc_double_4(double *****a, int d1, int d2, int d3, int d4)
{
  // allocs double[d1][d2][d3][d4]
  // returns pointer to the array
  
  int i, j, k;
  int l;
  
  // Allocate memory:
  *a = (double ****)malloc(d1 * sizeof(double***));
  (*a)[0] = (double ***)malloc(d1*d2 * sizeof(double**));
  (*a)[0][0] = (double **)malloc(d1*d2*d3 * sizeof(double*));
  (*a)[0][0][0] = (double *)malloc(d1*d2*d3*d4 * sizeof(double));
  if( (*a)[0][0][0] == NULL) 
    {
      printf("alloc_double_4:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d3; i++)
    (*a)[0][0][i] = (*a)[0][0][i-1] + d4;

  for(i=1; i<d2; i++)
    {
      (*a)[0][i] = (*a)[0][i-1] + d3;
      (*a)[0][i][0] = (*a)[0][i-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[0][i][j] = (*a)[0][i-1][j-1] + d4;
    }

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      (*a)[i][0][0] = (*a)[i-1][d2-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[i][0][j] = (*a)[i][0][j-1] + d4;

      for(j=1; j<d2; j++)
	{
	  (*a)[i][j] = (*a)[i][j-1] + d3;
	  (*a)[i][j][0] = (*a)[i][j-1][d3-1] + d4;
	  for(k=1; k<d3; k++)
	    (*a)[i][j][k] = (*a)[i][j][k-1] + d4;
	}
    }
  
  // Initialize:
  for(i=0; i<d1; i++)
    for(j=0; j<d2; j++)
      for(k=0; k<d3; k++)
	for(l=0; l<d4; l++)
	  (*a)[i][j][k][l] = 0.0;
}

//*************************************
//  alloc_double_5
//*************************************
void alloc_double_5(double ******a, int d1, int d2, int d3, int d4, int d5)
{
  // allocs double[d1][d2][d3][d4]
  // returns pointer to the array
  
  int i, j, k;
  int l, m;
  
  // Allocate memory:
  *a = (double *****)malloc(d1 * sizeof(double***));
  (*a)[0] = (double ****)malloc(d1*d2 * sizeof(double**));
  (*a)[0][0] = (double ***)malloc(d1*d2*d3 * sizeof(double*));
  (*a)[0][0][0] = (double **)malloc(d1*d2*d3*d4 * sizeof(double));
  (*a)[0][0][0][0] = (double *)malloc(d1*d2*d3*d4*d5 * sizeof(double));
  if( (*a)[0][0][0][0] == NULL) 
    {
      printf("alloc_double_5:  Problems allocating memory\n");
      exit(0);
    }

  fprintf(stderr, "alloc_double_5:  this routine is incomplete.  ABorting.\n");
  exit(0);
  // Arrange pointers:
  for(i=1; i<d3; i++)
    (*a)[0][0][i] = (*a)[0][0][i-1] + d4;

  for(i=1; i<d2; i++)
    {
      (*a)[0][i] = (*a)[0][i-1] + d3;
      (*a)[0][i][0] = (*a)[0][i-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[0][i][j] = (*a)[0][i-1][j-1] + d4;
    }

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      (*a)[i][0][0] = (*a)[i-1][d2-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[i][0][j] = (*a)[i][0][j-1] + d4;

      for(j=1; j<d2; j++)
	{
	  (*a)[i][j] = (*a)[i][j-1] + d3;
	  (*a)[i][j][0] = (*a)[i][j-1][d3-1] + d4;
	  for(k=1; k<d3; k++)
	    (*a)[i][j][k] = (*a)[i][j][k-1] + d4;
	}
    }
  
  // Initialize:
  for(i=0; i<d1; i++)
    for(j=0; j<d2; j++)
      for(k=0; k<d3; k++)
	for(l=0; l<d3; l++)
	  for(m=0; m<d3; m++)
	  (*a)[i][j][k][l][m] = 0.0;
}

//*************************************
//  alloc_long_2
//*************************************
void alloc_long_2(long ***a, int d1, int d2)
{
  // allocs long[d1][d2]
  // d1 = rows
  // d2 = columns
  // returns pointer to the array
    
  int i;
    
  // Allocate memory:
  *a = (long **)malloc(d1 * sizeof(long*));
  if( ((*a)[0] = (long *)malloc(d1*d2 * sizeof(long))) == NULL) 
    {
      printf("alloc_long_2:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d1; i++)
    (*a)[i] = (*a)[i-1] + d2;
    
}
  
//*************************************
//  alloc_long_3
//*************************************
void alloc_long_3(long ****a, int d1, int d2, int d3)
{
  // allocs long[d1][d2][d3]
  // d1 = layers?
  // d2 = rows
  // d3 = columns
  // returns pointer to the array
    
  int i, j;
    
  // Allocate memory:
  *a = (long ***)malloc(d1 * sizeof(long**));
  (*a)[0] = (long **)malloc(d1*d2 * sizeof(long*));
  (*a)[0][0] = (long *)malloc(d1*d2*d3 * sizeof(long));
  if( (*a)[0][0] == NULL) 
    {
      printf("alloc_long_3:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d2; i++)
    (*a)[0][i] = (*a)[0][i-1] + d3;

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      for(j=1; j<d2; j++)
	(*a)[i][j] = (*a)[i][j-1] + d3;
    }
  
    
}
  
//*************************************
//  alloc_long_4
//*************************************
void alloc_long_4(long *****a, int d1, int d2, int d3, int d4)
{
  // allocs long[d1][d2][d3][d4]
  // returns pointer to the array
    
  int i, j, k;
    
  // Allocate memory:
  *a = (long ****)malloc(d1 * sizeof(long***));
  (*a)[0] = (long ***)malloc(d1*d2 * sizeof(long**));
  (*a)[0][0] = (long **)malloc(d1*d2*d3 * sizeof(long*));
  (*a)[0][0][0] = (long *)malloc(d1*d2*d3*d4 * sizeof(long));
  if( (*a)[0][0][0] == NULL) 
    {
      printf("alloc_long_4:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d3; i++)
    (*a)[0][0][i] = (*a)[0][0][i-1] + d4;

  for(i=1; i<d2; i++)
    {
      (*a)[0][i] = (*a)[0][i-1] + d3;
      (*a)[0][i][0] = (*a)[0][i-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[0][i][j] = (*a)[0][i-1][j-1] + d4;
    }

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      (*a)[i][0][0] = (*a)[i-1][d2-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[i][0][j] = (*a)[i][0][j-1] + d4;

      for(j=1; j<d2; j++)
	{
	  (*a)[i][j] = (*a)[i][j-1] + d3;
	  (*a)[i][j][0] = (*a)[i][j-1][d3-1] + d4;
	  for(k=1; k<d3; k++)
	    (*a)[i][j][k] = (*a)[i][j][k-1] + d4;
	}
    }  
  
}  

//*************************************
//  alloc_char_2
//*************************************
void alloc_char_2(char ***a, int d1, int d2)
{
  // allocs char[d1][d2]
  // d1 = rows
  // d2 = columns
  // returns pointer to the array
  
  int i;
  
  // Allocate memory:
  *a = (char **)malloc(d1 * sizeof(char*));
  if( ((*a)[0] = (char *)malloc(d1*d2 * sizeof(char))) == NULL) 
    {
      printf("alloc_char_2:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d1; i++)
    (*a)[i] = (*a)[i-1] + d2;
  
}

//*************************************
//  alloc_char_3
//*************************************
void alloc_char_3(char ****a, int d1, int d2, int d3)
{
  // allocs char[d1][d2][d3]
  // d1 = layers?
  // d2 = rows
  // d3 = columns
  // returns pointer to the array
  
  int i, j;
    
  // Allocate memory:
  *a = (char ***)malloc(d1 * sizeof(char**));
  (*a)[0] = (char **)malloc(d1*d2 * sizeof(char*));
  (*a)[0][0] = (char *)malloc(d1*d2*d3 * sizeof(char));
  if( (*a)[0][0] == NULL) 
    {
      printf("alloc_char_3:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d2; i++)
    (*a)[0][i] = (*a)[0][i-1] + d3;

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      for(j=1; j<d2; j++)
	(*a)[i][j] = (*a)[i][j-1] + d3;
    }
  
}

//*************************************
//  alloc_char_4
//*************************************
void alloc_char_4(char *****a, int d1, int d2, int d3, int d4)
{
  // allocs char[d1][d2][d3][d4]
  // returns pointer to the array
  
  int i, j, k;
  
  // Allocate memory:
  *a = (char ****)malloc(d1 * sizeof(char***));
  (*a)[0] = (char ***)malloc(d1*d2 * sizeof(char**));
  (*a)[0][0] = (char **)malloc(d1*d2*d3 * sizeof(char*));
  (*a)[0][0][0] = (char *)malloc(d1*d2*d3*d4 * sizeof(char));
  if( (*a)[0][0][0] == NULL) 
    {
      printf("alloc_char_4:  Problems allocating memory\n");
      exit(0);
    }
  // Arrange pointers:
  for(i=1; i<d3; i++)
    (*a)[0][0][i] = (*a)[0][0][i-1] + d4;

  for(i=1; i<d2; i++)
    {
      (*a)[0][i] = (*a)[0][i-1] + d3;
      (*a)[0][i][0] = (*a)[0][i-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[0][i][j] = (*a)[0][i-1][j-1] + d4;
    }

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      (*a)[i][0][0] = (*a)[i-1][d2-1][d3-1] + d4;
      for(j=1; j<d3; j++)
	(*a)[i][0][j] = (*a)[i][0][j-1] + d4;

      for(j=1; j<d2; j++)
	{
	  (*a)[i][j] = (*a)[i][j-1] + d3;
	  (*a)[i][j][0] = (*a)[i][j-1][d3-1] + d4;
	  for(k=1; k<d3; k++)
	    (*a)[i][j][k] = (*a)[i][j][k-1] + d4;
	}
    }
  
  
}  



