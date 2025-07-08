#include <stdlib.h>
#include <stdio.h>


//*************************************
//  alloc_int_2
//*************************************
int **alloc_int_2(int ***a, int d1, int d2)
{
    // allocs int[d1][d2]
    // d1 = rows
    // d2 = columns
    // returns pointer to the array

  //int **a;
  int i;
  
   // Allocate memory:
  *a = (int **)calloc(d1, sizeof(int*));
  if( ((*a)[0] = (int *) calloc(d1 * d2, sizeof(int))) == NULL) 
    {
      printf("alloc_int_2:  Problems allocating memory\n");
      exit(0);
    }
  // Accomodate pointers:
  for(i=1; i<d1; i++)
    (*a)[i] = (*a)[i-1] + d2;
  
  return *a;

}

//*************************************
//  alloc_int_3
//*************************************
int ***alloc_int_3(int d1, int d2, int d3)
{
    // allocs int[d1][d2][d3]
    // d1 = layers?
    // d2 = rows
    // d3 = columns
    // returns pointer to the array

    int ***a;
    int i, j;


    //allocate memory for layers
  if( (a = (int***)malloc(d1*sizeof(int**))) == NULL )
  {
    printf("alloc_int_3:  Problems allocating memory\n");
    exit(0);
  }
  
    //for each layer allocate memory for rows 
  for(i=0; i<d1; i++)
  {
    if( (*(a+i) = (int**)malloc(d2*sizeof(int*))) == NULL )
    {
      printf("alloc_int_3:  Problems allocating memory\n");
      exit(0);
    }
        //for each row allocate memory for columns
    for(j=0; j<d2; j++)
    {
      if( (*(*(a+i)+j) = (int*)malloc(d3*sizeof(int))) == NULL )
      {
        printf("alloc_int_3:  Problems allocating memory\n");
        exit(0);
      }
    }
  }
  
  return a;

}

//*************************************
//  alloc_int_4
//*************************************
int ****alloc_int_4(int d1, int d2, int d3, int d4)
{
    // allocs int[d1][d2][d3][d4]
    // returns pointer to the array

    int ****a;
    int i, j, k;


    //allocate memory for d1
  if( (a = (int****)malloc(d1*sizeof(int***))) == NULL )
  {
    printf("alloc_int_4:  Problems allocating memory\n");
    exit(0);
  }
    
    //for each layer allocate memory d2 
  for(i=0; i<d1; i++)
  {
    if( (*(a+i) = (int***)malloc(d2*sizeof(int**))) == NULL )
    {
      printf("alloc_int_4:  Problems allocating memory\n");
      exit(0);
    }
        //for each row allocate memory d3
    for(j=0; j<d2; j++)
    {
      if( (*(*(a+i)+j) = (int**)malloc(d3*sizeof(int*))) == NULL )
      {
        printf("alloc_int_4:  Problems allocating memory\n");
        exit(0);
      }
      for(k=0; k<d3; k++)
      {
        if( (*(*(*(a+i)+j)+k) = (int*)malloc(d4*sizeof(int))) == NULL )
        {
          printf("alloc_int_4:  Problems allocating memory\n");
          exit(0);
        }
      }
    }
  }
  
    return a;

}

//*************************************
//  alloc_float_2
//*************************************
float **alloc_float_2(int d1, int d2)
{
    // allocs float[d1][d2]
    // d1 = rows
    // d2 = columns
    // returns pointer to the array
  
  float **a;
  int i;
 
  // Allocate memory:
  a = (float **)malloc(d1 * sizeof(float*));
  if( (a[0] = (float *) malloc(d1 * d2 * sizeof(float))) == NULL) 
    {
      printf("alloc_float_2:  Problems allocating memory\n");
      exit(0);
    }
  // Accomodate pointers:
  for(i=1; i<d1; i++)
    a[i] = a[i-1] + d2;
  
  return a;
  
}

//*************************************
//  alloc_float_3
//*************************************
float ***alloc_float_3(int d1, int d2, int d3)
{
    // allocs float[d1][d2][d3]
    // d1 = layers?
    // d2 = rows
    // d3 = columns
    // returns pointer to the array
  
  float ***a;
  int i, j;
  
  
    //allocate memory for layers
  if( (a = (float***)malloc(d1*sizeof(float**))) == NULL )
  {
  printf("alloc_float_3:  Problems allocating memory\n");
    exit(0);
  }
  
    //for each layer allocate memory for rows 
  for(i=0; i<d1; i++)
  {
    if( (*(a+i) = (float**)malloc(d2*sizeof(float*))) == NULL )
    {
    printf("alloc_float_3:  Problems allocating memory\n");
      exit(0);
    }
        //for each row allocate memory for columns
    for(j=0; j<d2; j++)
    {
      if( (*(*(a+i)+j) = (float*)malloc(d3*sizeof(float))) == NULL )
      {
      printf("alloc_float_3:  Problems allocating memory\n");
        exit(0);
      }
    }
  }
  
  return a;
  
}

//*************************************
//  alloc_float_4
//*************************************
float ****alloc_float_4(int d1, int d2, int d3, int d4)
{
    // allocs float[d1][d2][d3][d4]
    // returns pointer to the array
  
  float ****a;
  int i, j, k;
  
  
    //allocate memory for d1
  if( (a = (float****)malloc(d1*sizeof(float***))) == NULL )
  {
  printf("alloc_float_4:  Problems allocating memory\n");
    exit(0);
  }
  
    //for each layer allocate memory d2 
  for(i=0; i<d1; i++)
  {
    if( (*(a+i) = (float***)malloc(d2*sizeof(float**))) == NULL )
    {
    printf("alloc_float_4:  Problems allocating memory\n");
      exit(0);
    }
        //for each row allocate memory d3
    for(j=0; j<d2; j++)
    {
      if( (*(*(a+i)+j) = (float**)malloc(d3*sizeof(float*))) == NULL )
      {
      printf("alloc_float_4:  Problems allocating memory\n");
        exit(0);
      }
      for(k=0; k<d3; k++)
      {
        if( (*(*(*(a+i)+j)+k) = (float*)malloc(d4*sizeof(float))) == NULL )
        {
        printf("alloc_float_4:  Problems allocating memory\n");
          exit(0);
        }
      }
    }
  }
  
  return a;
  
}

//*************************************
//  alloc_double_2
//*************************************
double **alloc_double_2(double ***a, int d1, int d2)
{
    // allocs double[d1][d2]
    // d1 = rows
    // d2 = columns
    // returns pointer to the array
  
  //double **a;
  int i;
  
   // Allocate memory:
  //a = (double **)malloc(d1 * sizeof(double*));
  *a = (double **)calloc(d1, sizeof(double*));
  //if( (a[0] = (double *) malloc(d1 * d2 * sizeof(double))) == NULL) 
  if( ((*a)[0] = (double *) calloc(d1 * d2, sizeof(double))) == NULL) 
    {
      printf("alloc_double_2:  Problems allocating memory\n");
      exit(0);
    }
  // Accomodate pointers:
  for(i=1; i<d1; i++)
    (*a)[i] = (*a)[i-1] + d2;
  
  return *a;
  
}

//*************************************
//  alloc_double_3
//*************************************
double ***alloc_double_3(double ****a, int d1, int d2, int d3)
{
    // allocs double[d1][d2][d3]
    // d1 = layers?
    // d2 = rows
    // d3 = columns
    // returns pointer to the array
  
  //double ***a;
  int i, j;
  
  // Allocate memory:
  *a = (double ***)malloc(d1 * sizeof(double**));
  (*a)[0] = (double **)malloc(d1*d2 * sizeof(double*));
  (*a)[0][0] = (double *)malloc(d1*d2*d3 * sizeof(double));
  if( (*a)[0][0] == NULL) 
    {
      printf("alloc_double_3:  Problems allocating memory\n");
      exit(0);
    }
  // Accomodate pointers:
  for(i=1; i<d2; i++)
    (*a)[0][i] = (*a)[0][i-1] + d3;

  for(i=1; i<d1; i++)
    {
      (*a)[i] = (*a)[i-1] + d2;
      (*a)[i][0] = (*a)[i-1][d2-1] + d3;
      for(j=1; j<d2; j++)
	(*a)[i][j] = (*a)[i][j-1] + d3;
    }

  return *a;
  
}

//*************************************
//  alloc_double_4
//*************************************
double ****alloc_double_4(double *****a, int d1, int d2, int d3, int d4)
{
    // allocs double[d1][d2][d3][d4]
    // returns pointer to the array
  
  //double ****a;
  int i, j, k;
  
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
  // Accomodate pointers:
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

  return *a;
  
}

//*************************************
//  alloc_long_2
//*************************************
  long **alloc_long_2(int d1, int d2)
  {
    // allocs long[d1][d2]
    // d1 = rows
    // d2 = columns
    // returns pointer to the array
    
    long **a;
    int i;
    
    //allocate memory for rows
    if( (a = (long**)malloc(d1*sizeof(long*))) == NULL )
    {
    printf("alloc_long_2:  Problems allocating memory\n");
      exit(0);
    }
    //for each row allocate memory for columns
    for(i=0; i<d1; i++)
    {
      if( (*(a+i) = (long*)malloc(d2*sizeof(long))) == NULL )
      {
      printf("alloc_long_2:  Problems allocating memory\n");
        exit(0);
      }
    }
    
    return a;
    
  }
  
//*************************************
//  alloc_long_3
//*************************************
  long ***alloc_long_3(int d1, int d2, int d3)
  {
    // allocs long[d1][d2][d3]
    // d1 = layers?
    // d2 = rows
    // d3 = columns
    // returns pointer to the array
    
    long ***a;
    int i, j;
    
    
    //allocate memory for layers
    if( (a = (long***)malloc(d1*sizeof(long**))) == NULL )
    {
    printf("alloc_long_3:  Problems allocating memory\n");
      exit(0);
    }
    
    //for each layer allocate memory for rows 
    for(i=0; i<d1; i++)
    {
      if( (*(a+i) = (long**)malloc(d2*sizeof(long*))) == NULL )
      {
      printf("alloc_long_3:  Problems allocating memory\n");
        exit(0);
      }
        //for each row allocate memory for columns
      for(j=0; j<d2; j++)
      {
        if( (*(*(a+i)+j) = (long*)malloc(d3*sizeof(long))) == NULL )
        {
        printf("alloc_long_3:  Problems allocating memory\n");
          exit(0);
        }
      }
    }
    
    return a;
    
  }
  
//*************************************
//  alloc_long_4
//*************************************
  long ****alloc_long_4(int d1, int d2, int d3, int d4)
  {
    // allocs long[d1][d2][d3][d4]
    // returns pointer to the array
    
    long ****a;
    int i, j, k;
    
    
    //allocate memory for d1
    if( (a = (long****)malloc(d1*sizeof(long***))) == NULL )
    {
    printf("alloc_long_4:  Problems allocating memory\n");
      exit(0);
    }
    
    //for each layer allocate memory d2 
    for(i=0; i<d1; i++)
    {
      if( (*(a+i) = (long***)malloc(d2*sizeof(long**))) == NULL )
      {
      printf("alloc_long_4:  Problems allocating memory\n");
        exit(0);
      }
        //for each row allocate memory d3
      for(j=0; j<d2; j++)
      {
        if( (*(*(a+i)+j) = (long**)malloc(d3*sizeof(long*))) == NULL )
        {
        printf("alloc_long_4:  Problems allocating memory\n");
          exit(0);
        }
        for(k=0; k<d3; k++)
        {
          if( (*(*(*(a+i)+j)+k) = (long*)malloc(d4*sizeof(long))) == NULL )
          {
          printf("alloc_long_4:  Problems allocating memory\n");
            exit(0);
          }
        }
      }
    }
    
    return a;
    
  }  

//*************************************
//  alloc_char_2
//*************************************
char **alloc_char_2(int d1, int d2)
{
    // allocs char[d1][d2]
    // d1 = rows
    // d2 = columns
    // returns pointer to the array
  
  char **a;
  int i;
  
    //allocate memory for rows
  if( (a = (char**)malloc(d1*sizeof(char*))) == NULL )
  {
  printf("alloc_char_2:  Problems allocating memory\n");
    exit(0);
  }
    //for each row allocate memory for columns
  for(i=0; i<d1; i++)
  {
    if( (*(a+i) = (char*)malloc(d2*sizeof(char))) == NULL )
    {
    printf("alloc_char_2:  Problems allocating memory\n");
      exit(0);
    }
  }
  
  return a;
  
}

//*************************************
//  alloc_char_3
//*************************************
char ***alloc_char_3(int d1, int d2, int d3)
{
    // allocs char[d1][d2][d3]
    // d1 = layers?
    // d2 = rows
    // d3 = columns
    // returns pointer to the array
  
  char ***a;
  int i, j;
  
  
    //allocate memory for layers
  if( (a = (char***)malloc(d1*sizeof(char**))) == NULL )
  {
  printf("alloc_char_3:  Problems allocating memory\n");
    exit(0);
  }
  
    //for each layer allocate memory for rows 
  for(i=0; i<d1; i++)
  {
    if( (*(a+i) = (char**)malloc(d2*sizeof(char*))) == NULL )
    {
    printf("alloc_char_3:  Problems allocating memory\n");
      exit(0);
    }
        //for each row allocate memory for columns
    for(j=0; j<d2; j++)
    {
      if( (*(*(a+i)+j) = (char*)malloc(d3*sizeof(char))) == NULL )
      {
      printf("alloc_char_3:  Problems allocating memory\n");
        exit(0);
      }
    }
  }
  
  return a;
  
}

//*************************************
//  alloc_char_4
//*************************************
char ****alloc_char_4(int d1, int d2, int d3, int d4)
{
    // allocs char[d1][d2][d3][d4]
    // returns pointer to the array
  
  char ****a;
  int i, j, k;
  
  
    //allocate memory for d1
  if( (a = (char****)malloc(d1*sizeof(char***))) == NULL )
  {
  printf("alloc_char_4:  Problems allocating memory\n");
    exit(0);
  }
  
    //for each layer allocate memory d2 
  for(i=0; i<d1; i++)
  {
    if( (*(a+i) = (char***)malloc(d2*sizeof(char**))) == NULL )
    {
    printf("alloc_char_4:  Problems allocating memory\n");
      exit(0);
    }
        //for each row allocate memory d3
    for(j=0; j<d2; j++)
    {
      if( (*(*(a+i)+j) = (char**)malloc(d3*sizeof(char*))) == NULL )
      {
      printf("alloc_char_4:  Problems allocating memory\n");
        exit(0);
      }
      for(k=0; k<d3; k++)
      {
        if( (*(*(*(a+i)+j)+k) = (char*)malloc(d4*sizeof(char))) == NULL )
        {
        printf("alloc_char_4:  Problems allocating memory\n");
          exit(0);
        }
      }
    }
  }
  
  return a;
  
}  



