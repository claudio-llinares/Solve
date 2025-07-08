#include <stdlib.h>
#include <stdio.h>

//*************************************
//  alloc_int_2
//*************************************
void free_int_2(int **a, int d1, int d2)
{
  // dealocs memory allocated by the apropiated function.
  
  int i;
  
  for(i=0; i<d1; i++)
    free(a[i]);

  free(a);
  
}

//*************************************
//  free_int_3
//*************************************
void free_int_3(int ***a, int d1, int d2, int d3)
{
  // dealocs memory allocated by the apropiated function.
  
  int i, j;
  
  for(i=0; i<d1; i++)
  {
    for(j=0; j<d2; j++)
      free(a[i][j]);
    free(a[i]);
  }
  
  free(a);
  
}

//*************************************
//  free_int_4
//*************************************
void free_int_4(int ****a, int d1, int d2, int d3, int d4)
{
  // dealocs memory allocated by the apropiated function.
  
  int i, j, k;
  
  for(i=0; i<d1; i++)
  {
    for(j=0; j<d2; j++)
    {
      for(k=0; k<d3; k++)
        free(a[i][j][k]);
      free(a[i][j]);
    }    
    free(a[i]);
  }
  
  free(a);
  
}


//*************************************
//  alloc_float_2
//*************************************
void free_float_2(float **a, int d1, int d2)
{
  // dealocs memory allocated by the apropiated function.
  
  int i;
  
  for(i=0; i<d1; i++)
    free(a[i]);
  
  free(a);
  
}

//*************************************
//  free_float_3
//*************************************
void free_float_3(float ***a, int d1, int d2, int d3)
{
  // dealocs memory allocated by the apropiated function.
  
  int i, j;
  
  for(i=0; i<d1; i++)
  {
    for(j=0; j<d2; j++)
      free(*(*(a+i)+j));
    free(a[i]);
  }
  
  free(a);
  
}

//*************************************
//  free_float_4
//*************************************
void free_float_4(float ****a, int d1, int d2, int d3, int d4)
{
  // dealocs memory allocated by the apropiated function.
  
  int i, j, k;
  
  for(i=0; i<d1; i++)
  {
    for(j=0; j<d2; j++)
    {
      for(k=0; k<d3; k++)
        free(a[i][j][k]);
      free(a[i][j]);
    }
    free(a[i]);
  }
  
  free(a);
  
}


//*************************************
//  alloc_double_2
//*************************************
void free_double_2(double **a, int d1, int d2)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0]);
  free(a);
 
}

//*************************************
//  free_double_3
//*************************************
void free_double_3(double ***a, int d1, int d2, int d3)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0]);
  free(a[0]);
  free(a);
    
}

//*************************************
//  free_double_4
//*************************************
void free_double_4(double ****a, int d1, int d2, int d3, int d4)
{
  // dealocs memory allocated by the apropiated function.
  
  int i, j, k;
  
  for(i=0; i<d1; i++)
  {
    for(j=0; j<d2; j++)
    {
      for(k=0; k<d3; k++)
        free(a[i][j][k]);
      free(a[i][j]);
    }
    free(a[i]);
  }
  
  free(a);
  
}

//*************************************
//  alloc_long_2
//*************************************
void free_long_2(long **a, int d1, int d2)
{
  // dealocs memory allocated by the apropiated function.
  
  int i;
  
  for(i=0; i<d1; i++)
    free(a[i]);
  
  free(a);
  
}

//*************************************
//  free_long_3
//*************************************
void free_long_3(long ***a, int d1, int d2, int d3)
{
  // dealocs memory allocated by the apropiated function.
  
  int i, j;
  
  for(i=0; i<d1; i++)
  {
    for(j=0; j<d2; j++)
      free(a[i][j]);
    free(a[i]);
  }
  
  free(a);
  
}

//*************************************
//  free_long_4
//*************************************
void free_long_4(long ****a, int d1, int d2, int d3, int d4)
{
  // dealocs memory allocated by the apropiated function.
  
  int i, j, k;
  
  for(i=0; i<d1; i++)
  {
    for(j=0; j<d2; j++)
    {
      for(k=0; k<d3; k++)
        free(a[i][j][k]);
      free(a[i][j]);
    }
    free(a[i]);
  }
  
  free(a);
  
}


//*************************************
//  alloc_char_2
//*************************************
void free_char_2(char **a, int d1, int d2)
{
  // dealocs memory allocated by the apropiated function.
  
  int i;
  
  for(i=0; i<d1; i++)
    free(a[i]);
  
  free(a);
  
}

//*************************************
//  free_char_3
//*************************************
void free_char_3(char ***a, int d1, int d2, int d3)
{
  // dealocs memory allocated by the apropiated function.
  
  int i, j;
  
  for(i=0; i<d1; i++)
  {
    for(j=0; j<d2; j++)
      free(a[i][j]);
    free(a[i]);
  }
  
  free(a);
  
}

//*************************************
//  free_char_4
//*************************************
void free_char_4(char ****a, int d1, int d2, int d3, int d4)
{
  // dealocs memory allocated by the apropiated function.
  
  int i, j, k;
  
  for(i=0; i<d1; i++)
  {
    for(j=0; j<d2; j++)
    {
      for(k=0; k<d3; k++)
        free(a[i][j][k]);
      free(a[i][j]);
    }
    free(a[i]);
  }
  
  free(a);
  
}
