#include <stdlib.h>
#include <stdio.h>

//*************************************
//  alloc_int_2
//*************************************
void free_int_2(int **a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0]);
  free(a);

  
}

//*************************************
//  free_int_3
//*************************************
void free_int_3(int ***a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0]);
  free(a[0]);
  free(a);

}

//*************************************
//  free_int_4
//*************************************
void free_int_4(int ****a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0][0]);
  free(a[0][0]);
  free(a[0]);
  free(a);
  
}


//*************************************
//  alloc_float_2
//*************************************
void free_float_2(float **a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0]);
  free(a);
  
}

//*************************************
//  free_float_3
//*************************************
void free_float_3(float ***a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0]);
  free(a[0]);
  free(a);
  
}

//*************************************
//  free_float_4
//*************************************
void free_float_4(float ****a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0][0]);
  free(a[0][0]);
  free(a[0]);
  free(a);
  
}

//*************************************
//  alloc_double_2
//*************************************
void free_double_2(double **a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0]);
  free(a);
 
}

//*************************************
//  free_double_3
//*************************************
void free_double_3(double ***a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0]);
  free(a[0]);
  free(a);
    
}

//*************************************
//  free_double_4
//*************************************
void free_double_4(double ****a)
{
  // dealocs memory allocated by the apropiated function.

  free(a[0][0][0]);
  free(a[0][0]);
  free(a[0]);
  free(a);
    
}

//*************************************
//  alloc_long_2
//*************************************
void free_long_2(long **a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0]);
  free(a);

  
}

//*************************************
//  free_long_3
//*************************************
void free_long_3(long ***a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0]);
  free(a[0]);
  free(a);

}

//*************************************
//  free_long_4
//*************************************
void free_long_4(long ****a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0][0]);
  free(a[0][0]);
  free(a[0]);
  free(a);
  
}

//*************************************
//  alloc_char_2
//*************************************
void free_char_2(char **a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0]);
  free(a);

  
}

//*************************************
//  free_char_3
//*************************************
void free_char_3(char ***a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0]);
  free(a[0]);
  free(a);

}

//*************************************
//  free_char_4
//*************************************
void free_char_4(char ****a)
{
  // dealocs memory allocated by the apropiated function.
  
  free(a[0][0][0]);
  free(a[0][0]);
  free(a[0]);
  free(a);
  
}

