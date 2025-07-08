/// \file read_params.c
/// \brief Reads the parameters file which is assumed to be called run.input

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../solve.h"
#include "../grid/init_grid.h"

#include "utils.h"

void init_keywords(struct params *par, char keywords[50][1000], void *values[50], char type[50][10], int *found);
void check_keyword(char keywords[50][1000], void *values[50], char type[50][10], int *found, char *key, char *val);
void delete_spaces(char *str);


//======================
// init_keywords
//======================
void init_keywords(struct params *par, char keywords[50][1000], void *values[50], char type[50][10], int *found)
{

  int i;

  strcpy(keywords[0],  "file_in_parts");
  strcpy(keywords[1],  "file_out");
  strcpy(keywords[2],  "box");
  strcpy(keywords[3],  "grid");
  strcpy(keywords[4],  "grid_p");
  strcpy(keywords[5],  "omegam");
  strcpy(keywords[6],  "omegal");
  strcpy(keywords[7],  "hsmall");
  strcpy(keywords[8],  "z_init");
  strcpy(keywords[9],  "z_end");
  strcpy(keywords[10], "nthreads");
  strcpy(keywords[11], "n_outputs");
  strcpy(keywords[12], "nsteps");
  strcpy(keywords[13], "verbose");
  strcpy(keywords[14], "static_test");
  strcpy(keywords[15], "gravity");
  strcpy(keywords[16], "mm_a0");
  strcpy(keywords[17], "mm_M");
  strcpy(keywords[18], "mm_KB");
  strcpy(keywords[19], "verbose_mm");
  strcpy(keywords[20], "verbose_mm_coarse");
  keywords[21][0] = '\0';
  
  values[0]  = &par->file_in_parts;
  values[1]  = &par->file_out;
  values[2]  = &par->box;
  values[3]  = &par->grid;
  values[4]  = &par->grid_p;
  values[5]  = &par->omegam;
  values[6]  = &par->omegal;
  values[7]  = &par->hsmall;
  values[8]  = &par->z_init;
  values[9]  = &par->z_end;
  values[10] = &par->nthreads;
  values[11] = &par->n_outputs;
  values[12] = &par->nsteps;
  values[13] = &par->verbose;
  values[14] = &par->static_test;
  values[15] = &par->gravity;
  values[16] = &par->mm_a0;
  values[17] = &par->mm_M;
  values[18] = &par->mm_KB;
  values[19] = &par->verbose_mm;
  values[20] = &par->verbose_mm_coarse;

  strcpy(type[0],  "char");
  strcpy(type[1],  "char");
  strcpy(type[2],  "double");
  strcpy(type[3],  "int");
  strcpy(type[4],  "int");
  strcpy(type[5],  "double");
  strcpy(type[6],  "double");
  strcpy(type[7],  "double");
  strcpy(type[8],  "double");
  strcpy(type[9],  "double");
  strcpy(type[10], "int");
  strcpy(type[11], "int");
  strcpy(type[12], "int");
  strcpy(type[13], "int");  
  strcpy(type[14], "int"); 
  strcpy(type[15], "int"); 
  strcpy(type[16], "double");
  strcpy(type[17], "double");
  strcpy(type[18], "double"); 
  strcpy(type[19], "int"); 
  strcpy(type[20], "int"); 

  for(i=0; i<50; i++)
    found[i] = 0;

}


//=====================
// read_params
//=====================
void read_params(struct params *par)
{

  FILE *fp;

  char line_orig[1000];    // original line in input file
  char line_orig_clean[1000];  // original line without spaces
  char line[1000];
  char *lineep;
  char *line_val;

  char *ret_fget;

  char keywords[50][1000];
  void *values[50];
  char type[50][10];
  int found[50];  // flag for founded keyword

  int i;

  // Initializes the list of keywords
  init_keywords(par, keywords, values, type, found);

  fp = fopen_check("run.input", "r");
  while(!feof(fp))
    {
      // Get a line from input file
      ret_fget = fgets(line_orig, 1000, fp);
      if(feof(fp))
	break;
      line_orig[strlen(line_orig)-1] = '\0';  // Get rid of '\n'
      strcpy(line_orig_clean, line_orig);
      delete_spaces(line_orig_clean);
            
      if(line_orig[0]!='#')   // if is not a comment
	{
	  // Look for the "="
	  strcpy(line, line_orig_clean);
	  lineep = strstr(line, "=");
	  if(lineep == NULL)
	    {
	      //fprintf(stderr, "Ignored line in input file\n");
	      continue;
	    }
	  else
	    {
	      // Cut the line to identify keyword
	      *lineep = '\0';
	      line_val = line_orig_clean;
	      // Obtain the value
	      strsep(&line_val, "=");
	      // Check if this keywor exist in the list
	      check_keyword(keywords, values, type, found, line, line_val);
	      fprintf(stderr, "%s\n", line_orig);
	    }
	}
    }
  fclose(fp);
  fprintf(stderr, "\n");

  // Check that we found all the keywords:
  i = 0;
  while( keywords[i][0] != '\0')
    {
      if(!found[i])
	{
	  fprintf(stderr, "Keyword \"%s\" is not present in input file.  Aborting.\n", keywords[i]);
	  exit(0);
	}
      i++;
    }
}

//==================================
// check_keyword
//==================================
void check_keyword(char keywords[50][1000], void *values[50], char type[50][10], int *found, char *key, char *val)
{

  int i;
  int founded;

  i = 0;
  founded = 0;

  while( keywords[i][0] != '\0')
    {
      if(!strcmp(keywords[i], key))
	{
	  if(found[i])  // we already have it
	    {
	      fprintf(stderr, "Duplicated keyword: %s. Aborting.\n", key);
	      exit(0);
	    }
	  if(!strcmp(type[i], "char"))
	    strcpy((char *)values[i], val);
	  else
	    if(!strcmp(type[i], "int"))
	      *((int *)values[i]) = atoi(val);
	    else
	      if(!strcmp(type[i], "double"))
		*((double *)values[i]) = atof(val);
	      else
		{
		  fprintf(stderr, "check_keyword: unrecognized type.  Aborting.\n");
		  exit(0);
		}
	  // Change all the flags
	  found[i] = 1;
	  founded = 1;
	  break;
	}

      i++;
    }

  if(!founded)
    {
      fprintf(stderr, "Unrecognized keyword in input file: %s. Aborting...\n", key);
      exit(0);
    }
  
}


//========================
// delete_spaces
//========================
void delete_spaces(char *str)
{

  //  Delete spaces of string str and put the new strng in the same place.

  int i, j;

  i = 0;
  j = 0;

  while(str[i]!='\0')
    {
      if(str[i]!=' ')
	{
	  str[j] = str[i];
	  j++;
	}
      i++;
    }
  
  str[j] = '\0';
  

}


