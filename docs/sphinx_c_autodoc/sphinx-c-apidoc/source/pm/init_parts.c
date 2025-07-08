/**
 * Here we allocate memory for the particles and read initial conditions.
 **/

#include <stdio.h>
#include <stdlib.h>

#include "../solve.h"
#include "init_parts.h"

#include "input_particles.h"

#include "../libs/allocator.h"
#include "../libs/releaser.h"

//====================
// init_parts
//====================
void init_parts(struct params *par, struct particles *part, struct units *unit)
{
  int i;

  // Allocates particles:
  //=====================
  part->pos_1 = calloc(par->nparts, sizeof(double));
  part->pos_2 = calloc(par->nparts, sizeof(double));
  part->pos_3 = calloc(par->nparts, sizeof(double));
  part->vel_1 = calloc(par->nparts, sizeof(double));
  part->vel_2 = calloc(par->nparts, sizeof(double));
  part->vel_3 = calloc(par->nparts, sizeof(double));
  
  // Read initial conditions for particles:
  //=======================================
  read_particles(par, unit, 
		 part->pos_1, part->pos_2, part->pos_3, 
		 part->vel_1, part->vel_2, part->vel_3, 
		 par->file_in_parts);

}


