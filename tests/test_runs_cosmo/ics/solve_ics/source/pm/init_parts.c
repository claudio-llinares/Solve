#include <stdio.h>
#include <stdlib.h>

#include "../solve.h"
#include "init_parts.h"

#include "../libs/allocator.h"
#include "../libs/releaser.h"

#include "kernels.h"

//====================
// init_parts
//====================
void init_parts(struct params *par, struct particles *part, struct units *unit)
{


  // Allocates particles:
  //=====================
  part->pos_1 = calloc(par->nparts, sizeof(double));
  part->pos_2 = calloc(par->nparts, sizeof(double));
  part->pos_3 = calloc(par->nparts, sizeof(double));
  part->vel_1 = calloc(par->nparts, sizeof(double));
  part->vel_2 = calloc(par->nparts, sizeof(double));
  part->vel_3 = calloc(par->nparts, sizeof(double));
  
  // Calculates particle mass:
  //==========================
  fprintf(stderr, "Calculating particle mass from mean density.\n");
  // My version
  par->mass = (par->omegam * unit->rho_crit * pow3(par->box)) / pow3(par->grid_p); // / 1.0e+10;
  fprintf(stderr, "Particle's mass (mine) = %lg Msun = %lg Msun/h\n" , 
	  par->mass, par->mass * par->hsmall);
  // grafics
  par->mass   = 2.776e7*pow2(100.0*par->hsmall) * par->omegam *pow3(par->box/par->grid_p);
  fprintf(stderr, "Particle's mass (grafic2) = %lg Msun = %lg Msun/h\n" , 
	  par->mass, par->mass * par->hsmall);

  // Initialize pointer to kernels:
  //===============================
  part->kernel = &cic_kernel;
  
}


