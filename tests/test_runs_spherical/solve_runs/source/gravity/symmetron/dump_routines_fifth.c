#include <stdio.h>
#include <math.h>

#include "../../solve.h"

#include "../commons_multigrid/commons_multigrid.h"
#include "solve_fifth_multi.h"

//======================================
// write_1d
//======================================
void write_1d(struct params *par, struct grids_gravity *gg, int n, int grid, double a)
{

  // For debuging.

  int i, j, k;
  FILE *fp;
  char name[100];
  static int l = 0;
  float out[3];

  FILE *fp_hola;

  return;

  fp_hola = fopen("hola.txt", "a");

  // Phi:
  //=====
  sprintf(name, "pot_1d.%05d.%05d.%03ld.dat", l, n, gg->n[grid]);
  fprintf(fp_hola, "VAPORphi raw2vdf -varname phi phi_%03ld.vdf $prefix/%s\n", gg->n[grid], name);
  fp = fopen(name, "w");
  for(i=0; i<=gg->n[grid]-1; i++)
    for(j=0; j<=gg->n[grid]-1; j++)
      for(k=0; k<=gg->n[grid]-1; k++)
	{
	  out[0] = (float)gg->fi[grid][i][j][k];
	  fwrite(out, sizeof(float), 1, fp);
	}
  fclose(fp);

  // Rho:
  //=====
  sprintf(name, "rho_1d.%05d.%05d.%03ld.dat", l, n, gg->n[grid]);
  fprintf(fp_hola, "VAPORrho raw2vdf -varname rho phi_%03ld.vdf $prefix/%s\n", gg->n[grid], name);
  fp = fopen(name, "w");
  for(i=0; i<=gg->n[grid]-1; i++)
    for(j=0; j<=gg->n[grid]-1; j++)
      for(k=0; k<=gg->n[grid]-1; k++)
	{
	  out[0] = (float)gg->rho[grid][i][j][k];
	  fwrite(out, sizeof(float), 1, fp);
	}
  fclose(fp);

  // Res:
  //=====
  sprintf(name, "res_1d.%05d.%05d.%03ld.dat", l, n, gg->n[grid]);
  fprintf(fp_hola, "VAPORres raw2vdf -varname res phi_%03ld.vdf $prefix/%s\n", gg->n[grid], name);
  fp = fopen(name, "w");
  for(i=0; i<=gg->n[grid]-1; i++)
    for(j=0; j<=gg->n[grid]-1; j++)
      for(k=0; k<=gg->n[grid]-1; k++)
	{
	  out[0] = (float)gg->res[grid][i][j][k];
	  fwrite(out, sizeof(float), 1, fp);
	}
  fclose(fp);

  fclose(fp_hola);

  l++;
    
}

//======================================
// write_pot
//======================================
void write_pot_fifth(struct params *par, struct grids_gravity *gg, int n, int grid, double a)
{

  // For debuging.

  int i, j, k;
  FILE *fp;
  char name[100];
  static int l = 0;
  double rho_crit, rho_mean;

  return;

  //write_1d(par, gg, n, grid, a);
  
  //if(gg->n[grid]!=64)
  l++;
  //fprintf(stderr, "l = %d\n", l);
  //if( l!=466 && l!=2 && l!=450 )
  //  return;
  
  /* if(!(gg->n[grid]==64)) */
  /*   { */
  /*     l++; */
  /*     return; */
  /*   } */
  /* if(!(n==4 || n==1))  */
  /*   { */
  /*     l++; */
  /*     return; */
  /*   } */
  /* if( !(  */
  /* 	(n==4 && (l==5 || l==81 || l==201 || l==557 || l== 773 || l==1517 || l==2057 || l==2517 || l>1000)) ||  */
  /* 	(n==1 && l==1)  */
  /* 	 )) */
  /*   { */
  /*     l++; */
  /*     return; */
  /*   } */

  /* fprintf(stderr, "write_pot_fifth: WARNING:  Writting only specific parts of specific files...\n"); */

  rho_crit = 2.776e+7 * pow2(100.0*par->hsmall);
  rho_mean = par->omegam * rho_crit;
  
  sprintf(name, "pot.%05d.%05d.%03ld.dat", l-1, n, gg->n[grid]);
  fprintf(stderr, "Writing %s...\n", name);
  fp = fopen(name, "w");
  
  int ip1, im1, jp1, jm1, kp1, km1;
  double nb_sum;
  int nm1 = gg->n[grid]-1;
  double norm_rho;

  double mean_nabla = 0.0;
  double mean_source = 0.0;

  for(i=0; i<=gg->n[grid]-1; i++)
    {
      ip1 = (i==nm1 ? 0   : (i+1));  // periodic stuff
      im1 = (i==0   ? nm1 : (i-1));
      for(j=0; j<=gg->n[grid]-1; j++)
	{
	  jp1 = (j==nm1 ? 0   : (j+1));  // periodic stuff
	  jm1 = (j==0   ? nm1 : (j-1));
	  for(k=0; k<=gg->n[grid]-1; k++)
	    {
	      kp1 = (k==nm1 ? 0   : (k+1));  // periodic stuff
	      km1 = (k==0   ? nm1 : (k-1));
	      
	      nb_sum = 
		gg->fi[grid][i  ][j  ][km1] + gg->fi[grid][i  ][j  ][kp1] +
		gg->fi[grid][i  ][jm1][k  ] + gg->fi[grid][i  ][jp1][k  ] +
		gg->fi[grid][im1][j  ][k  ] + gg->fi[grid][ip1][j  ][k  ];
	      norm_rho = 
		sqrt(pow2((gg->rho[gg->g][im1][j  ][k  ] - 
			   gg->rho[gg->g][i][j][k])/gg->rho[gg->g][i][j][k]) + 
		     pow2((gg->rho[gg->g][i  ][jm1][k  ] - 
			   gg->rho[gg->g][i][j][k])/gg->rho[gg->g][i][j][k]) + 
		     pow2((gg->rho[gg->g][i  ][j  ][km1] - 
			   gg->rho[gg->g][i][j][k])/gg->rho[gg->g][i][j][k]));
	      if(!grid)
		{
		  fprintf(fp, "%g  %g  %g  %g  %g  %g  %g  %g  %g  %g\n", 
			  (((double)i * gg->h[grid] + gg->h[grid]/2.0)) * par->hsmall, 
			  (((double)j * gg->h[grid] + gg->h[grid]/2.0)) * par->hsmall, 
			  (((double)k * gg->h[grid] + gg->h[grid]/2.0)) * par->hsmall, 
			  gg->rho[grid][i][j][k]/rho_mean,
			  gg->fi[grid][i][j][k],
			  gg->res[grid][i][j][k], 
			  0.0,
			  gg->trunc[grid][i][j][k], 
			  nb_sum, 
			  norm_rho); // - 6.0*gg->fi[grid][i][j][k])/pow2(gg->h[grid]));
		}
	      else
		{
		  fprintf(fp, "%g  %g  %g  %g  %g  %g  %g  %g  %g  %g\n", 
			  (((double)i * gg->h[grid] + gg->h[grid]/2.0)) * par->hsmall, 
			  (((double)j * gg->h[grid] + gg->h[grid]/2.0)) * par->hsmall, 
			  (((double)k * gg->h[grid] + gg->h[grid]/2.0)) * par->hsmall, 
			  gg->rho[grid][i][j][k]/rho_mean,
			  gg->fi[grid][i][j][k],
			  gg->res[grid][i][j][k], 
			  gg->source[grid][i][j][k],
			  gg->trunc[grid][i][j][k], 
			  nb_sum, 
			  norm_rho); // - 6.0*gg->fi[grid][i][j][k])/pow2(gg->h[grid]));
		}
	      mean_nabla += (nb_sum - 6.0*gg->fi[grid][i][j][k])/pow2(gg->h[grid]);
	      // Linearized around 1:
	      //mean_source += 0.5*pow2(a/par->symm_range) * 
	      //	(pow3(par->symm_assb/a)*gg->rho[grid][i][j][k]/
	      // rho_mean*(1+0.0)) ;// + 2.0*gg->fi[grid][i][j][k]);
	      mean_source = 0.5*pow2(a/par->symm_range)*gg->fi[grid][i][j][k] * 
		(-1.0 + pow2(gg->fi[grid][i][j][k]) +
		 pow3(par->symm_assb/a)*gg->rho[grid][i][j][k]/rho_mean);
	    }
	}
    }

  fclose(fp);

  fprintf(stderr, "Mean nabla = %g    mean source = %g\n", 
	  mean_nabla, mean_source);
	  

  //l++;

}

//======================================
// write_temp_fifth
//======================================
/* void write_temp_fifth(int n, int grid) */
/* { */

/*   // For debuging. */

/*   int i, j, k; */
/*   FILE *fp; */
/*   char name[100]; */

/*   return; */

/*   sprintf(name, "temp.%05d.%03ld.dat", n, gg.n[grid]); */
/*   fp = fopen(name, "w"); */
    
/*   for(i=0; i<=gg.n[grid]-1; i++) */
/*     { */
/*       for(j=0; j<=gg.n[grid]-1; j++) */
/* 	{ */
/* 	  for(k=0; k<=gg.n[grid]-1; k++) */
/* 	    { */
/* 	      fprintf(fp, "%g  %g  %g\n",  */
/* 		      ((double)i * gg.h[grid] + gg.h[grid]/2.0),  */
/* 		      ((double)j * gg.h[grid] + gg.h[grid]/2.0),  */
/* 		      gg.temp[grid][i][j][k]); */
/* 	    } */
/* 	} */
/*     } */

/*   // The columns are: */
/*   //----------------- */
/*   //  1 = x */
/*   //  2 = y */
/*   //  3 = temp */

/*   fclose(fp); */
    
/* } */
