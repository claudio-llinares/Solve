#include "../../solve.h"
#include "../../grid/init_grid.h"

//=====================================
//  store_fields_newton
//=====================================
void store_fields_newton(struct params *par, struct grids *grid)
{
  //
  // Copy rho, phi and f into rho_n, phi_n, f_n. //

  int i, j, k;

  for(i=0; i<par->grid; i++)
    for(j=0; j<par->grid; j++)
      for(k=0; k<par->grid; k++)
	{
	  //grid->delta_n[i][j][k] = grid->delta[i][j][k];
	  //grid->phi_n[i][j][k] = grid->phi[i][j][k];
	  grid->force_n_1[i][j][k] = grid->force_1[i][j][k];
	  grid->force_n_2[i][j][k] = grid->force_2[i][j][k];
	  grid->force_n_3[i][j][k] = grid->force_3[i][j][k];
	}

}

//=====================================
//  store_fields_fifth
//=====================================
void store_fields_fifth(struct params *par, struct grids *grid)
{
  //
  // Copy rho, phi and f into rho_n, phi_n, f_n.
  //

  int i, j, k;

  for(i=0; i<par->grid; i++)
    for(j=0; j<par->grid; j++)
      for(k=0; k<par->grid; k++)
	{
	  //rho_f(ind_cell(i)) = rho(ind_cell(i))        
	  //phi_f[i][j][k] = phi[i][j][k];
	  grid->force_s_1[i][j][k] = grid->force_1[i][j][k];
	  grid->force_s_2[i][j][k] = grid->force_2[i][j][k];
	  grid->force_s_3[i][j][k] = grid->force_3[i][j][k];
	}
}             


//=====================================
//  recover_fields_newton
//=====================================
void recover_fields_newton(struct params *par, struct grids *grid)
{
  //
  // Copy rho, phi and f into rho_n, phi_n, f_n.
  //

  int i, j, k;

  for(i=0; i<par->grid; i++)
    for(j=0; j<par->grid; j++)
      for(k=0; k<par->grid; k++)
	{
	  //rho[i][j][k] = rho_n[i][j][k];
	  //phi[i][j][k] = phi_n[i][j][k];
	  grid->force_1[i][j][k] = grid->force_n_1[i][j][k];
	  grid->force_2[i][j][k] = grid->force_n_2[i][j][k];
	  grid->force_3[i][j][k] = grid->force_n_3[i][j][k];
	}

}


//=====================================
//  recover_fields_fifth
//=====================================
void recover_fields_fifth(struct params *par, struct grids *grid)
{
  //
  // Copy rho, phi and f into rho_n, phi_n, f_n.
  //

  int i, j, k;

  for(i=0; i<par->grid; i++)
    for(j=0; j<par->grid; j++)
      for(k=0; k<par->grid; k++)
	{
	  //phi[i][j][k] = phi_f[i][j][k];
	  grid->force_1[i][j][k] = grid->force_s_1[i][j][k];
	  grid->force_2[i][j][k] = grid->force_s_2[i][j][k];
	  grid->force_3[i][j][k] = grid->force_s_3[i][j][k];
	}

}



