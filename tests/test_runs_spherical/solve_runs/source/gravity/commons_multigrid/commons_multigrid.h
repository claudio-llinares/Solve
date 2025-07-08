#define MAX_GRIDS          20
#define MAX_ITERATION_SYMM 5     // maximun number of iterations/.
#define MAX_ITERATION      200    // maximun number of iterations/.

struct grids_gravity
{
  long num_grids;  // Number of grids according to the finnest.
  long g;          // Grid in which we are working.
  //     0 = finest    
  long n[MAX_GRIDS];  // size of each grid.
  double h[MAX_GRIDS];  // spacing of each grid (UA).
  double ***mu[MAX_GRIDS];  // for the first fifth
  double ***rho[MAX_GRIDS];  // Grids with density (for symmetron)
  double ***source[MAX_GRIDS];  // Grids with density.
                                // In symmetron case,this have something like 0.
  double ***fi[MAX_GRIDS];      // Grids with fi.

  double ***res[MAX_GRIDS];   // residual in all cells all grids
  double ***trunc[MAX_GRIDS]; // truncation error in all cells all grids

  double ***temp[MAX_GRIDS];   // to make the interpolation to fine grids.
  double ***temp2[MAX_GRIDS];  // to make the interpolation to fine grids.

  double old_resid[MAX_GRIDS], cur_resid[MAX_GRIDS], max_resid[MAX_GRIDS];  // integrated values.
  double old_trunc[MAX_GRIDS], cur_trunc[MAX_GRIDS], max_trunc[MAX_GRIDS];

  // MOND with mass arrays:
  double ***chi_bar_mm[MAX_GRIDS],      ***phi_bar_mm[MAX_GRIDS];
  double ***lap_chi_bar_mm[MAX_GRIDS],  ***lap_phi_bar_mm[MAX_GRIDS];
  double ***res_chi_bar_mm[MAX_GRIDS],  ***res_phi_bar_mm[MAX_GRIDS];
  double ***temp_chi_bar_mm[MAX_GRIDS], ***temp_phi_bar_mm[MAX_GRIDS];

  double (*l_operator_mm)(struct params *par, struct grids_gravity *gg, int grid, double source, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1);
  double (*l_prime_operator_mm)(struct params *par, double phi, double rho, double source, double h, double a);
  double (*L_operator_mm)(int i, int j, int k, double ***phi, double h, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int n);

  // Symmetron operators:
  double (*l_fifth_operator)(struct params *par, struct grids_gravity *gg, int grid, double source, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1);
  double (*l_prime_fifth_operator)(struct params *par, double phi, double rho, double source, double h, double a);
  double (*L_operator)(int i, int j, int k, double ***phi, double h, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int n);
};



