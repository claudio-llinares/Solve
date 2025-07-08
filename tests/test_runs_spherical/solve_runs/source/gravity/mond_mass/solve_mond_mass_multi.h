// Number of pointers to grids (not necessarily allocated).
#define MAX_GRIDS 20   

// For evaluating the MOND operator:
#define X 0
#define Y 1
#define Z 2
#define MINUS 0
#define PLUS  1

#define NGS   2              // Number of Gauss-Seidel iterations per level.  
#define MAX_ITERATION 200    // Maximun number of complete V cycles.

struct grids_gravity_mm
{
  struct params *par;          // Pointer to the parameters of the main code
                               // so we have access to things like a0 and a.

  int verbose;
  int verbose_coarse;

  double cte_poisson_over_a;
  double aa0;                  // MOND parameters
  double M2;                   // (we copy them from par for gaining speed).
  double KB;

  int num_grids;              // Number of grids according to the finnest.
  int g;                      // Grid in which we are working.
                              //     0 = finest (so what people look V cycle, is actually a Lambda).   

  int    n[MAX_GRIDS];        // Size of each grid.
  double h[MAX_GRIDS];        // Spacing of each grid (UA).
  double fourh[MAX_GRIDS];    // This one is optimization

  double ***rho[MAX_GRIDS];   // Density.

  //double old_resid[MAX_GRIDS], cur_resid[MAX_GRIDS], max_resid[MAX_GRIDS];  // integrated values.
  //double old_trunc[MAX_GRIDS], cur_trunc[MAX_GRIDS], max_trunc[MAX_GRIDS];

  // MOND with mass arrays:
  double ***phi[MAX_GRIDS];        // For phi
  double ***lap_phi[MAX_GRIDS];
  double ***res_phi[MAX_GRIDS];
  double ***source_phi[MAX_GRIDS];
  double ***temp_phi[MAX_GRIDS];
  double ***temp2_phi[MAX_GRIDS];
  double ***trunc_phi[MAX_GRIDS];  // truncation error in all cells all grids

  double ***chi[MAX_GRIDS];        // For chi
  double ***lap_chi[MAX_GRIDS];
  double ***res_chi[MAX_GRIDS];
  double ***source_chi[MAX_GRIDS];
  double ***temp_chi[MAX_GRIDS];
  double ***temp2_chi[MAX_GRIDS];
  double ***trunc_chi[MAX_GRIDS];  // truncation error in all cells all grids


};



// Functions:
//===========
extern int solve_mm_multigrid(struct params *par, double ***dens, double ***chi, double ***phi, double box, double grid, double a);
