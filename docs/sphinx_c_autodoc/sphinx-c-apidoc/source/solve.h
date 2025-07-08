#include <stdio.h>

/** Structure with input parameters and fields */
struct params      ///< Input parameters and some physical constants (including expansion factor).
{
  char file_in_parts[1000];
  char file_out[1000];
  double box;     /**< Mpc/h and then converted to Mpc in init_physics. */
  int grid;       ///< Grid size (the code is not tested with this number being different than the number of particles).
  int grid_p;     ///< Grid of particles (cubic root of total number of particles).
  double omegam;
  double omegal;
  double hsmall;
  double z_init;    
  double z_end;
  int nsteps;       ///< Number of time steps.
  int nthreads;     ///< Number of threads. (interesting only if OPENMP).
  int n_outputs;
  int verbose;

  int nparts;       ///< Total number of particles.
  
  int static_test;

  // Gravity
  int gravity;     /**<\n 
		      0 -> Newton
		      1 -> MOND with mass term.   */
  
  // MOND with mass parameters:
  double mm_a0;    ///< MOND a0 in km/sec^2 and converted to Mpc/Gyr^2.
  double mm_M;     ///< MOND mass in h/Mpc and converted to 1/Mpc.
  double mm_KB;    ///< MOND constant.
  int verbose_mm;
  int verbose_mm_coarse;

  double a;        ///< Expansion factor.

  double H0;

  double t_of_a[1000];

  // Derived quantities:
  double cte_poisson;

};

struct units
{
  double pi;
  double twopi;
  double clight;
  double kg2msun;
  double cm2mpc;
  double m2mpc;
  double mpc2m;
  double sec2gyr;
  double gyr2sec;
  double km2mpc;
  double mpc2km;

  double mpc2km_over_gyr2sec;

  double H0;        // H0 real = 100 hsmall  gyr^-1
  double cte_poisson;
  
};

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define pow5(x) ((x)*(x)*(x)*(x)*(x))
#define pow6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define pow7(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x))
