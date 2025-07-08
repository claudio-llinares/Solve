

struct params
{
  char file_in[1000];
  char file_out[1000];
  double box;  // mpc/h
  int grid;  // grid to use (a finer one).
  int grid_p;  // number of particles per dimenssion.
  double omegam;
  double omegar;
  double omegal;
  double hsmall;
  double z_init;  // if -1 -> use value from file (for calculations with output of this program.  Otherwise, p2rho_2d doesn't know about time).
  double z_end;
  int nthreads;    // number of threads. (interesting only if OPENMP).
  long Nseed;         // Initial seed for the Bob realizations.
                      // To be used with numerical recipies, so it must be <0.
  int nsteps_nbody;     // for getting non-sinchronized potential and q.
  double sign_delta;

  double mass;
  int nparts;
  double H0;
  double clight;
  double a_init;
  double a_end;
};


struct units
{
  double pi;
  double twopi;
  double G_mks;
  double clight;
  double kg2msun;
  double m2mpc;
  double mpc2m;
  double sec2gyr;
  double gyr2sec;
  double km2mpc;
  double mpc2km;
  double cte_poisson;  // 3/2 omega0 H0
  double fourpiG;     // 4 \pi G (for phantom stuff)
  double H0;      // H0 real = 100 hsmall  gyr^-1
  double rho_crit;  // critical density (msun mpc^-3)
  double rho_mean;   // mean density (msun mpc^-3)
};

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define pow5(x) ((x)*(x)*(x)*(x)*(x))
#define pow6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define pow7(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x))
