
extern struct param_p2rho
{
  char prefix[1000];
  double box; // Box size (Mpc/h)
  int grid;     // Number of grid cells per dimension
    
  double omega_0;     // Omega matter
  double omega_l;   // OmLambda0  (derived from curvature)
  double hsmall; // hsmall
  int model;   // model for the power spectrum
               //     0 -> Formulas from Klypin.
               //     1 -> Power law.
               //     2 -> Table with the transfter function given by cmbfast.
               //     4 -> Table with normalized power spectrum given by CAMB.
               //     5 -> Table with normalized power spectrum given by grafic.
  long Nseed;  // Random seed (Integer  1-2^31)
  int show_spectr;  // flag for showing the spectrum.
  int test;    // if 1-> output the calculated power spectrum
  int mean_delta;
  int show_sigma_8;
  int nthreads;
  // Derived from these:
  double z;
  double a;
  double t;
  double Delta;
  int test_corr; // if 1-> output the crrelation function
} par;


extern void p2rho_lib(double ***gr, double ***gi, struct param_p2rho par, fftw_complex *data, fftw_plan p_forward, fftw_plan p_backward);
extern void spectr(double ***gr, double ***gi, struct param_p2rho par);
extern void get_delta(double ***gr, double ***gi, fftw_complex *data, fftw_plan p_forward, struct param_p2rho par);
