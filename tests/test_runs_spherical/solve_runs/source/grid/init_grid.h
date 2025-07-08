
#include <fftw3.h>

struct grids
{
  // Overdensity:
  double ***delta;

  // Newtonia fields:
  double ***source;
  double ***phi;
  double ***force_1;
  double ***force_2;
  double ***force_3;

  // MOND with mass fields:
  double ***chi_bar_mm;
  double ***phi_bar_mm;
  double ***force_chi_bar_mm_1;
  double ***force_chi_bar_mm_2;
  double ***force_chi_bar_mm_3;
  double ***force_phi_bar_mm_1;
  double ***force_phi_bar_mm_2;
  double ***force_phi_bar_mm_3;

  // For doing FFTs:
  fftw_complex *data;
  fftw_plan p_backward;
  fftw_plan p_forward;

};

extern void init_grid(struct params *par, struct grids *grid);



