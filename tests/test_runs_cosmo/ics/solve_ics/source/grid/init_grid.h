
#include <fftw3.h>

struct grids
{
  double ***delta;
  double ***source;
  double ***phi;
  double ***lap_phi;
  double ***force_1;
  double ***force_2;
  double ***force_3;

  double ***q;
  
  fftw_complex *data;
  fftw_plan p_backward;
  fftw_plan p_forward;

};

extern void init_grid(struct params *par, struct grids *grid);



