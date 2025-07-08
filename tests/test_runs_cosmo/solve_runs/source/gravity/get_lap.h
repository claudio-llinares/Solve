
// Routines to calculate laplacian of potential for many gravities.

//extern void get_lap(double ***source, double ***phi, double ***lap_phi, double a, fftw_complex *data, fftw_plan p_forward, fftw_plan p_backward, double power_a0);

extern void calc_lap(struct params *par, double ***phi, double ***lap_phi);


