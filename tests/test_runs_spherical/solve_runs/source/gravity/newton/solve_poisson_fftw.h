
extern int solve_poisson_fftw(struct params *par, struct grids *grid, struct units *unit);

extern void transforms(double ***gr, double ***gi, double box, int grid, fftw_complex *data, fftw_plan p_backward);
extern void anti_transforms(double ***gr, double ***gi, double box, int grid, fftw_complex *data, fftw_plan p_forward);
