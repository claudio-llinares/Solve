

extern void gs_mm(struct grids_gravity_mm *gg);
extern void residual_mm(struct grids_gravity_mm *gg, int grid, double ***res_phi, double ***res_chi);
extern void trunc_err_mm(struct grids_gravity_mm *gg, int grid, double ***trunc_phi, double ***trunc_chi);
