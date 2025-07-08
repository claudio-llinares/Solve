

// Second order:
//==============
extern double l_fifth(struct params *par, struct grids_gravity *gg, int grid, double source, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1);
extern double l_prime_fifth(struct params *par, double phi, double rho, double source, double h, double a);
extern double L(int i, int j, int k, double ***phi, double h, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int n);

// Fourth order:
//==============
extern double l_fifth_fourth(struct params *par, struct grids_gravity *gg, int grid, double source, int i, int j, int k, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1);
extern double l_prime_fifth_fourth(struct params *par, double phi, double rho, double source, double h, double a);
extern double L_fourth(int i, int j, int k, double ***phi, double h, double a, int ip1, int im1, int jp1, int jm1, int kp1, int km1, int n);

// Background:
//============
extern double phi_b(struct params *par, double a);   // symmetron
