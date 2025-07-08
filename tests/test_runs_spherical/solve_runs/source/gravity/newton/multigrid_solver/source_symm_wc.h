
extern void diff_operator_symm_wc(struct params *par, double ***phi, double h, int i, int j, int k, int ip1, int im1, int jp1, int jm1, int kp1, int km1, long int g, double *operator, double *doperator_dw);

extern double l_symm_wc(struct params *par, double nabla_phi, double phi, double rho, double source, double a);
extern double l_prime_symm_wc(struct params *par, double phi, double rho, double source, double h, double a);

extern double dchange_domega_symm_wc(double change, double omega);
extern double change_symm_wc(double change, double omega);
extern double inverse_change_symm_wc(double change, double omega);
