

extern double hubble(struct params *par, struct units *unit, double a);
extern double hubble2(struct params *par, struct units *unit, double a);
extern double hubble_dot(struct params *par, struct units *unit, double a);
extern double a_dot(struct params *par, struct units *unit, double a);

extern double integrate_chi_of_a(struct params *par, double a);
extern double integrand_chi(double a, void * params);


