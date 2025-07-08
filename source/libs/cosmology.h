

extern double hubble(struct params *par, struct units *unit, double a);
extern double hubble2(struct params *par, struct units *unit, double a);
extern double hubble_dot(struct params *par, struct units *unit, double a);
extern double a_dot(struct params *par, struct units *unit, double a);
extern double F(struct params *par, struct units *unit, double a);
extern double A(struct params *par, struct units *unit, double a);

extern double integrate_time_of_a(struct params *par, double a);
extern double integrand_time(double a, void * params);


extern void tabulate_t_of_a(struct params *par);
extern double a_of_t(struct params *par, double t);


