

struct particles
{
  double *pos_1;
  double *pos_2;
  double *pos_3;
  double *vel_1;
  double *vel_2;
  double *vel_3;
  
  // Smoothing kernel
  double (*kernel)(double, double, double, double, double, double, struct params *);

};

extern void init_parts(struct params *par, struct particles *part, struct units *unit);


