
extern void copy_q(struct params *par, struct grids *grid);
extern void copy_force_phi(struct params *par, struct grids *grid);

extern FILE *fopen_check(char *file, char *type);

extern double modulus(double x, double y);
extern double euclides2(double x, double y, double z);
extern double euclides(double x, double y, double z);
extern double dist(double x, double y, double z);
extern void change_system(double vec[3][3], double *x, double *y, double *z);
extern void sort_3_values(double *a);
extern void index_1d_to_3d(double q, int *i, int *j, int *k, int n);

#define max(a,b) ((a) > (b) ? (a) : (b))

