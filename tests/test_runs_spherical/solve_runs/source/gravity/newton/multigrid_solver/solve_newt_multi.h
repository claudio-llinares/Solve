
#define MAX_GRIDS 20
#define W_SOR  1.34

#define X 0
#define Y 1
#define Z 2

#define MINUS 0
#define PLUS  1

#define NGS   2

// two colors
//int sx_first[8] = {0, 1, 1, 0,  1, 0, 0, 1};
//int sy_first[8] = {0, 1, 0, 1,  0, 1, 0, 1};
//int sz_first[8] = {0, 0, 1, 1,  0, 0, 1, 1};

// eight colors
extern int sx_first[8];
extern int sy_first[8];
extern int sz_first[8];

// lexicographic
extern int sx_first_lexi[8];
extern int sy_first_lexi[8];
extern int sz_first_lexi[8];

// Functions:
//===========
extern int solve_fifth_multigrid(struct params *par, double ***dens, double ***phi, double box, double grid, double a);
extern void allocate_grids_fifth(struct params *par, struct grids_gravity *gg, double ***dens, double ***phi);
extern void deallocate_grids_fifth(struct params *par, struct grids_gravity *gg, double ***dens, double ***phi);
extern void go_down_fifth(struct grids_gravity *gg, double a);

