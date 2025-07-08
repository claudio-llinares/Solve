#define MAX_GRIDS 20        /**< Number of pointers to grids (not necessarily allocated). (i.e. not important unless you want to have more than 20 levels). */

#define X 0                 /**< For evaluating the MOND operator. */
#define Y 1                 /**< For evaluating the MOND operator. */
#define Z 2                 /**< For evaluating the MOND operator. */
#define MINUS 0             /**< For evaluating the MOND operator. */
#define PLUS  1             /**< For evaluating the MOND operator. */

#define NGS   2             /**< Number of Gauss-Seidel iterations per level.   */
#define MAX_ITERATION 200   /**< Maximun number of complete V cycles before accepting that the solver will not converge. */

#define GO_UP_FIRST  1      /**< default=1.   If 0 then solver does not go up in the first multigrid level (i.e. one grid scheme). */
#define GO_UP_SECOND 1      /**< default=1.   If 0 then solver does not go up in the second multigrid level (i.e. two grids scheme). **/


/** 
 * Main data structure of the multigrid solver.
 * 
 * Comments:
 *     - Note that here we have a copy of the main parameters of the code.  
 *       That is to avoid having to accessing too many structures in some of the 
 *       loops in the solver.
 **/
struct grids_gravity_mm
{
  struct params *par;          /**< Pointer to the parameter structure of the main code (see solve.h for details).*/
                              
  int verbose;
  int verbose_coarse;

  double cte_poisson_over_a;
  double aa0;                 /**< MOND parameter a0 */
  double a2M2;                /**< MOND parameter M^2 (multipled by a^2). */
  double KB;                  /**< MOND parameter KB */

  int num_grids;              /**< Number of grids according to the finnest. */
  int g;                      /**< Grid in which we are working.  0 = finest (so what people call V cycle, is actually a Lambda). */  

  int    n[MAX_GRIDS];        /**< Size of each grid. */
  double h[MAX_GRIDS];        /**< Spacing of each grid. */
  double fourh[MAX_GRIDS];    /**< This one is optimization */

  double ***rho[MAX_GRIDS];   /**< Density (or overdensity) (see how it is used in gs routine)).  */

  // MOND with mass arrays:
  double ***phi[MAX_GRIDS];        /**< Solution for the Newtonia potential \f$\tilde{\Phi}\f$.*/
  double ***chi[MAX_GRIDS];        /**< Solution for the MOND equation (\f$\chi\f$). */
  double ***res_phi[MAX_GRIDS];    /**< Residual of the equation for phi. */
  double ***res_chi[MAX_GRIDS];    /**< Residual of the equation for chi. */
  double ***source_phi[MAX_GRIDS]; /**< Source of the equation for coarse grids in the FAS algorithm */
  double ***source_chi[MAX_GRIDS];
  double ***temp_phi[MAX_GRIDS];
  double ***temp_chi[MAX_GRIDS];
  double ***temp2_phi[MAX_GRIDS];
  double ***temp2_chi[MAX_GRIDS];
  double ***trunc_phi[MAX_GRIDS];  /**< truncation error in all cells all grids. */
  double ***trunc_chi[MAX_GRIDS];  /**< truncation error in all cells all grids. */


};


// Functions:
//===========
extern int solve_multigrid_mm(struct params *par, double ***dens, double ***chi, double ***phi, double box, double grid, double a);
