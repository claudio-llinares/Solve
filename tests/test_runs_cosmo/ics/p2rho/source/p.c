#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include "libs/cuad_min.h"

struct param_p
{
  // For lcdm
  double lcdm_n;   // slope of te primordial power spectrum
  double lcdm_t_cmb;   // T_{CMB} in formula 24.
  double lcdm_omega_b;   // Omega barion
  double lcdm_norm;    // normalization.  assuming D(1/1001)=1.
  // For power law
  double power_slope;  // slope
  double power_norm;   // value at k=.1 (h/Mpc)^-3
  // For the table
  char table[1000];   // file name of with the table.
  int col;         // column where to find the transfer function to use.
                      // cmbfast gives:
                      //     2 -> cold dark matter
                      //     3 -> baryon
                      //     4 -> foton
                      //     5 -> massless neutrino
                      //     6 -> massive neutrino
  double norm;     // To change normalization.
                   // Use 1 to get the standard case.
} par_p;

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))

double model_lcdm(double wk, double omega_0, double hsmall);
double model_power_law(double wk);
double model_table(double wk, double hsmall, char *prefix, int model, double box, double grid);
double model_table_camb(double wk, double hsmall, char *prefix, int model, double box, double grid);
double model_table_grafic(double wk, double hsmall, char *prefix, int model, double box, double grid);

int get_index(double wk, double *k, int n, double hsmall);

//=============================
// P
//=============================
double P(double wk, double omega_0, double hsmall, char *prefix, int model, double box, double grid)
{
  // Calculates P according to what you give in model.
  // 
  //
  // Arguments:
  //            wk = ????????????
  //            model = model to use.
  //                    0 -> Formula 24 in Anatoly.
  //                    1 -> Power law  with variable slope and normalization.
  //                           The normalization is the value at w=.01
  //                    2 -> Table with the transfter function given by cmbfast.
  //                         (the correct normalization should come in the 
  //                         table).  Used for calculating CMB accelerations 
  //                         for the thesis.
  //                    3 -> Table given already in the vectors (no files
  //                         are read.  The normalization shoudl be the correct
  //                         one.
  //                    4 -> Table with normalized power spectrum given by CAMB.
  //                          
  //
  // Comments:
  //            - Reads parameters of P from file P.input or P_table.input
  //              There are parameters for each model.  In case you don't 
  //              need for a certain model, put 0 (or whatever, will not
  //              be used).
  // CHECK WK!!!!!!!!!!!!!!!!!!!!!
  
  static int first = 1;
  static FILE *fp_in;

  char file[1000];

  int return_fscanf;

  // Initialize models:
  //===================
  if(first)
    {
      // Read parameters:
      //=================
      if(model==0 || model==1)
	{
	  fprintf(stderr, "Reading power spectrum parameters for model %d from file %s\n", model, file);
	  sprintf(file, "p.input");
	  if( (fp_in=fopen(file, "r")) == NULL)
	    {
	      fprintf(stderr, "Problems opening p.input\n");
	      exit(0);
	    }
	  return_fscanf = fscanf(fp_in, "%lf", &par_p.lcdm_n);
	  return_fscanf = fscanf(fp_in, "%lf", &par_p.lcdm_t_cmb);
	  return_fscanf = fscanf(fp_in, "%lf", &par_p.lcdm_omega_b);
	  return_fscanf = fscanf(fp_in, "%lf", &par_p.lcdm_norm);
	  return_fscanf = fscanf(fp_in, "%lf", &par_p.power_slope);
	  return_fscanf = fscanf(fp_in, "%lf", &par_p.power_norm);
	  fclose(fp_in);
	}
      else if(model==2 || model==3)
	{
	  fprintf(stderr, "Reading power spectrum parameters for model %d from file p_table.input\n", model);
	  //sprintf(file, "%s.p_table.input", prefix);
	  sprintf(file, "p_table.input");
	  if( (fp_in=fopen(file, "r")) == NULL)
	    {
	      fprintf(stderr, "Problems opening p_table.input\n");
	      exit(0);
	    }
	  return_fscanf = fscanf(fp_in, "%s", par_p.table);
	  return_fscanf = fscanf(fp_in, "%d", &par_p.col);
	  return_fscanf = fscanf(fp_in, "%lf", &par_p.norm);
	  fprintf(stderr, "par_p.table = %s\n", par_p.table);
	  fprintf(stderr, "par_p.col   = %d\n", par_p.col);
	  fprintf(stderr, "par_p.norm  = %lg\n", par_p.norm);
	  fclose(fp_in);
	}
      else if(model==4 || model==5 )
	{
	  fprintf(stderr, "Reading power spectrum parameters for model %d from file p_table.input\n", model);
	  //sprintf(file, "%s.p_table.input", prefix);
	  sprintf(file, "p_table.input");
	  if( (fp_in=fopen(file, "r")) == NULL)
	    {
	      fprintf(stderr, "Problems opening p_table.input\n");
	      exit(0);
	    }
	  return_fscanf = fscanf(fp_in, "%s", par_p.table);
	  fprintf(stderr, "par_p.table = %s\n", par_p.table);
	  fclose(fp_in);
	}
      first = 0;
    }

  
  //return model_table_camb(wk, hsmall, prefix, model, box, grid);

  switch(model)
    {
    case 0:
      return model_lcdm(wk, omega_0, hsmall);
      break;
    case 1:
      return model_power_law(wk);
      break;
    case 2:
      return model_table(wk, hsmall, prefix, model, box, grid);
      break;
    case 3:
      return model_table(wk, hsmall, prefix, model, box, grid);
    case 4:
      return model_table_camb(wk, hsmall, prefix, model, box, grid);
      break;
    case 5:
      return model_table_grafic(wk, hsmall, prefix, model, box, grid);
      break;
    default:
      fprintf(stderr, "Wrong model for P(k)... aborting\n");
      exit(0);
    }
  
}

//============================
// model_lcdm
//============================
double model_lcdm(double wk, double omega_0, double hsmall)
{
  // Suport routine for LCDM model.
  // Return P(k) for LCDM.
  
  static double a1;
  static double a2;
  static double alpha;
  static double c;

  double q, t, p;
  
  static int first = 1;

  if(first)
    {
      // Define some constants:
      //=======================
      a1 = pow(46.9 * omega_0 * pow2(hsmall), 0.67) * 
	(1.0 + pow(32.1*omega_0 * pow2(hsmall), -.532));
      
      a2 = pow(12.0*omega_0*pow2(hsmall), .424) *
	(1.0 + pow(45.0*omega_0*pow2(hsmall), -.582));
      
      alpha = pow(a1, -par_p.lcdm_omega_b/omega_0) *
	pow(a2, - pow3(par_p.lcdm_omega_b/omega_0));
      
      c = pow2(par_p.lcdm_t_cmb/2.7) / 
	(omega_0*pow2(hsmall) * sqrt(alpha) * 
	 pow(1.0 - par_p.lcdm_omega_b/omega_0, 0.6));

      first = 0;

    }

  q = wk * c;

  t = log(1.0 + 2.34 * q)/(2.34*q) *
    pow((1.0 + q*(13.0 + q*(pow2(10.5) + q*(pow3(10.4) + q*pow4(6.51))))), -.25);

  p = par_p.lcdm_norm * pow(wk, par_p.lcdm_n) * pow2(t);

  return p;


}

//=========================
// model_power_law
//=========================
double model_power_law(double wk)
{
  // Suport routine for power law model.
  // Returns:
  //          P(wk) = a * wk^n
  //          par_p.power_norm = value at wk=.1
  //

  static double a;   // normalization
 
  static int first = 1;

  
  if(first)
    {
      a = par_p.power_norm / pow(.1, par_p.power_slope);
      first = 0;
    }


  return a * pow(wk, par_p.power_slope);

}

//=========================
// model_table
//=========================
double model_table(double wk, double hsmall, char *prefix, int model, double box, double grid)
{
  // Suport routine for model given in a table.
  //
  // Units:
  //     - The expected table should be output of cmbfast with a transfer 
  //       function with units of h/Mpc (for T and k).
  //     - We get rid of the small h in this function when reading.
  //
  // Comments:
  //     - Normalization is taken from the table (no further changes are made).
  // 
  // Returns:
  //     - P in Mpc^-3 (without the small h).
  //

  static int first = 1;

  static double *k, *p;   // the table
  static int n;  // number of elements in the table
  static double k_min, k_max;   // minimun and maximun frequencies in the table
  static double pi;

  int index;

  double pr;  // p to return.

  FILE *fp, *fp_out_table;
  char file_name[1000];
  char line[1000];
  double x[6];     // to read the table

  double lim1, li, lip1, lip2;
  double wkl, klm1, kl, klp1, klp2;  // for logaritm values

  // To use extrapolated values:
  int j, l;
  double x_fit[100], y_fit[100];
  static double a_fit, b_fit;

  char *return_fgets;

  // Reads the table:
  //=================
  if(first)
    {
      pi = 4.0*atan(1.0);

      // Count lines:
      if( (fp=fopen(par_p.table, "r")) == NULL)
	{
	  fprintf(stderr, "Problems opening %s\n", par_p.table);
	  exit(0);
	}
      n = 0;
      while(!feof(fp))
	{
	  return_fgets = fgets(line, 1000, fp);
	  if(feof(fp))
	    break;
	  n++;
	}
      fclose(fp);

      // Allocates table
      k = (double*)calloc(n, sizeof(double));
      p = (double *)calloc(n, sizeof(double));

      // Reads the table
      fprintf(stderr, "Reading transfer function from file %s\n", par_p.table);
      if( (fp=fopen(par_p.table, "r")) == NULL)
	{
	  fprintf(stderr, "Problems opening %s\n", par_p.table);
	  exit(0);
	}
      if( (fp_out_table=fopen("p_input.dat", "w")) == NULL)
	{
	  fprintf(stderr, "Problems opening %s\n", file_name);
	  exit(0);
	}
      n = 0;
      k_min = 1.0e+50;
      k_max = -1.0;
      while(!feof(fp))
	{
	  return_fgets = fgets(line, 1000, fp);
	  if(feof(fp))
	    break;
	  sscanf(line, "%lf%lf%lf%lf%lf%lf", 
		 &x[0], &x[1], &x[2], &x[3], &x[4], &x[5]);
	  // Asume the file has:
	  // [k] = h/Mpc
	  k[n] = x[0] * hsmall;
	  if(model==2)
	    p[n] = pow(x[par_p.col],2) * (x[0]*hsmall) * (2.0*pow2(pi)) * 2.4e-9 * par_p.norm;
	  else
	    p[n] = x[1] / pow(hsmall, 3);
	  fprintf(fp_out_table, "%lg  %lg\n", 
	    k[n]/hsmall, 
	    p[n] * pow3(hsmall));
	  // get max/min
	  if(k_min>k[n])
	    k_min = k[n];
	  if(k_max<k[n])
	    k_max = k[n];
	  n++;
	}
      fclose(fp);
      fclose(fp_out_table);

      fprintf(stderr, "model_table: k_min  = %lf\n", k_min/hsmall);
      fprintf(stderr, "model_table: k_max  = %lf\n", k_max/hsmall);
      fprintf(stderr, "model_table: k[n-1] = %lf\n", k[n-1]/hsmall);
      fprintf(stderr, "model_table: n      = %d\n", n);

      // Fit power law to high frequencies:
      //===================================
      if(k[n-1] < (sqrt(3.0)*pi*grid)/box)
	{
	  // copy into new arrays:
	  l = 0;
	  fprintf(stderr, "Fitting power law to this data:\n");
	  // 314
	  for(j=n-10; j<n; j++)
	    {
	      x_fit[l] = log10(k[j]);
	      y_fit[l] = log10(p[j]);
	      fprintf(stderr, "%le  %le\n", x_fit[l]/hsmall, y_fit[l]);
	      l++;
	    }
	  // get parameters:
	  cuad_min(x_fit, y_fit, 10, &a_fit, &b_fit);
	  fprintf(stderr, "Power law fitted to input power spectrum.  a = %lf   b = %lf\n", a_fit, b_fit);
	}

      first = 0;
    }

  //fprintf(stderr, "%g\n", wk/hsmall);

  if(wk==0)
    {
      return 0.0;
    }
  else
    {
      // n-2 because we need one more point for the interpolation
      if(wk>k[n-2])
	{
	  //fprintf(stderr, "1 %lf  %lf  %d\n", wk/hsmall, k[n-1]/hsmall, n);
	  // Gets extrapolated value:
	  //=========================
	  //fprintf(stderr, "Using extrapolated value\n");
	  pr = pow(10, a_fit) * pow(wk, b_fit);
	}
      else
	{
	  //fprintf(stderr, "0 %lf  %lf  %d\n", wk/hsmall, k[n-1]/hsmall, n);
	  // Get the index for the frequency:
	  //=================================
	  // This is if frequencies are equi-spaced
	  //index = (int)((wk - k_min) * (n-1) / (k_max-k_min) + .5);
	  index = get_index(wk, k, n, hsmall);
	  
	  // Interpolate the value:
	  //=======================
	  // Make linear interpolation (pr = p to return)
	  // Third order interpolation:
	  wkl  = log10(wk);
	  klm1 = log10(k[index-1]);
	  kl   = log10(k[index]);
	  klp1 = log10(k[index+1]);
	  klp2 = log10(k[index+2]);
	  
	  lim1 = ((wkl-kl) *
		  (wkl-klp1) *
		  (wkl-klp2)) /
	    ((klm1-kl) *
	     (klm1-klp1) *
	     (klm1-klp2));
	  
	  li = ((wkl-klm1) *
		(wkl-klp1) *
		(wkl-klp2)) /
	    ((kl-klm1) *
	     (kl-klp1) *
	     (kl-klp2));
	  
	  lip1 = ((wkl-klm1) *
		  (wkl-kl) *
		  (wkl-klp2)) /
	    ((klp1-klm1) *
	     (klp1-kl) *
	     (klp1-klp2));
	  
	  lip2 = ((wkl-klm1) *
		  (wkl-kl) *
		  (wkl-klp1)) /
	    ((klp2-klm1) *
	     (klp2-klp1) *
	     (klp2-kl));
	  
	  pr = lim1 * log10(p[index-1]) +
	    li   * log10(p[index]) +
	    lip1 * log10(p[index+1]) +
	    lip2 * log10(p[index+2]);
	  pr = pow(10, pr);
	  
	  // If is negative, use linear interpolation.
	  if(pr<=0)
	    {
	      //pr = (wk - k[index]) / (k[index+1] - k[index]) * 
	      //(p[index+1] - p[index]) + p[index];
	      // Makes interpolation in log space:
	      pr = (log10(wk) - log10(k[index])) / (log10(k[index+1]) - log10(k[index])) * 
		(log10(p[index+1]) - log10(p[index])) + log10(p[index]);
	      pr = pow(10, pr);
	    }
	}
      
      //pr = p[index];
      
      //fprintf(stderr, "%lg  %lg  %lg  %lg\n", 
      //	      k[index]/hsmall,
      //	      wk/hsmall, 
      //	      k[index+1]/hsmall,
      //	      pr * pow3(hsmall));
      
      return pr;
    }

}

//=========================
// model_table_camb
//=========================
double model_table_camb(double wk, double hsmall, char *prefix, int model, double box, double grid)
{
  // Suport routine for model given in a table.
  //
  // Units:
  //     - The expected table should be output of cmbfast with a transfer 
  //       function with units of h/Mpc (for T and k).
  //     - We get rid of the small h in this function when reading.
  //
  // Comments:
  //     - Normalization is taken from the table (no further changes are made).
  // 
  // Returns:
  //     - P in Mpc^-3 (without the small h).
  //

  static int first = 1;

  static double *k, *p;   // the table
  static double *log10_k, *log10_p;
  static int n;  // number of elements in the table
  static double k_min, k_max;   // minimun and maximun frequencies in the table
  static double pi;

  int index;

  double pr;  // p to return.

  FILE *fp, *fp_out_table;
  char file_name[1000];
  char line[1000];
  double x[6];     // to read the table

  double lim1, li, lip1, lip2;
  double wkl, klm1, kl, klp1, klp2;  // for logaritm values

  // To use extrapolated values:
  int j, l;
  double x_fit[100], y_fit[100];
  static double a_fit, b_fit;

  char *return_fgets;

  // Reads the table:
  //=================
  if(first)
    {
      pi = 4.0*atan(1.0);

      // Count lines:
      if( (fp=fopen(par_p.table, "r")) == NULL)
	{
	  fprintf(stderr, "Problems opening %s\n", par_p.table);
	  exit(0);
	}
      n = 0;
      while(!feof(fp))
	{
	  return_fgets = fgets(line, 1000, fp);
	  if(feof(fp))
	    break;
	  n++;
	}
      fclose(fp);

      // Allocates table
      k = (double*)calloc(n, sizeof(double));
      p = (double *)calloc(n, sizeof(double));
      log10_k = (double*)calloc(n, sizeof(double));
      log10_p = (double *)calloc(n, sizeof(double));

      // Reads the table
      fprintf(stderr, "Reading transfer function from file %s\n", par_p.table);
      if( (fp=fopen(par_p.table, "r")) == NULL)
	{
	  fprintf(stderr, "Problems opening %s\n", par_p.table);
	  exit(0);
	}
      if( (fp_out_table=fopen("p_input.dat", "w")) == NULL)
	{
	  fprintf(stderr, "Problems opening %s\n", file_name);
	  exit(0);
	}
      n = 0;
      k_min = 1.0e+50;
      k_max = -1.0;
      while(!feof(fp))
	{
	  return_fgets = fgets(line, 1000, fp);
	  if(feof(fp))
	    break;
	  sscanf(line, "%lf%lf", 
		 &x[0], &x[1]);
	  // Asume the file has:
	  // [k] = h/Mpc
	  k[n] = x[0] * hsmall;
	  p[n] = x[1] / pow(hsmall, 3);
	  // Calculate logarithms:
	  log10_k[n] = log10(k[n]);
	  log10_p[n] = log10(p[n]);
	  fprintf(fp_out_table, "%lg  %lg\n", 
	    k[n]/hsmall, 
	    p[n] * pow3(hsmall));
	  // get max/min
	  if(k_min>k[n])
	    k_min = k[n];
	  if(k_max<k[n])
	    k_max = k[n];
	  n++;
	}
      fclose(fp);
      fclose(fp_out_table);

      fprintf(stderr, "model_table: k_min  = %lf\n", k_min/hsmall);
      fprintf(stderr, "model_table: k_max  = %lf\n", k_max/hsmall);
      fprintf(stderr, "model_table: k[n-1] = %lf\n", k[n-1]/hsmall);
      fprintf(stderr, "model_table: n      = %d\n", n);

      // Fit power law to high frequencies:
      //===================================
      if(k[n-1] < (sqrt(3.0)*pi*grid)/box)
	{
	  // copy into new arrays:
	  l = 0;
	  fprintf(stderr, "Fitting power law to this data:\n");
	  // 314
	  for(j=n-10; j<n; j++)
	    {
	      x_fit[l] = log10(k[j]);
	      y_fit[l] = log10(p[j]);
	      fprintf(stderr, "%le  %le\n", x_fit[l]/hsmall, y_fit[l]);
	      l++;
	    }
	  // get parameters:
	  cuad_min(x_fit, y_fit, 10, &a_fit, &b_fit);
	  fprintf(stderr, "Power law fitted to input power spectrum.  a = %lf   b = %lf\n", a_fit, b_fit);
	}

      first = 0;
      fprintf(stderr, "Finish reading.\n");
    }  // if first

  //fprintf(stderr, "%g\n", wk/hsmall);

  if(wk==0)
    {
      return 0.0;
    }
  else
    {
      // n-2 because we need one more point for the interpolation
      if(wk>k[n-2])
	{
	  //fprintf(stderr, "1 %lf  %lf  %d\n", wk/hsmall, k[n-1]/hsmall, n);
	  // Gets extrapolated value:
	  //=========================
	  //fprintf(stderr, "Using extrapolated value\n");
	  pr = pow(10, a_fit) * pow(wk, b_fit);
	}
      else
	{
	  //fprintf(stderr, "0 %lf  %lf  %d\n", wk/hsmall, k[n-1]/hsmall, n);
	  // Get the index for the frequency:
	  //=================================
	  // This is if frequencies are equi-spaced
	  //index = (int)((wk - k_min) * (n-1) / (k_max-k_min) + .5);
	  index = get_index(wk, k, n, hsmall);
	  
	  // Interpolate the value:
	  //=======================
	  // Make linear interpolation (pr = p to return)
	  // Third order interpolation:
	  wkl  = log10(wk);
	  //klm1 = log10(k[index-1]);
	  //kl   = log10(k[index]);
	  //klp1 = log10(k[index+1]);
	  //klp2 = log10(k[index+2]);
	  // Oprimized version
	  klm1 = log10_k[index-1] ;
	  kl   = log10_k[index] ;
	  klp1 = log10_k[index+1];
	  klp2 = log10_k[index+2] ;
	  
	  lim1 = ((wkl-kl) *
		  (wkl-klp1) *
		  (wkl-klp2)) /
	    ((klm1-kl) *
	     (klm1-klp1) *
	     (klm1-klp2));
	  
	  li = ((wkl-klm1) *
		(wkl-klp1) *
		(wkl-klp2)) /
	    ((kl-klm1) *
	     (kl-klp1) *
	     (kl-klp2));
	  
	  lip1 = ((wkl-klm1) *
		  (wkl-kl) *
		  (wkl-klp2)) /
	    ((klp1-klm1) *
	     (klp1-kl) *
	     (klp1-klp2));
	  
	  lip2 = ((wkl-klm1) *
		  (wkl-kl) *
		  (wkl-klp1)) /
	    ((klp2-klm1) *
	     (klp2-klp1) *
	     (klp2-kl));
	  
	  //pr = lim1 * log10(p[index-1]) +
	  //  li   * log10(p[index]) +
	  //  lip1 * log10(p[index+1]) +
	  //  lip2 * log10(p[index+2]);
	  // Optimized version
	  pr = lim1 * log10_p[index-1] +
	    li   * log10_p[index] +
	    lip1 * log10_p[index+1] +
	    lip2 * log10_p[index+2];
	  pr = pow(10, pr);
	  
	  // If is negative, use linear interpolation.
	  if(pr<=0)
	    {
	      //pr = (wk - k[index]) / (k[index+1] - k[index]) * 
	      //(p[index+1] - p[index]) + p[index];
	      // Makes interpolation in log space:
	      //pr = (log10(wk) - log10(k[index])) / (log10(k[index+1]) - log10(k[index])) * 
	      //	(log10(p[index+1]) - log10(p[index])) + log10(p[index]);
	      // Optimized version
	      pr = (wkl - log10_k[index]) / (log10_k[index+1] - log10_k[index]) * 
		(log10_p[index+1] - log10_p[index]) + log10_p[index];
	      pr = pow(10, pr);
	    }
	}
      
      //pr = p[index];
      
      //fprintf(stderr, "%lg  %lg  %lg  %lg\n", 
      //	      k[index]/hsmall,
      //	      wk/hsmall, 
      //	      k[index+1]/hsmall,
      //	      pr * pow3(hsmall));
      
      return pr;
    }

}


//=========================
// model_table_camb
//=========================
double model_table_grafic(double wk, double hsmall, char *prefix, int model, double box, double grid)
{
  // Suport routine for model given in a table.
  //
  // Units:
  //     - The expected table should be output of cmbfast with a transfer 
  //       function with units of h/Mpc (for T and k).
  //     - We get rid of the small h in this function when reading.
  //
  // Comments:
  //     - Normalization is taken from the table (no further changes are made).
  // 
  // Returns:
  //     - P in Mpc^-3 (without the small h).
  //

  static int first = 1;

  static double *k, *p;   // the table
  static double *log10_k, *log10_p;
  static int n;  // number of elements in the table
  static double k_min, k_max;   // minimun and maximun frequencies in the table
  static double pi;

  int index;

  double pr;  // p to return.

  FILE *fp, *fp_out_table;
  char file_name[1000];
  char line[1000];
  double x[6];     // to read the table

  double lim1, li, lip1, lip2;
  double wkl, klm1, kl, klp1, klp2;  // for logaritm values

  // To use extrapolated values:
  int j, l;
  double x_fit[100], y_fit[100];
  static double a_fit, b_fit;

  char *return_fgets;

  // Reads the table:
  //=================
  if(first)
    {
      pi = 4.0*atan(1.0);

      // Count lines:
      if( (fp=fopen(par_p.table, "r")) == NULL)
	{
	  fprintf(stderr, "Problems opening %s\n", par_p.table);
	  exit(0);
	}
      n = 0;
      while(!feof(fp))
	{
	  return_fgets = fgets(line, 1000, fp);
	  if(feof(fp))
	    break;
	  n++;
	}
      fclose(fp);

      // Allocates table
      k = (double*)calloc(n, sizeof(double));
      p = (double *)calloc(n, sizeof(double));
      log10_k = (double*)calloc(n, sizeof(double));
      log10_p = (double *)calloc(n, sizeof(double));

      // Reads the table
      fprintf(stderr, "Reading transfer function from file %s\n", par_p.table);
      if( (fp=fopen(par_p.table, "r")) == NULL)
	{
	  fprintf(stderr, "Problems opening %s\n", par_p.table);
	  exit(0);
	}
      if( (fp_out_table=fopen("p_input.dat", "w")) == NULL)
	{
	  fprintf(stderr, "Problems opening %s\n", file_name);
	  exit(0);
	}
      n = 0;
      k_min = 1.0e+50;
      k_max = -1.0;
      while(!feof(fp))
	{
	  return_fgets = fgets(line, 1000, fp);
	  if(feof(fp))
	    break;
	  sscanf(line, "%lf%lf", 
		 &x[0], &x[1]);
	  // Asume the file has:
	  // [k] = h/Mpc
	  k[n] = x[0] * hsmall;
	  p[n] = x[1]*pow3(2.0*pi) / pow(hsmall, 3);
	  // Calculate logarithms:
	  log10_k[n] = log10(k[n]);
	  log10_p[n] = log10(p[n]);
	  fprintf(fp_out_table, "%lg  %lg\n", 
	    k[n], 
	    p[n]);
	  // get max/min
	  if(k_min>k[n])
	    k_min = k[n];
	  if(k_max<k[n])
	    k_max = k[n];
	  n++;
	}
      fclose(fp);
      fclose(fp_out_table);

      fprintf(stderr, "model_table: k_min  = %lf\n", k_min/hsmall);
      fprintf(stderr, "model_table: k_max  = %lf\n", k_max/hsmall);
      fprintf(stderr, "model_table: k[n-1] = %lf\n", k[n-1]/hsmall);
      fprintf(stderr, "model_table: n      = %d\n", n);

      // Fit power law to high frequencies:
      //===================================
      if(k[n-1] < (sqrt(3.0)*pi*grid)/box)
	{
	  // copy into new arrays:
	  l = 0;
	  fprintf(stderr, "Fitting power law to this data:\n");
	  // 314
	  for(j=n-10; j<n; j++)
	    {
	      x_fit[l] = log10(k[j]);
	      y_fit[l] = log10(p[j]);
	      fprintf(stderr, "%le  %le\n", x_fit[l]/hsmall, y_fit[l]);
	      l++;
	    }
	  // get parameters:
	  cuad_min(x_fit, y_fit, 10, &a_fit, &b_fit);
	  fprintf(stderr, "Power law fitted to input power spectrum.  a = %lf   b = %lf\n", a_fit, b_fit);
	}

      first = 0;
      fprintf(stderr, "Finish reading.\n");
    }  // if first

  //fprintf(stderr, "%g\n", wk/hsmall);

  if(wk==0)
    {
      return 0.0;
    }
  else
    {
      // n-2 because we need one more point for the interpolation
      if(wk>k[n-2])
	{
	  //fprintf(stderr, "1 %lf  %lf  %d\n", wk/hsmall, k[n-1]/hsmall, n);
	  // Gets extrapolated value:
	  //=========================
	  //fprintf(stderr, "Using extrapolated value\n");
	  pr = pow(10, a_fit) * pow(wk, b_fit);
	}
      else
	{
	  //fprintf(stderr, "0 %lf  %lf  %d\n", wk/hsmall, k[n-1]/hsmall, n);
	  // Get the index for the frequency:
	  //=================================
	  // This is if frequencies are equi-spaced
	  //index = (int)((wk - k_min) * (n-1) / (k_max-k_min) + .5);
	  index = get_index(wk, k, n, hsmall);
	  
	  // Interpolate the value:
	  //=======================
	  // Make linear interpolation (pr = p to return)
	  // Third order interpolation:
	  wkl  = log10(wk);
	  //klm1 = log10(k[index-1]);
	  //kl   = log10(k[index]);
	  //klp1 = log10(k[index+1]);
	  //klp2 = log10(k[index+2]);
	  // Oprimized version
	  klm1 = log10_k[index-1] ;
	  kl   = log10_k[index] ;
	  klp1 = log10_k[index+1];
	  klp2 = log10_k[index+2] ;
	  
	  lim1 = ((wkl-kl) *
		  (wkl-klp1) *
		  (wkl-klp2)) /
	    ((klm1-kl) *
	     (klm1-klp1) *
	     (klm1-klp2));
	  
	  li = ((wkl-klm1) *
		(wkl-klp1) *
		(wkl-klp2)) /
	    ((kl-klm1) *
	     (kl-klp1) *
	     (kl-klp2));
	  
	  lip1 = ((wkl-klm1) *
		  (wkl-kl) *
		  (wkl-klp2)) /
	    ((klp1-klm1) *
	     (klp1-kl) *
	     (klp1-klp2));
	  
	  lip2 = ((wkl-klm1) *
		  (wkl-kl) *
		  (wkl-klp1)) /
	    ((klp2-klm1) *
	     (klp2-klp1) *
	     (klp2-kl));
	  
	  //pr = lim1 * log10(p[index-1]) +
	  //  li   * log10(p[index]) +
	  //  lip1 * log10(p[index+1]) +
	  //  lip2 * log10(p[index+2]);
	  // Optimized version
	  pr = lim1 * log10_p[index-1] +
	    li   * log10_p[index] +
	    lip1 * log10_p[index+1] +
	    lip2 * log10_p[index+2];
	  pr = pow(10, pr);
	  
	  // If is negative, use linear interpolation.
	  if(pr<=0)
	    {
	      //pr = (wk - k[index]) / (k[index+1] - k[index]) * 
	      //(p[index+1] - p[index]) + p[index];
	      // Makes interpolation in log space:
	      //pr = (log10(wk) - log10(k[index])) / (log10(k[index+1]) - log10(k[index])) * 
	      //	(log10(p[index+1]) - log10(p[index])) + log10(p[index]);
	      // Optimized version
	      pr = (wkl - log10_k[index]) / (log10_k[index+1] - log10_k[index]) * 
		(log10_p[index+1] - log10_p[index]) + log10_p[index];
	      pr = pow(10, pr);
	    }
	}
      
      //pr = p[index];
      
      //fprintf(stderr, "%lg  %lg  %lg  %lg\n", 
      //	      k[index]/hsmall,
      //	      wk/hsmall, 
      //	      k[index+1]/hsmall,
      //	      pr * pow3(hsmall));
      
      return pr;
    }

}


//=======================
// get_index
//=======================
int get_index(double wk, double *k, int n, double hsmall)
{

  // Looks for index (left) corresponding to wk in the array k
  // of size n.
  //
  // Returns:  
  //        - the index.
  //

  int i;
  int il, ir;  // i left, i right
  int count;   // in case does not converge.
  
  il = 0;
  ir = n-1;
  i = n/2;
  count = 0;
  while( ! (wk>=k[i] && wk<=k[i+1]) )
    {
      //if(count>10)
      //fprintf(stderr, "%d  %d  %d  %lg  %lg  %lg\n", 
      //il, i, ir,
      //	k[i]/hsmall, wk/hsmall, k[i+1]/hsmall);
      if(k[i]>=wk)
	ir = i;
      else
	{
	  if(k[i]<=wk)
	    il = i;
	}
      i = il + (ir-il) / 2;
      //if(count>10)
      //fprintf(stderr, "           %d  %d  %d  %lg  %lg  %lg\n", 
      //	il, i, ir,
      //	k[i]/hsmall, wk/hsmall, k[i+1]/hsmall);
      count ++;
      if(count>1000)
	{
	  fprintf(stderr, "get_index does not converge. Aborting.\n");
	  exit(0);
	}
    } 

  ///fprintf(stderr, "====================================\n");

  return i;


}
