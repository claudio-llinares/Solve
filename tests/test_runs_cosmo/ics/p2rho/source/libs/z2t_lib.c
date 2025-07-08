#include <stdio.h>
#include<math.h>
#include<stdlib.h>

#define FUNC(x)     ((*func)(x))

double calct(double a);
double dtda(double a);
double INTEGRATE(double (*func)(double), double a, double b, double step, double eps);
void   RUNGE5VAR (double *y, double dydx, double *x, double htry, double eps, 
                  double yscale, double *dxnext, double (*func)(double));
void   RUNGE     (double y, double dydx, double x, double h,double(*func)(double),
                  double *yout, double *yerr);


#define FUNC(x)     ((*func)(x))
#define MAXSTEPS     1000000000
#define MAX(A,B) ((A)>(B)?(A):(B))
#define pow2(x) ((x)*(x))

double omega0;
double lambda0;

//double tofz(double omega0, double lambda0, double z)
double z2t(double omega0, double lambda0, double H0, double z)
{
    
    // Calculates time from z and cosmology.
    /// Arguments:
    //            omega0
    //            lambda0
    //            H0 in km/sec Mpc
    //            z
    // Output: time from a=0 to z=datum in Gyr.
   
    double a;

    //double km2mpc  = 1.0 / 3.085678e+19;
    double mpc2km  = 3.085678e+19;
    double sec2gyr = 1.0 / 3.15576000e+16;

    a = 1.0 / (1.0 + z);

    //printf("a = %lf\n", a);
    //printf("t = %lf\n",calct(a)*mpc2km/H0*sec2gyr);

    return calct(a)*mpc2km/H0*sec2gyr;

}

/*==========================================================================
 * calc_t: calculate t(a) via integration
 *==========================================================================*/
double calct(double a)
{
   double  t;
   double  dtda(double);
   
   t = (double) INTEGRATE(dtda, (double)0.0, a, (double)0.1*a, (double)1.e-5);
   
   return(t);
}

/*===========================================================================
 * da/dt  as given by Friedman equation
 c
 c     Find dt/da given a, from the Friedmann Equation.  This is exact for
 c     any isotropic-metric cosmology consistent with General Relativity.
 c     Here, "t" is understood to be in units of the inverse Hubble constant
 c     (i.e. "t" = H0*t).
 c
 c     Definitions for parameters are as in Peebles 1993, eqn (5.53).
 c
 *===========================================================================*/
double dtda(double a)
{
   if(a > 1.0e-7)
      return(1.0/sqrt(lambda0*(pow2(a)-1.) + omega0/a - omega0 + 1.));
   else
      return(0.0);
}


/*===============================================================================
* INTEGRATE:  integrates a supplied function func(x) from a to b accurately !
*===============================================================================*/


double INTEGRATE (double (*func)(double), double a, double b, double step,
                  double eps)
{
   double x, dx, y, dydx, yscale, dxnext;
   int    Nstep;
   
   
   x     = a;
   dx    = step;
   y     = (double)0.0;
   Nstep = 0;
   
   while( ((x-b)*(b-a) < (double)0.0) && (Nstep < MAXSTEPS))
     {
      Nstep  = Nstep + 1;
      dydx   = FUNC(x);
      yscale = MAX(fabs(y) + fabs(dx*dydx), (double)1.e-8);
      
      if ((x+dx-b)*(x+dx-a) > (double)0.0)
        {
         dx = b - x;
        }
      
      RUNGE5VAR(&y,dydx,&x,dx,eps,yscale,&dxnext,func);
      
      dx = dxnext;
     }
   
   return(y);
}

#define safety  (double) 0.9e0
#define pgrow   (double) -0.2e0
#define pshrink (double) -0.25e0
#define errcon  (double)  1.89e-4

void RUNGE5VAR(double *y, double dydx, double *x, double htry, double eps, 
               double yscale, double *hnext, double (*func)(double))
{
   double errmax,h,hold,htemp,xnew,yerr,ytemp,temp;
   
   h      = htry;
   errmax = (double) 10.0;
   while(errmax > (double) 1.0)
     {
      RUNGE(*y,dydx,*x,h,func,&ytemp,&yerr);
      
      errmax = fabs(yerr/yscale)/eps;
      if(errmax > (double) 1.0)
        { 
         htemp = safety*h * pow(errmax,pshrink);
         hold  = h;
         temp  = MAX(fabs(htemp),(double)0.1*fabs(h));
         
         if(h  >= 0)
            h = fabs(temp);
         else
            h = -fabs(temp);
         
         xnew = *x + h;
         
         if (fabs(xnew-*x) < (double)1e-5)
           {
            /*    	      printf("WARNING: Stepsize underflow in RUNGE5VAR()\n"); */
            h      = hold;
            errmax = (double)0.0;
           }
        }
     }
   if (errmax > errcon)
     {
      *hnext = safety*h * pow(errmax,pgrow);
     }
   else
     {
      *hnext = (double) 5.0 * h;
     }
   
   *x = *x + h;
   *y = ytemp;
   
   return;
}


#define a2      ((double) 0.2e0)
#define a3      ((double) 0.3e0)
#define a4      ((double) 0.6e0)
#define a5      ((double) 1.e0)
#define a6      ((double) 0.875e0)
#define c1      ((double) 37.e0/378.e0)
#define c3      ((double) 250.e0/621.e0)
#define c4      ((double) 125.e0/594.e0)
#define c6      ((double) 512.e0/1771.e0)
#define dc1     (c1 -  (double) 2825.e0/27648.e0)
#define dc3     (c3 -  (double) 18575.e0/48384.e0)
#define dc4     (c4 -  (double) 13525.e0/55296.e0)
#define dc5     ((double) -277.e0/14336.e0)
#define dc6     (c6 -  (double) 0.25e0)

void RUNGE(double y, double dydx, double x, double h, double(*func)(double),
           double *yout, double *yerr)
{
   double ak3, ak4, ak5, ak6;
   
   ak3  = FUNC((double)(x+a3*h));
   ak4  = FUNC((double)(x+a4*h));
   ak5  = FUNC((double)(x+a5*h));
   ak6  = FUNC((double)(x+a6*h));
   
   *yout = y + h*(c1*dydx + c3*ak3 + c4*ak4  + c6*ak6);
   *yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6);
}


