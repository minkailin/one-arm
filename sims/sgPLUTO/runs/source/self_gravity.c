/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief initialization for self-gravity [assume spherical]
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
//#include "mkl.h"
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>

#define   Re     0
#define   Im     1

#define   LTEST  6
#define   MTEST  0
#define   SIG    1.0

static double AnalyticDensity(const double r, const double theta, const double phi);
static double AnalyticPhilm(const double r, const double r1, const double r2);
static double AnalyticPotential(const double r, const double theta, const double phi, const double r1, const double r2);

int LMAX, MMAX;
int nrad, nrhs, ngst;

double sg_switchon;
double *work_real, *work_imag;
double *rho_lm_real_global, *rho_lm_imag_global;
double *raxis;

gsl_vector *rhs_real;
gsl_vector *rhs_imag;
gsl_vector *diag;
gsl_vector *upper;
gsl_vector *lower;
gsl_vector *re_sol;
gsl_vector *im_sol;

static double IntegrateChi(const double K, const double beta0, const double xi_final);
static double Getbeta0(const double K, const double nmax);

struct getbeta0_params
{
         double K, n;
};
     
double getbeta0_F (double x, void *params);
double getbeta0_dF (double x, void *params);
void   getbeta0_FdF (double x, void *params, 
                         double *y, double *dy);

/* ********************************************************************* */
void StartupSG (Data *d, Grid *G)
/*! 
 *
 *
 *
 *
 *********************************************************************** */
{
  int i, j, k, l, m, n, counter, nstar, nprev, two_n;
  int Nlm;
  double x2, x3, cos_theta, m_phi, test;
  double x, xl, xr, xp1, xm1;
  struct GRID *GX, *GY, *GZ;

  GX = G;
  GY = G + 1;
  GZ = G + 2;

  print1 ("> Initializing self-gravity ...\n");
//  TotalMass(d, G);

  LMAX = (int) g_inputParam[Lmax];
  MMAX = (int) g_inputParam[Mmax];
  sg_switchon = g_inputParam[sg_on]*2.0*CONST_PI;

  /*
    ---allocate memory for the Ylm array
    Plm is array[theta, (l,m)]. the l,m part is stored as 1D. only l+m
    = even are stored, and only l,m >= 0 stored, and only l<=LMAX,
    m<=MMAX. there are Nlm number of such elements. 
    in the (l,m) plane these elements form the locus of points
    ******************************
    m = l, for l=0....LMAX
    m = l-2, for l=2...LMAX
    .
    .
    .
    m = l-2n, for l=2n...LMAX
    .
    .
    .
    m = l-LMAX, for l=LMAX (i.e. n=LMAX/2)
    ******************************
    But m<= MMAX so not every line accepts the fill l range
    for n < nstar, some higher l modes are disgarded (because they
    imply m > MMAX)
    
    so there are nstar lines with cut-off, each has MMAX + 1 elements
    for n>=nstar is no cut-off, full range of l applies
  UPDATE: now store the Ylm's. insert trig factor

  */
  nstar = (LMAX - MMAX)/2;
  Nlm = nstar*(MMAX + 1) + (LMAX/2 + 1 - nstar)*(LMAX/2 + 1 - nstar);
  d->Ylm = ARRAY_4D(2, NX3_TOT, NX2_TOT, Nlm, double);
 
  /* --------------------------------------------------------------
   Fill Ylm with spherical harmonics. these are used for
   constructing potential from harmonic coefficients
   -------------------------------------------------------------- */
  double cosmphi, sinmphi, phi, Plm;
  KTOT_LOOP(k){ 
  phi = GZ->x[k];
  JTOT_LOOP(j) { 
    x2 = GY->x[j];
    cos_theta = cos(x2);
    
    for(n=0; n<=nstar-1; n++){/*go through each line with cut-off*/
      nprev = n*(MMAX+1);
      two_n = 2*n;
      
      for(l=two_n; l<=MMAX+two_n; l++){   
	m = l - two_n;
	counter = nprev + m;
	
        cosmphi = cos(phi*(double)m);
        sinmphi = sin(phi*(double)m);

        Plm = gsl_sf_legendre_sphPlm(l, m, cos_theta);

	d->Ylm[Re][k][j][counter]  = Plm*cosmphi;
        d->Ylm[Im][k][j][counter]  = Plm*sinmphi;
      }
    }
    
    for(n=nstar; n<=LMAX/2; n++){/*go through each line without cut-off*/
      nprev = nstar*(MMAX+1) + (n-nstar)*(LMAX + 2 - nstar - n);
      two_n = 2*n;
      
      for(l=two_n; l<=LMAX; l++){   
	m = l - two_n;
	counter = nprev + m;
	

        cosmphi = cos(phi*(double)m);
        sinmphi = sin(phi*(double)m);
        
        Plm = gsl_sf_legendre_sphPlm(l, m, cos_theta);

        d->Ylm[Re][k][j][counter]  = Plm*cosmphi;
        d->Ylm[Im][k][j][counter]  = Plm*sinmphi;

      }
    }
  }
 }
  
  /*
    Allocate array for weights in the phi integration
    for each m. 
    quad interpolation for now. 
    integration across the jth cell (from j-1/2 to j+1/2) is done by 
    cubic interpolation of the integrand across (j-1, j, j+1). 
    so need one ghost on either side
    for phi integration, integrand includes an exp(-i*m*phi)
    factor. this factor is integrated by hand. 
  */
  double dphi, th, costh, sincth, ReI0, ImI0, ReI1, ImI1, ReI2, ImI2, dblem, Re0, Im0, Re1, th2;
  double filon_corr; 
  dphi = GZ->dx[KBEG]; /*assume uniform phi spacing*/
  
  d->phi_weights_m = ARRAY_3D(2, MMAX+1, NX3_TOT, double);
  
  for(m=0;m<=MMAX;m++){ /*zero the weights array*/
    KTOT_LOOP(k){
      d->phi_weights_m[Re][m][k] = 0.0;
      d->phi_weights_m[Im][m][k] = 0.0;
    }
  }
  
  m = 0 ; /*axisymmetric mode, no trig factor, only weights from quad
	    interpolation. for uniform spacing these expressions
	    actually simplify a lot w_(j-1) is dphi/24, w_j is
	    11.dphi/12.0 and w_j(j+1) = dphi/24. 
	  */
  KDOM_LOOP(k){
    d->phi_weights_m[Re][m][k-1] += dphi/24.0;
    d->phi_weights_m[Re][m][k]   += 11.0*dphi/12.0;
    d->phi_weights_m[Re][m][k+1] += dphi/24.0;

    /* below is for filon-trapezoidal rule */
    /*
      d->phi_weights_m[Re][m][k] = dphi;
      d->phi_weights_m[Im][m][k] = 0.0;
    */
  }

  for(m=1;m<=MMAX;m++){/*non-axisymmetric modes, there are terms
			 accounting for the trig factor. see notes*/
    KDOM_LOOP(k){
      x   = GZ->x[k];
      xp1 = GZ->x[k+1];
      xm1 = GZ->x[k-1];

      dblem = (double)m;
      m_phi = x*dblem;
      cosmphi = cos(m_phi);
      sinmphi = sin(m_phi);

      th = dphi*dblem/2.0;
      costh = cos(th);
      sincth= sin(th)/th;

/*
      ReI0 = cosmphi*dphi*sincth;
      ImI0 =-sinmphi*dphi*sincth;

      ReI1 = cosmphi*x*sincth + sinmphi*(1.0/dblem)*(costh - sincth);
      ReI1*= dphi;
      ImI1 = cosmphi*(1.0/dblem)*(costh - sincth) - sinmphi*x*sincth;
      ImI1*= dphi;

      ReI2 = cosmphi*( 2.0*costh/(dblem*dblem) + sincth*(
							 (th*th-2.0)/(dblem*dblem)
							 + x*x) )
	+ sinmphi*2.0*x*(1.0/dblem)*(costh - sincth);
      ReI2*= dphi;
      ImI2 = cosmphi*2.0*(x/dblem)*(costh - sincth) 
	- sinmphi*(2.0*costh/(dblem*dblem) + sincth*(
						     (th*th-2.0)/(dblem*dblem) + x*x)); 
      ImI2*= dphi;

      d->phi_weights_m[Re][m][k-1]
      	+=( x*xp1*ReI0 - (x+xp1)*ReI1 + ReI2 )/(2.0*dphi*dphi);
      d->phi_weights_m[Re][m][k]
      	+=( (xm1+xp1)*ReI1 - xm1*xp1*ReI0 - ReI2)/(dphi*dphi);
      d->phi_weights_m[Re][m][k+1]
      	+=( ReI2 - (x+xm1)*ReI1 + x*xm1*ReI0)/(2.0*dphi*dphi);

      
      d->phi_weights_m[Im][m][k-1]
      	+=( x*xp1*ImI0 - (x+xp1)*ImI1 + ImI2 )/(2.0*dphi*dphi);
      d->phi_weights_m[Im][m][k]
      	+=( (xm1+xp1)*ImI1 - xm1*xp1*ImI0 - ImI2)/(dphi*dphi);
      d->phi_weights_m[Im][m][k+1]
      	+=( ImI2 - (x+xm1)*ImI1 + x*xm1*ImI0)/(2.0*dphi*dphi);
*/


     th2 = th*th;
     Re0 = 2.0*costh/th2 + (th2 - 2.0)*sincth/th2;
     Im0 = (2.0/th)*(sincth - costh);
     Re1 = costh/th2 - 0.5*(3.0*th2+2.0)*sincth/th2;

      d->phi_weights_m[Re][m][k-1]
        +=(dphi/8.0)*(cosmphi*Re0 + sinmphi*Im0);
      d->phi_weights_m[Re][m][k]
        +=-(dphi/2.0)*cosmphi*Re1;
      d->phi_weights_m[Re][m][k+1]
        +=(dphi/8.0)*(cosmphi*Re0 - sinmphi*Im0);

      d->phi_weights_m[Im][m][k-1]
        +=(dphi/8.0)*(cosmphi*Im0 - sinmphi*Re0);
      d->phi_weights_m[Im][m][k]
        +=(dphi/2.0)*sinmphi*Re1;
      d->phi_weights_m[Im][m][k+1]
        +=(dphi/8.0)*(-cosmphi*Im0 - sinmphi*Re0);
 

      /* below is for filon-trapezoidal rule*/     
/*
       filon_corr = sincth*sincth; 
       d->phi_weights_m[Re][m][k] = cosmphi*filon_corr*dphi; 
       d->phi_weights_m[Im][m][k] =-sinmphi*filon_corr*dphi; 
*/
    }
  }

  /*
    Allocate array for weights in the theta integration
    for each (l,m). 
    quad interpolation for now. 
    integration across the jth cell (from j-1/2 to j+1/2) is done by 
    quad interpolation of the integrand across (j-1, j, j+1).
    so need one ghost on either side
    the weights simplify for uniform spacing but here it is for general spacing
  */

  double dx, p, pm1, pp1;

  d->theta_weights_lm = ARRAY_2D(Nlm, NX2_TOT, double);
  
  for(n=0; n<=nstar-1; n++){/*go through each line with cut-off*/
    nprev = n*(MMAX+1);
    two_n = 2*n;
    
    for(l=two_n; l<=MMAX+two_n; l++){
      m = l - two_n;
      counter = nprev + m;
      
      JTOT_LOOP(j){/*zero the weights array*/
	d->theta_weights_lm[counter][j] = 0.0;
      }
      
      JDOM_LOOP(j){/*go through each active zone, updating the weights
		     for j-1, j and j+1 */
	x   = GY->x[j];
	xm1 = GY->x[j-1];
	xp1 = GY->x[j+1];
	xl  = GY->xl[j];
	xr  = GY->xr[j];
	dx  = GY->dx[j];
	
	p = gsl_sf_legendre_sphPlm(l, m, cos(x))*sin(x);
	pm1 = gsl_sf_legendre_sphPlm(l, m, cos(xm1))*sin(xm1);
	pp1 = gsl_sf_legendre_sphPlm(l, m, cos(xp1))*sin(xp1);

	/*below: for cos theta = x formulation, need a minus sign at
	  the front, but note that density is not at center of such
	  grid cell*/
	
	/* x = cos(x); */
	/* xm1=cos(xm1); */
	/* xp1=cos(xp1); */
	/* xl= cos(xl); */
	/* xr= cos(xr); */
	/* dx= xr - xl; */

	/* p = -gsl_sf_legendre_sphPlm(l, m, x); */
	/* pm1 = -gsl_sf_legendre_sphPlm(l, m, xm1); */
	/* pp1 = -gsl_sf_legendre_sphPlm(l, m, xp1); */
	/*above: for cos theta = x formulation*/


	
	d->theta_weights_lm[counter][j-1] 
           += pm1*(x*xp1*dx - 0.5*(x+xp1)*(xr*xr-xl*xl)+(1.0/3.0)*(xr*xr*xr - xl*xl*xl))
                   /( (x-xm1)*(xp1-xm1) );

	d->theta_weights_lm[counter][j] 
           +=p*(0.5*(xp1+xm1)*(xr*xr-xl*xl) - (1.0/3.0)*(xr*xr*xr - xl*xl*xl) - xm1*xp1*dx )
                /( (x-xm1)*(xp1-x) );

	d->theta_weights_lm[counter][j+1] 
           +=pp1*((1.0/3.0)*(xr*xr*xr -  xl*xl*xl) - 0.5*(x+xm1)*(xr*xr - xl*xl) + x*xm1*dx)
                  /( (xp1-x)*(xp1-xm1) );

      }
    } 
  }
  
  for(n=nstar; n<=LMAX/2; n++){/*go through each line without cut-off*/
    nprev = nstar*(MMAX+1) + (n-nstar)*(LMAX + 2 - nstar - n);
    two_n = 2*n;
    
    for(l=two_n; l<=LMAX; l++){
      m = l - two_n;
      counter = nprev + m;
      
      JTOT_LOOP(j){
	d->theta_weights_lm[counter][j] = 0.0;
      }
      
      JDOM_LOOP(j){
	x = GY->x[j];
	xm1 = GY->x[j-1];
	xp1 = GY->x[j+1];
	xl= GY->xl[j];
	xr= GY->xr[j];
	dx= GY->dx[j];
	
	p = gsl_sf_legendre_sphPlm(l, m, cos(x))*sin(x);
	pm1 = gsl_sf_legendre_sphPlm(l, m, cos(xm1))*sin(xm1);
	pp1 = gsl_sf_legendre_sphPlm(l, m, cos(xp1))*sin(xp1);
	
	/*below: for cos theta = x formulation, need a minus sign at
	  the front, but note that density is not at center of such
	  grid cell*/
	/* x = cos(x); */
	/* xm1=cos(xm1); */
	/* xp1=cos(xp1); */
	/* xl= cos(xl); */
	/* xr= cos(xr); */
	/* dx= xr-xl; */

	/* p = -gsl_sf_legendre_sphPlm(l, m, x); */
	/* pm1 = -gsl_sf_legendre_sphPlm(l, m, xm1); */
	/* pp1 = -gsl_sf_legendre_sphPlm(l, m, xp1); */
	/*above: for cos theta = x formulation*/


	d->theta_weights_lm[counter][j-1]
           += pm1*(x*xp1*dx -  0.5*(x+xp1)*(xr*xr-xl*xl)+(1.0/3.0)*(xr*xr*xr  - xl*xl*xl))
                   /( (x-xm1)*(xp1-xm1) );
	
	d->theta_weights_lm[counter][j]
           +=p*(0.5*(xp1+xm1)*(xr*xr-xl*xl) - (1.0/3.0)*(xr*xr*xr - xl*xl*xl) - xm1*xp1*dx )
                /( (x-xm1)*(xp1-x) );

	d->theta_weights_lm[counter][j+1]
           +=pp1*((1.0/3.0)*(xr*xr*xr - xl*xl*xl) - 0.5*(x+xm1)*(xr*xr - xl*xl) + x*xm1*dx)
                  /( (xp1-x)*(xp1-xm1) );

	/*if uniform spacing...*/
	/* d->theta_weights_lm[counter][j-1] += pm1*dx/24.0; */
	/* d->theta_weights_lm[counter][j]   += 11.0*p*dx/12.0; */
	/* d->theta_weights_lm[counter][j+1] += pp1*dx/24.0; */

      }
    }      
  }
  
  /*
    ---allocate memory for rho_m(r, theta) =
    integral[rho*exp(i*m*phi)] dphi. 
    we need one ghost on either side of the theta range for
    theta-integration step after phi-integration is done.
    first index is re/im part
  */ 
  d->rho_m = ARRAY_4D(2, MMAX+1, NX2_TOT, NX1_TOT, double);

 /*
 *     ---allocate memory for
 *     local copy of the Phi_lm(r) array. need one ghost on either side of the radial grid
 *                       */

  d->Philm = ARRAY_3D(2, Nlm, NX1_TOT, double);

/*
 *     ---allocate memory for
 *     local copy of the Phi array. fill all active zones plus one
 *     ghost on all boundaries. initialize to zero in case we're doing
 *     non-sg run
 */
  d->Phi = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  TOT_LOOP(k,j,i){
    d->Phi[k][j][i] = 0.0;
  }

  /*
    ---allocate memory for work arrays when solving 1D poisson
  */ 
  nrad = GX->np_tot_glob;
  nrhs = GX->np_int_glob + 2;
  ngst = GX->nghost;

  work_real = ARRAY_1D(nrad, double);
  work_imag = ARRAY_1D(nrad, double);  
  rho_lm_real_global = ARRAY_1D(nrad, double);
  rho_lm_imag_global = ARRAY_1D(nrad, double);
  raxis    =  ARRAY_1D(nrhs, double);
  for(i=0; i<nrhs; i++){/*initialize radial grid*/
    raxis[i]    = GX->x_glob[ngst-1+i];
  }
  
  rhs_real  = gsl_vector_alloc(nrhs);
  rhs_imag  = gsl_vector_alloc(nrhs);
  diag  = gsl_vector_alloc(nrhs);
  upper = gsl_vector_alloc(nrhs-1);
  lower = gsl_vector_alloc(nrhs-1);
  re_sol = gsl_vector_alloc(nrhs);
  im_sol = gsl_vector_alloc(nrhs);


 d->beta2d = ARRAY_2D(NX2_TOT, NX1_TOT, double); //2d array to store initial density enhancement due to SG
  TOT_LOOP(k,j,i){
    d->beta2d[j][i] = 1.0;
  }
}


void GetRhom (Data *d, Grid *G)
/*! 
 * Given the density cube rho(r,theta,phi) perform the integral
 *                      rho.exp(-i.m.phi)dphi
 *over local active phi zones
 *for all active radial zones
 *for all active theta-zones PLUS one ghost on either side
 *currently use quad interpolation so need one phi ghost on either
 *side (loop is from KBEG-1 to KEND+1)
 *if reverting to filon-trapezoidal, then don't need ghosts (loop over
 *active zones only)
 *********************************************************************** */
{
  int i, j, k, m;
  struct GRID *GX, *GY, *GZ;
 
  GX = G;
  GY = G + 1; 
  GZ = G + 2;

  /*axisymmetric mode*/ 
  m = 0;
  k = KBEG-1; /*initialize rho_m(r,theta) with contributions from phi=phi_start*/
  for(j=JBEG-1; j<=JEND+1; j++){
    IDOM_LOOP(i){
      d->rho_m[Re][m][j][i] = d->Vc[RHO][k][j][i]*d->phi_weights_m[Re][m][k];
      d->rho_m[Im][m][j][i] = 0.0;
    }
  }
  for(k=KBEG;k<=KEND+1;k++){/*add contributions from other (r,theta) slices. no imaginary part from m=0*/
    for(j=JBEG-1; j<=JEND+1; j++){
      IDOM_LOOP(i){
	d->rho_m[Re][m][j][i] += d->Vc[RHO][k][j][i]*d->phi_weights_m[Re][m][k];
      }
    }
  }
  /*non-axisymmetric modes*/
  for(m=1; m<=MMAX; m++){
    k=KBEG-1;
    for(j=JBEG-1; j<=JEND+1; j++){
      IDOM_LOOP(i){
  	d->rho_m[Re][m][j][i] = d->Vc[RHO][k][j][i]*d->phi_weights_m[Re][m][k];
  	d->rho_m[Im][m][j][i] = d->Vc[RHO][k][j][i]*d->phi_weights_m[Im][m][k];
      }
    }
    for(k=KBEG;k<=KEND+1;k++){
      for(j=JBEG-1; j<=JEND+1; j++){
  	IDOM_LOOP(i){
  	  d->rho_m[Re][m][j][i] += d->Vc[RHO][k][j][i]*d->phi_weights_m[Re][m][k];
  	  d->rho_m[Im][m][j][i] += d->Vc[RHO][k][j][i]*d->phi_weights_m[Im][m][k];
 	}
      }
    }
  }
}

int GetCounter(const int l, const int m)
{
  int n, nstar, nprev;
  /* Given (l,m), work out the position of storage in 1D array of
     spherical harmonics
   */

  nstar = (LMAX - MMAX)/2;
  n = (l-m)/2;
  
  if(n<=nstar-1){/*we're on a line with cut off*/
    nprev = n*(MMAX+1);
  } else {/*we're on a line without cutoff*/
    nprev = nstar*(MMAX+1) + (n-nstar)*(LMAX + 2 - nstar - n);
  }
  return nprev + m;
}

void GetPhilm (Data *d, Grid *G, const int l, const int m)
/*! 
 * given the l,m (corresonding to a unique counter in the 1D storage of the spherical harmonics),
 * solve the 1D poisson equation to get the local portion of phi_lm(r) --> all active plus one ghost on either side
 *
 * step1 rho_m --> rho_lm(r) via theta-integration 
 * step2 copy rho_lm(r) into appropriate portion of rho_lm_global(r) (which covers the global grid plus all ghosts)
 * step3 MPI reduce rho_lm_global(r) -->  but only the active zones are non zero 
 * step4 using rho_lm_global to get rhs of poisson --> active zones plus one ghost on either side
 * step5 solve 1D poisson ---> gives phi_lm(r) of the global grid. (size active plus only 1 ghost on either side). 
 * step6 phi_lm_local ---> only fill active zones plus one ghost on either size 
 ************************************************************************ */
{
  int i,j, counter;
  double rho_lm_real[NX1_TOT], rho_lm_imag[NX1_TOT];
  struct GRID *GY, *GX, *GZ;
  
  GZ = G+2;
  counter = GetCounter(l, m);

  /* Step 1: theta-integration (requires one ghost on either side)
********************************/
  GY = G + 1;

  j=JBEG-1;/*initialize for all r*/
  IDOM_LOOP(i){
    rho_lm_real[i] =
      d->rho_m[Re][m][j][i]*d->theta_weights_lm[counter][j];
    rho_lm_imag[i] =
      d->rho_m[Im][m][j][i]*d->theta_weights_lm[counter][j];
  }
  
  for(j=JBEG;j<=JEND+1;j++){/*add contribution from other theta*/
    IDOM_LOOP(i){
      rho_lm_real[i] +=
      	d->rho_m[Re][m][j][i]*d->theta_weights_lm[counter][j];
      rho_lm_imag[i] +=
      	d->rho_m[Im][m][j][i]*d->theta_weights_lm[counter][j];
    }
  }

  /* Step 2-3: add up contributions to rho_lm from all procs
********************************/
  int rcoord;
  
  GX = G;
  rcoord = GX->beg; /*start index in global radial grid*/
   
  for(i=0; i<nrad; i++){/*initialize*/
    work_real[i] = 0.0;
    work_imag[i] = 0.0;

    rho_lm_real_global[i] = 0.0;
    rho_lm_imag_global[i] = 0.0;
  }
  
  IDOM_LOOP(i){/*double to account for lower half disk (z<0)*/
    work_real[rcoord + i - IBEG] = rho_lm_real[i]*2.0;
    work_imag[rcoord + i - IBEG] = rho_lm_imag[i]*2.0;
  }
  MPI_Barrier (MPI_COMM_WORLD);

  MPI_Allreduce (&work_real[0], &rho_lm_real_global[0], nrad, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  MPI_Allreduce (&work_imag[0], &rho_lm_imag_global[0], nrad, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   
  /* Step 4A: construct rhs for 1D poisson
     rho_lm now covers the entire global radial grid, including
     all ghosts on either side. 
     poisson rhs is all the active zones plus ONE ghost on either side
  ********************************/
  double re, im;
 
  for(i=0; i<nrhs; i++){/*initialize rhs grid*/
    gsl_vector_set(rhs_real, i,   0.0);
    gsl_vector_set(rhs_imag, i,   0.0);
  }
  
  for(i=1; i<=nrhs-2; i++){/*fill the active zones of rhs with
			     rho_lm. assume G=1. */
    re = 4.0*CONST_PI*rho_lm_real_global[ngst+i-1];
    im = 4.0*CONST_PI*rho_lm_imag_global[ngst+i-1];

    gsl_vector_set(rhs_real, i,   re);
    gsl_vector_set(rhs_imag, i,   im);
  }

   /* STEP  4B: construct matrix for 1D poisson. get
      diagonal, with nrhs elements
      lower diagonal, with nrhs-1 elements
      upper diagonal, with nrhs-1 elements
      in future, should store some of these in memory
********************************/
   double r1, r2, idx, x, xp1, xm1, dd, ud, ld;
   double dblel=(double)l; 

   r1 = GX->xl_glob[GX->gbeg];
   r2 = GX->xr_glob[GX->gend];
  
   /* /\*inner boundary condition*\/ */
   i=0;
   dd = -( 1.0/(raxis[i+1] - raxis[i]) + 0.5*dblel/r1);
   ud = 1.0/(raxis[i+1] - raxis[i]) - 0.5*dblel/r1;

   gsl_vector_set(diag, i, dd);
   gsl_vector_set(upper,i, ud);
   /* /\*outer boundary condition*\/ */
   i=nrhs-1;
   dd = 1.0/(raxis[i] - raxis[i-1]) + 0.5*(dblel+1.0)/r2;
   ld = -1.0/(raxis[i] - raxis[i-1]) + 0.5*(dblel+1.0)/r2;
   
   gsl_vector_set(diag, i,   dd);
   gsl_vector_set(lower,i-1, ld);
   /* /\*interior active zones*\/ */
   for(i=1; i<=nrhs-2; i++){
     idx = 1.0/GX->dx_glob[ngst-1+i];
     x  = raxis[i];
     xm1= raxis[i-1];
     xp1= raxis[i+1];

     dd  = -( 1.0/(xp1-x) + 1.0/(x-xm1) )*idx;
     dd += -dblel*(dblel + 1.0)/(x*x);
     
     ud  = idx/(xp1-x);
     ud += (2.0/x)/(xp1-xm1);

     ld  = idx/(x-xm1);
     ld +=-(2.0/x)/(xp1-xm1);

     gsl_vector_set(diag, i, dd);
     gsl_vector_set(upper,i, ud);
     gsl_vector_set(lower,i-1, ld);
   }

   /* STEP  4C: solve tri-diagonal system
********************************/
   /*real part*/
   gsl_linalg_solve_tridiag(diag, upper, lower, rhs_real, re_sol); 
   /*imag part*/
   gsl_linalg_solve_tridiag(diag, upper, lower, rhs_imag, im_sol);


   /* STEP  5: extract relevant portion of global solution phi_lm to
      local storage.
      local copy of phi_lm has size NX1_TOT (all active plus ghosts)
      but only active plus one ghost either side is filled
 *********************************/
   int pos;
   for(i=IBEG-1; i<=IEND+1; i++){
     pos = rcoord - (ngst-1) - 1  + i -(IBEG-1);
     d->Philm[Re][counter][i] = gsl_vector_get(re_sol, pos );  
     d->Philm[Im][counter][i] = gsl_vector_get(im_sol, pos );
   }
      
}

void GetPhi(Data *d, Grid *G)
/*!
  Direct sum to get potential in each cell, self-potential then indirect potential
  ************************************************************************ */
{
  int i,j,k,n,nstar,two_n,l,nprev,m,counter;
  
  nstar = (LMAX - MMAX)/2;
  
  for(k=KBEG-1; k<=KEND+1; k++){
    for(j=JBEG-1; j<=JEND+1; j++){
     

      ITOT_LOOP(i){/*initialize potential to zero*/
        d->Phi[k][j][i] = 0.0;
      }

      for(n=0; n<=nstar-1; n++){/*go through each line with cut-off */
	nprev = n*(MMAX+1);
	two_n = 2*n;
	
	for(i=IBEG-1; i<=IEND+1; i++){
	  d->Phi[k][j][i] +=
	    d->Philm[Re][nprev][i]*d->Ylm[Re][k][j][nprev];/*contribution
							     from m=0, which is the first element of this line*/ 
	}
	
	for(l=two_n+1; l<=MMAX+two_n; l++){/*contribution from rest of
					     line (nonaxisymmetric modes)*/
	  m = l - two_n;
	  counter = nprev + m;
	  
	  for(i=IBEG-1; i<=IEND+1; i++){
	    d->Phi[k][j][i]+=
	      2.0*(d->Philm[Re][counter][i]*d->Ylm[Re][k][j][counter] -
		   d->Philm[Im][counter][i]*d->Ylm[Im][k][j][counter]);
	      }
	}       
      }
     
      for(n=nstar; n<=LMAX/2; n++){/*go through each line without cut-off*/
	nprev = nstar*(MMAX+1) + (n-nstar)*(LMAX + 2 - nstar - n);
	two_n = 2*n;
	
	for(i=IBEG-1; i<=IEND+1; i++){
	  d->Phi[k][j][i] +=
	    d->Philm[Re][nprev][i]*d->Ylm[Re][k][j][nprev];/*contribution
							     from m=0, which is the first element of this line*/ 
	}
	
	for(l=two_n+1; l<=LMAX; l++){
	  m = l - two_n;
	  counter = nprev + m;
	  
	  for(i=IBEG-1; i<=IEND+1; i++){
	    d->Phi[k][j][i]+=
	      2.0*(d->Philm[Re][counter][i]*d->Ylm[Re][k][j][counter] -
		   d->Philm[Im][counter][i]*d->Ylm[Im][k][j][counter]);
	      } 
	}
      }
      
    }
  }

//indirect potential. assume equitorial symm so no z-contrib. assume G=1
 double I_loc[2], I_glob[2];  
 double dm, dvol, cosphi, sinphi, sintheta_cosphi, sintheta_sinphi;
 struct GRID *GX, *GY, *GZ;

  GX = G;
  GY = G + 1;
  GZ = G + 2;
 
  I_loc[0] = I_loc[1] = 0.0;
  I_glob[0] = I_glob[1] = 0.0;

/*
  KDOM_LOOP(k){
   cosphi = cos(GZ->x[k]);
   sinphi = sin(GZ->x[k]);
  JDOM_LOOP(j){
   sintheta_cosphi = sin(GY->x[j])*cosphi;
   sintheta_sinphi = sin(GY->x[j])*sinphi;
  IDOM_LOOP(i){
    dvol  = GX->dV[i]*GY->dV[j]*GZ->dV[k];
    dm    = dvol*d->Vc[RHO][k][j][i];
    dm   /= GX->x[i]*GX->x[i];

    I_loc[0] += dm*sintheta_cosphi;
    I_loc[1] += dm*sintheta_sinphi;
}
}
}
*/

//global reduction
MPI_Allreduce (&I_loc[0], &I_glob[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

//add to self-potential
double x,y,indir;
for(k=KBEG-1;k<=KEND+1;k++){
cosphi = cos(GZ->x[k]);
sinphi = sin(GZ->x[k]);
for(j=JBEG-1;j<=JEND+1;j++){
sintheta_cosphi = sin(GY->x[j])*cosphi;
sintheta_sinphi = sin(GY->x[j])*sinphi;
for(i=IBEG-1;i<=IEND+1;i++){
   x = GX->x[i]*sintheta_cosphi;
   y = GX->x[i]*sintheta_sinphi;
//factor of 2 to account for lower half disk. no z-contribu due to symm
   indir = 2.0*( x*I_glob[0] + y*I_glob[1]);

   d->Phi[k][j][i] += indir;
}
}
}





}




void SolvePoisson(Data *d, Grid *G)
/*!
  main routine to solve poission's equation via spherical harmonic
  expansion:
  step 1: azimuthal interation of density field
  step 2: solve 1D poisson for each l,m
  step 3: re-construct potential
  ************************************************************************ */
{
  /*Step 1: azimuthal integration (Fourier transform)
   */
  GetRhom (d, G);
  MPI_Barrier (MPI_COMM_WORLD);

  /*Step 2: solve 1D poisson for each l,m
   */
  int l, m, n, nstar, nprev, counter, two_n;

  nstar = nstar = (LMAX - MMAX)/2;

  for(n=0; n<=nstar-1; n++){/*go through each line with cut-off*/
    nprev = n*(MMAX+1);
    two_n = 2*n;  
    for(l=two_n; l<=MMAX+two_n; l++){   
      m = l - two_n;
      GetPhilm (d, G, l, m);
      /* MPI_Barrier (MPI_COMM_WORLD); */
    }
  }
  
  for(n=nstar; n<=LMAX/2; n++){/*go through each line without cut-off*/
    nprev = nstar*(MMAX+1) + (n-nstar)*(LMAX + 2 - nstar - n);
    two_n = 2*n;
    
    for(l=two_n; l<=LMAX; l++){   
      m = l - two_n;
      GetPhilm (d, G, l, m);
      /* MPI_Barrier (MPI_COMM_WORLD); */
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);

  /*Step 3: construct potential
   */
  GetPhi(d, G);

  if(g_inputParam[sg_correction]<0.0)  { //then we are turning on the potential slowly
  //Step 4: taper off potential
  int i,j,k;
  double taper, arg;

  if(g_time <= sg_switchon){
    arg   = 0.5*CONST_PI*g_time/sg_switchon;
    taper = pow(sin(arg), 2.0);

   TOT_LOOP(k,j,i){
   d->Phi[k][j][i] *= taper;
   }
  }
}
} 

int func (double t, const double y[], double f[],
           void *params)
     {
       double K = *(double *)params;
       f[0] = y[1];
       f[1] = -K - exp(y[0]);
       return GSL_SUCCESS;
     }
     
int  jac (double t, const double y[], double *dfdy, 
          double dfdt[], void *params)
     {
       double K = *(double *)params;
       gsl_matrix_view dfdy_mat 
         = gsl_matrix_view_array (dfdy, 2, 2);
       gsl_matrix * m = &dfdy_mat.matrix; 
       gsl_matrix_set (m, 0, 0, 0.0);
       gsl_matrix_set (m, 0, 1, 1.0);
       gsl_matrix_set (m, 1, 0, -exp(y[0]));
       gsl_matrix_set (m, 1, 1, 0.0);
       dfdt[0] = 0.0;
       dfdt[1] = 0.0;
       return GSL_SUCCESS;
     }

static double IntegrateChi(const double K, const double beta0, const double xi_final)
{
double mu, chi_final;
mu = K;
  
gsl_odeiv2_system sys = {func, jac, 2, &mu};
     
gsl_odeiv2_driver * d = 
         gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
     				  1e-9, 1e-9, 0.0);
       int i;
       double t = 0.0, t1 = xi_final;
       double y[2] = { log(beta0), 0.0 };

       for (i = 1; i <= 100; i++)
         {
           double ti = i * t1 / 100.0;
           int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
     
           if (status != GSL_SUCCESS)
     	{
     	  printf ("error, return value=%d\n", status);
     	  break;
     	}
          
        }
       chi_final = y[0];
       gsl_odeiv2_driver_free (d);
       return chi_final; //chi(xi_final)
}

double getbeta0_F (double x, void *params)
{
struct getbeta0_params *p = (struct getbeta0_params *) params;
     
       double K = p->K;
       double n = p->n;

       double Chi_final, xi_final, Chiprime_final;

       xi_final = n/sqrt(K);

       Chi_final = IntegrateChi(K, x, xi_final);
       Chiprime_final = 2.0*(x + K*log(x) -K*Chi_final - exp(Chi_final));
     
       return sqrt(2.0*K/CONST_PI)*(-sqrt(K)*n + sqrt(Chiprime_final)) - gsl_sf_erf(n/sqrt(2.0));
}
     
double getbeta0_dF (double x, void *params)
     {
       struct getbeta0_params *p 
         = (struct getbeta0_params *) params;
       double K = p->K;
       double n = p->n;
       struct getbeta0_params params1 = {K, n};       

       gsl_function F;
       F.function = &getbeta0_F;
       F.params =   &params1;

       double dx =1.0e-6;
       double f, fp;
       f =  GSL_FN_EVAL(&F, x);
       fp=  GSL_FN_EVAL(&F, x+dx); 
     
       return (fp - f)/dx;
     }
     
void getbeta0_FdF (double x, void *params, 
                    double *y, double *dy)
     {
       struct getbeta0_params *p 
         = (struct getbeta0_params *) params;
       double K = p->K;
       double n = p->n;
       struct getbeta0_params params1 = {K, n}; 
     
       double f;
       gsl_function F;
       F.function = &getbeta0_F;
       F.params = & params1;
       f = GSL_FN_EVAL(&F, x);

       double df;    
       F.function = &getbeta0_dF;
       F.params = & params1;
       df = GSL_FN_EVAL(&F, x);
  
       *y = f;
       *dy= df;
     }

static double Getbeta0(const double K, const double nmax)
{
       int status;
       int iter = 0, max_iter = 1000;
       const gsl_root_fdfsolver_type *T;
       gsl_root_fdfsolver *s;
       double x0, x = 1.3;//initial guess for midplane enhancement

       gsl_function_fdf FDF;

       struct getbeta0_params params = {K, nmax};
     
       FDF.f = &getbeta0_F;
       FDF.df = &getbeta0_dF;
       FDF.fdf = &getbeta0_FdF;
       FDF.params = &params;
     
       T = gsl_root_fdfsolver_newton;
       s = gsl_root_fdfsolver_alloc (T);
       gsl_root_fdfsolver_set (s, &FDF, x);
     
       do
         {
           iter++;
           status = gsl_root_fdfsolver_iterate (s);
           x0 = x;
           x = gsl_root_fdfsolver_root (s);
           status = gsl_root_test_delta (x, x0, 0, 1e-9);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
     
       gsl_root_fdfsolver_free (s);
       return x;
}

void SGDensity(Data *d, Grid *G)
/*!
 *   ************************************************************************ */
{
  int i,j,k;
  double K, nmax;
  double r, theta, R, z;
  double beta, beta0, xi_final, chi_final; 
  struct GRID *GY, *GX, *GZ;
  GX = G;
  GY = G + 1;
  GZ = G + 2;
  //account for vertical SG on density field
    JTOT_LOOP(j){
      theta = GY->x[j];
      ITOT_LOOP(i){
      r = GX ->x[i];

      R = r*sin(theta);
      K = csq(R)/(4.0*CONST_PI*pow(bigH(R), 2.0));
      K/= surface_density(R)/( bigH(R)*sqrt(2.0*CONST_PI) );
 
      z = r*cos(theta);    
      xi_final = (z/bigH(R))/sqrt(K);

      nmax  = g_inputParam[max_H]*g_inputParam[smallh]*R/bigH(R);
      beta0 = Getbeta0(K, nmax);

      if(xi_final < 0.0) xi_final = -xi_final;

      chi_final = IntegrateChi(K, beta0, xi_final);
       
      beta = exp(chi_final + 0.5*K*xi_final*xi_final);
      d->beta2d[j][i] = beta;

     if(g_inputParam[sg_correction]>0.0)  {   
      KTOT_LOOP(k){
         d->Vc[RHO][k][j][i] *= beta;
         #if EOS == IDEAL
         d->Vc[PRS][k][j][i] *= beta;
         #endif
      }
      }
      }
    }

    if(g_inputParam[sg_correction]>0.0)  {
    SolvePoisson(d, G);
    
    double vphisq_nsg; 
    double vphisq_sg, dbeta_dr, dbeta_dth, dPhi_dr, dPhi_dth;
    double dbeta_dR, dPhi_dR;

    DOM_LOOP(k,j,i){
    r     = GX->x[i];
    theta = GY->x[j];   
    R = r*sin(theta);

    dbeta_dr  = (log(d->beta2d[j][i+1]) - log(d->beta2d[j][i-1]))/(GX->x[i+1] - GX->x[i-1]);
    dbeta_dth = (log(d->beta2d[j+1][i]) - log(d->beta2d[j-1][i]))/(GY->x[j+1] - GY->x[j-1]);
    dbeta_dR = sin(theta)*dbeta_dr + cos(theta)*dbeta_dth/r;

    dPhi_dr  = (d->Phi[k][j][i+1] - d->Phi[k][j][i-1])/(GX->x[i+1] - GX->x[i-1]);
    dPhi_dth = (d->Phi[k][j+1][i] - d->Phi[k][j-1][i])/(GY->x[j+1] - GY->x[j-1]);
    dPhi_dR = sin(theta)*dPhi_dr + cos(theta)*dPhi_dth/r;

    vphisq_sg =  R*(dPhi_dR + csq(R)*dbeta_dR);

    vphisq_nsg = pow(d->Vc[VX3][k][j][i], 2.0);

    d->Vc[VX3][k][j][i] = sqrt(vphisq_nsg + vphisq_sg);
   }

  //outer boundary correction /linear extrapolation for dphi/dr
    double vphi_sg1, vphi_sg2;

    X1_END_LOOP(k,j,i){
        r        = GX->x[i];
        vphi_sg1 = (d->Phi[k][j][IEND] - d->Phi[k][j][IEND-1])/(GX->x[IEND] - GX->x[IEND-1]);
        vphi_sg2 = (d->Phi[k][j][IEND+1] - d->Phi[k][j][IEND])/(GX->x[IEND+1] - GX->x[IEND]);
        vphisq_sg  = vphi_sg2 - (GX->xr[IEND] - r)*(vphi_sg2 - vphi_sg1)/GX->dx[IEND];
        vphisq_sg  = r*vphisq_sg;
                                          
        vphisq_nsg = pow(d->Vc[VX3][k][j][i], 2.0);

        d->Vc[VX3][k][j][i] = sqrt(vphisq_nsg + vphisq_sg);
   }
 }
}

static double AnalyticDensity(const double r, const double theta, const double phi)
{
 /*
 * prescribe non-axisymmetric density field
 */
 double nonaxi, plm, rho0;

 nonaxi  = cos(phi*(double)MTEST);
 nonaxi += sin(phi*(double)MTEST);	
 plm     = gsl_sf_legendre_sphPlm(LTEST, MTEST, cos(theta));
 rho0    = pow(1.0/r, SIG);
 
 return rho0*nonaxi*plm;
}

static double AnalyticPhilm(const double r, const double r1, const double r2)
{
 /*
 * get the real part of phi_lm. 
 * for above density field, imag part is -real part (for non axisymmetric mode) or zero (for axisymmetric mode)  
 */	
   double l = (double)LTEST;
   double philm_real, factor, p;

   p = SIG - 2.0;
   factor = 4.0*CONST_PI;
   factor/= p*(p-1.0) - l*(l + 1.0);
   if(MTEST != 0) factor /= 2.0;

   philm_real =( (SIG - l - 3.0)/(2.0*l + 1.0) )
       *pow(1.0/r2, p)*pow(r/r2, l)
       -( (p + l)/(2.0*l + 1.0) )
       *pow(1.0/r1, p)*pow(r1/r , l+1.0)
       +pow(1.0/r, p);
   philm_real*= factor;

 return philm_real;
}

static double AnalyticPotential(const double r, const double theta, const double phi, const double r1, const double r2)
{
 double philm_real, philm_imag, pot, plm, cosmphi, sinmphi;

 philm_real = AnalyticPhilm(r, r1, r2);
 plm        = gsl_sf_legendre_sphPlm(LTEST, MTEST, cos(theta));

 if(MTEST == 0){
  return philm_real*plm;
 } else {
  cosmphi    = cos(phi*(double)MTEST);
  sinmphi    = sin(phi*(double)MTEST);

  philm_imag = -philm_real; /*true for the chosen density field*/

  pot = philm_real*cosmphi - philm_imag*sinmphi;
  pot*= 2.0*plm;
 
  return pot;
 }
}


void TestPoisson(Data *d, Grid *G)
/*!
  ************************************************************************ */
{
  int i,j,k;
  struct GRID *GY, *GX, *GZ;
  GX = G;
  GY = G + 1;
  GZ = G + 2;

  /*over-write density field*/
  
  KTOT_LOOP(k){
    JTOT_LOOP(j){
      ITOT_LOOP(i){
	d->Vc[RHO][k][j][i] = AnalyticDensity(GX->x[i], GY->x[j], GZ->x[k]);
      }
    }
  }

 /* GetRhom (d, G); */
 /* GetPhilm (d, G, LTEST, MTEST); */


   /*Solve Poisson */
  SolvePoisson(d, G);
    
  /* print the potential */

  int jtest = (JEND + JBEG)/2, ktest = (KEND + KBEG)/2;
  double r1, r2, potential, error;
  double dpot;

  r1 = GX->xl_glob[GX->gbeg];
  r2 = GX->xr_glob[GX->gend];

  ITOT_LOOP(i){
   potential = AnalyticPotential(GX->x[i], GY->x[jtest], GZ->x[ktest], r1, r2);
   error     = fabs((potential - d->Phi[ktest][jtest][i])/potential);
 
//   dpot      = AnalyticPotential(GX->x[i+1], GY->x[JBEG], GZ->x[KBEG], r1, r2)-AnalyticPotential(GX->x[i-1], GY->x[JBEG], GZ->x[KBEG], r1, r2);
//   error     = fabs((dpot - (d->Phi[KBEG][JBEG][i+1]-d->Phi[KBEG][JBEG][i-1]))/dpot);

     printn (1, "%d, %.16lf , %.16lf, %.16lf, %.16lf \n", i,
   	    GX->x[i], d->Phi[ktest][jtest][i], potential, error);
   }

  if(prank == 1) exit(1);
}
