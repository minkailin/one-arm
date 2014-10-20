/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Initialize geometry-dependent grid quantities.

  Compute grid quantities (such as interface area, volume, geometrical
  center, etc..) that depend on the geometry.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 23, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void MakeGeometry (Grid *GXYZ)
/*!
 *
 * \param [in,out] GXYZ  Pointer to an array of Grid structures;
 *
 *********************************************************************** */
{
  int     i, j, k, idim, ngh, ileft;
  int     iright;
  int     iend, jend, kend;
  double  xiL, xiR, dxi, dvol;
  double  x, dx, xr, xl;
  double  y, dy, yr, yl;
  struct  GRID *GG;

  iend = GXYZ[0].lend + GXYZ[0].nghost;
  jend = GXYZ[1].lend + GXYZ[1].nghost;
  kend = GXYZ[2].lend + GXYZ[2].nghost;

/*  --------------------------------------------------------------
     Memory allocation. All values are defined at the cell center 
     with the exception of the area element which is defined on a
     staggered mesh and therefore starts at [-1].
    ----------------------------------------------------------- */

  for (idim = 0; idim < 3; idim++) {
    (GXYZ + idim)->A       = ARRAY_1D (GXYZ[idim].np_tot+1, double)+1;
    (GXYZ + idim)->xgc     = ARRAY_1D (GXYZ[idim].np_tot, double);
    (GXYZ + idim)->dV      = ARRAY_1D (GXYZ[idim].np_tot, double);
    (GXYZ + idim)->r_1     = ARRAY_1D (GXYZ[idim].np_tot, double);
    (GXYZ + idim)->ct      = ARRAY_1D (GXYZ[idim].np_tot, double);
    (GXYZ + idim)->inv_dx  = ARRAY_1D (GXYZ[idim].np_tot, double);
    (GXYZ + idim)->inv_dxi = ARRAY_1D (GXYZ[idim].np_tot, double);

  /* -- Interpolation coefficients and 
        default values on Cartesian meshes -- */

    (GXYZ + idim)->dfg     = ARRAY_1D (GXYZ[idim].np_tot, double);
    (GXYZ + idim)->df2g    = ARRAY_1D (GXYZ[idim].np_tot, double);
    (GXYZ + idim)->dfL     = ARRAY_1D (GXYZ[idim].np_tot, double);
    (GXYZ + idim)->dfR     = ARRAY_1D (GXYZ[idim].np_tot, double);
    for (i = 0; i < GXYZ[idim].np_tot; i++){
      GXYZ[idim].dfg[i]  = 1.0;    /* = dx/(xg(i+1) - xg(i))   */
      GXYZ[idim].df2g[i] =  0.5;   /* = dx/(xg(i+1) - xg(i-1)) */
      GXYZ[idim].dfL[i]  = -0.5;   /* = (xl - xg(i))/dx        */
      GXYZ[idim].dfR[i]  =  0.5;   /* = (xr - xg(i))/dx        */
    }
  }

/* ------------------------------------------------------------
    Define area (A), volume element (dV) and cell geometrical
    centers for each direction.
    
    Conventions:
    - dV is always positive 
    - xgc can be negative if x < 0 
    ----------------------------------------------------------- */  
 
  GG = GXYZ;
  for (i = 0; i <= iend; i++) {

    dx  = GG->dx[i];
    x   = GG->x[i];
    xr  = x + 0.5*dx;
    xl  = x - 0.5*dx;

    #if GEOMETRY == CARTESIAN
     GG->A[i]    = 1.0;
     if (i == 0) GG->A[-1] = 1.0;
     GG->dV[i]   = dx;
     GG->xgc[i]  = x;
 
    #elif GEOMETRY == CYLINDRICAL
     GG->A[i]    = fabs(xr); 
     if (i == 0) GG->A[-1] = fabs(xl);

     GG->dV[i]  = fabs(x)*dx;
     #ifdef NEW_GEOM
      GG->xgc[i] = x + dx*dx/(12.0*x); 
     #else
      GG->xgc[i]  = DSIGN(x)*sqrt(0.5*(xr*xr + xl*xl)); 
     #endif

     GG->r_1[i] = 1.0/x;
     
     #ifdef NEW_GEOM
      GG->dfL[i]  = (xl - GG->xgc[i])/dx; 
      GG->dfR[i]  = (xr - GG->xgc[i])/dx;
     #endif
     
    #elif GEOMETRY == POLAR
     GG->A[i]    = fabs(xr); 
     if (i == 0) GG->A[-1] = fabs(xl);

     GG->dV[i]  = fabs(x)*dx;
/*     GG->xgc[i] = x + dx*dx/(12.0*x);  */
GG->xgc[i]  = DSIGN(x)*sqrt(0.5*(xr*xr + xl*xl));

     GG->r_1[i] = 1.0/x;
/*
     GG->dfL[i]  = (xl - GG->xgc[i])/dx; 
     GG->dfR[i]  = (xr - GG->xgc[i])/dx;
*/
    #elif GEOMETRY == SPHERICAL
     GG->A[i]    = xr*xr;
     if (i == 0) GG->A[-1] = xl*xl;

     GG->dV[i]   = fabs(xr*xr*xr - xl*xl*xl)/3.0;
     GG->xgc[i]  = x + 2.0*x*dx*dx/(12.0*x*x + dx*dx);
GG->xgc[i]  = DSIGN(x)*pow(fabs(0.5*(xr*xr*xr + xl*xl*xl)), 1.0/3.0);

/*     GG->r_1[i]  = 1.5*(xr*xr - xl*xl)/(xr*xr*xr - xl*xl*xl); */
     GG->r_1[i] = 1.0/x;
/*
     GG->dfL[i]  = (xl - GG->xgc[i])/dx; 
     GG->dfR[i]  = (xr - GG->xgc[i])/dx;
*/
    #endif
  }

  GG = GXYZ + 1;
  for (j = 0; j <= jend; j++) {

    dx  = GG->dx[j];
    x   = GG->x[j];
    xr  = x + 0.5*dx;
    xl  = x - 0.5*dx;

    #if GEOMETRY != SPHERICAL
     GG->A[j]    = 1.0;
     if (j == 0) GG->A[-1] = 1.0;
     GG->dV[j]   = dx;
     GG->xgc[j]  = x;
    #else
     GG->A[j]    = fabs(sin(xr));
     if (j == 0) GG->A[-1] = fabs(sin(xl));

     GG->dV[j]  = fabs(cos(xl) - cos(xr));

     GG->xgc[j]  = DSIGN(x)*acos(0.5*(cos(xr) + cos(xl)));

/* -- using f(sin\theta) 
     GG->xgc[j] = 0.5*(sin(xl)*cos(xl) - sin(xp)*cos(xp) + dx)/GG->dV[j];*/

/* -- using f(\theta)  */
/*
     GG->xgc[j]  = (sin(xr) - sin(xl) + xl*cos(xl) - xr*cos(xr));
     GG->xgc[j] /= cos(xl) - cos(xr);
*/
     GG->ct[j]  = 1.0/tan(x);  /* (sin(xr) - sin(xl))/(cos(xl) - cos(xr)); */
/*
     GG->dfL[i]  = (xl - GG->xgc[i])/dx; 
     GG->dfR[i]  = (xr - GG->xgc[i])/dx;
*/
    #endif
  }
  
  GG = GXYZ + 2;
  for (k = 0; k <= kend; k++) {

    dx  = GG->dx[k];
    x   = GG->x[k];
    xr  = x + 0.5*dx;
    xl  = x - 0.5*dx;

    GG->A[k]   = 1.0;
    if (k == 0) GG->A[-1] = 1.0;
    GG->dV[k]  = dx;
    GG->xgc[k] = x;
  }

/* ---------------------------------------------------------
    compute and store the reciprocal of cell spacing
    between interface (inv_dx) and cell centers (inv_dxi)
   --------------------------------------------------------- */

  #if GEOMETRY == CYLINDRICAL && (defined NEW_GEOM)

   GG = GXYZ;
   for (i = 0; i < iend; i++) {
     dx = GG->dx[i];
     GG->dfg[i]  = dx/(GG->xgc[i+1] - GG->xgc[i]);
   }

   for (i = 1; i < iend; i++) {
     dx = GG->dx[i];
     GG->df2g[i] = dx/(GG->xgc[i+1] - GG->xgc[i-1]);
   }

  #endif

  #if GEOMETRY == SPHERICAL
/*
   GG = GXYZ+1;
   for (j = 0; j < jend; j++) {
     dx = GG->dx[j];
     GG->dfg[j]  = dx/(GG->xgc[j+1] - GG->xgc[j]);
   }

   for (j = 1; j < jend; j++) {
     dx = GG->dx[j];
     GG->df2g[j] = dx/(GG->xgc[j+1] - GG->xgc[j-1]);
   }
*/
  #endif

  for (idim = 0; idim < DIMENSIONS; idim++){
    for (i = 0; i < GXYZ[idim].np_tot; i++) {
      GXYZ[idim].inv_dx[i] = 1.0/(GXYZ[idim].dx[i]);
    }

    for (i = 0; i < GXYZ[idim].np_tot-1; i++) {
      GXYZ[idim].inv_dxi[i] = 2.0/(GXYZ[idim].dx[i] + GXYZ[idim].dx[i+1]);
    }
  }
  
}

/* ********************************************************************** */
real Length_1 (int i, int j, int k, Grid *grid)
/*
 *
 *
 *
 *
 *
 ************************************************************************ */
{
  return (grid[0].dx[i]);
}

/* ********************************************************************** */
real Length_2 (int i, int j, int k, Grid *grid)
/*
 *
 *
 *
 *
 *
 ************************************************************************ */
{
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
   return (grid[1].dx[j]);
  #elif GEOMETRY == POLAR ||  GEOMETRY == SPHERICAL
   return (fabs(grid[0].xgc[i])*grid[1].dx[j]);
  #endif
}

/* ********************************************************************** */
real Length_3 (int i, int j, int k, Grid *grid)
/*
 *
 *
 *
 *
 *
 ************************************************************************ */
{
  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
   return (grid[2].dx[k]);
  #elif GEOMETRY == CYLINDRICAL
   return (fabs(grid[0].xgc[i])*grid[2].dx[k]);
  #elif GEOMETRY == SPHERICAL
   return (fabs(grid[0].xgc[i]*sin(grid[1].xgc[j]))*grid[2].dx[k]);
  #endif
}

/* ******************************************************************* */
double *GetInverse_dl (const Grid *grid)
/*
 *
 *  Return an array of (inverse) physical cell lengths in the
 *  direction given by g_dir.
 *  For spherical coordinates, for instance this will be
 *
 *    {dr}_i                      if g_dir == IDIR
 *    {r_i*dtheta}_j              if g_dir == JDIR
 *    {r_i*sin(theta_j)*dphi}_k   if g_dir == KDIR
 *   
 *
 ********************************************************************* */
{
#if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL

  return grid[g_dir].inv_dx;

#elif GEOMETRY == POLAR

  if (g_dir != JDIR){
    return grid[g_dir].inv_dx;
  }else{
    int    j;
    double r_1;
    static double *inv_dl;
   
    if (inv_dl == NULL) inv_dl = ARRAY_1D(NX2_TOT, double);
    r_1 = grid[IDIR].r_1[*g_i];
    JTOT_LOOP(j) inv_dl[j] = grid[JDIR].inv_dx[j]*r_1;
    return inv_dl;
  }

#elif GEOMETRY == SPHERICAL

  int    j, k;
  double r_1, s;
  static double *inv_dl2, *inv_dl3;

  if (inv_dl2 == NULL) {
    inv_dl2 = ARRAY_1D(NX2_TOT, double);
    inv_dl3 = ARRAY_1D(NX3_TOT, double);
  }

  if (g_dir == IDIR){
    return grid[IDIR].inv_dx;
  }else if (g_dir == JDIR) {
    r_1 = grid[IDIR].r_1[*g_i];
    JTOT_LOOP(j) inv_dl2[j] = grid[JDIR].inv_dx[j]*r_1;
    return inv_dl2;
  }else if (g_dir == KDIR){
    r_1 = grid[IDIR].r_1[*g_i];
    s   = grid[JDIR].x[*g_j];
    s   = sin(s);
    KTOT_LOOP(k) inv_dl3[k] = grid[KDIR].inv_dx[k]*r_1/s;
    return inv_dl3;
  }

#endif

  return NULL;
}
