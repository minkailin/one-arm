#include "pluto.h"

#if CHAR_LIMITING == NO
static void PPM_COEFF (Grid *, real **, real **, real **, real **, real **);

/* ********************************************************************** */
void PPM_COEFF (Grid *grid, real **a, real **b, real **c, 
                real **d, real **e)
/* 
 *
 * PURPOSE
 *
 *   Compute PPM interpolation coefficients for arbitrary 
 *   mesh size. Coefficients are stored in the 2D arrays 
 *   a, b, c, d, e. The first index refers to the dimension,
 *   the second to the grid. 
 *   Coefficients for the steepening are computed in the 
 *   STEEPEN function.
 *
 *
 *
 ************************************************************************ */
{
  int    dim, i, beg, end;
  double scrh, *dx;

  for (dim = 0; dim < DIMENSIONS; dim++){

  /* ----------------------------------------------
      Initialize the coefficients everywhere to 
      the default values for uniform grids.
     ---------------------------------------------- */

    for (i = 0; i < NMAX_POINT; i++){ 
      a[dim][i] = 0.5;
      b[dim][i] = c[dim][i] = 1.0/6.0;
      d[dim][i] = e[dim][i] = 0.5;
    }
    
  /* -----------------------------------------------
      for non-uniform grids, the coefficients will
      depends on the grid spacing and are computed
      on a slightly smaller stencil. 
      Since initialization is done only once at the
      beginning, this should be avoided in Chombo
     ----------------------------------------------- */

    #ifndef CH_SPACEDIM
     beg = grid[dim].lbeg - grid[dim].nghost + 1;
     end = grid[dim].lend + grid[dim].nghost - 2;
     dx  = grid[dim].dx;
     for (i = beg; i <= end; i++) {
 
       scrh = 1.0/(dx[i - 1] + dx[i] + dx[i + 1] + dx[i + 2]);
 
       b[dim][i] = dx[i + 1]*scrh*(dx[i + 1] + dx[i + 2])/(dx[i] + 2.0*dx[i + 1]);
       c[dim][i] = dx[i]*scrh*(dx[i - 1] + dx[i])/(2.0*dx[i] + dx[i + 1]);

       a[dim][i] = dx[i]/(dx[i] + dx[i + 1]) + 
              + 2.0*(dx[i + 1]*c[dim][i] - dx[i]*b[dim][i])/(dx[i] + dx[i + 1]);

       scrh = dx[i]/(dx[i - 1] + dx[i] + dx[i + 1]);
       d[dim][i] = scrh*(2.0*dx[i - 1] + dx[i])/(dx[i + 1] + dx[i]);
       e[dim][i] = scrh*(2.0*dx[i + 1] + dx[i])/(dx[i - 1] + dx[i]);
     }
    #endif
  }
}

/* ********************************************************************** */
void States (const State_1D *state, int beg, int end, Grid *grid)
/*
 *  PURPOSE
 *    
 *    provide piece-wise parabolic reconstruction inside each 
 *    cell. Notice that wl and wr are left and right states
 *    with respect to the INTERFACE, while am and ap (later
 *    in this function) refer to the cell center, that is:
 *
 *                    vL-> <-vR
 *      |--------*--------|--------*--------|
 *       <-vm   (i)   vp->       (i+1)
 *
 ************************************************************************ */
{
  int   i, nv;
  double  **vm, **vp, **vc, *a, *b, *c, *d, *e;
  double *dvp, *dvm, dp, dm, d2, dc;
  static double  **dvF, **dvlim, **vR;
  static double **aa, **bb, **cc, **dd, **ee;

/* --------------------------------------------------- 
     PPM stencil is +- 2 zones:

     |---*---|---*---|---*---|---*---|---*---|
        i-2     i-1      i      i+1     i+2

     Thus, only 3 ghost zones are necessary.
     However, if FLATTENING or STEEPENING are added, 
     one more zone is required and the number of 
     ghost zones becomes 4.

     We define the leftmost and rightmost indexes  
     in the domain as 

      beg - s = 0
      end + s = grid[g_dir].np_tot - 1

     where s is the required stencil.
   ---------------------------------------------------- */

  if (dvF == NULL){
    dvF   = ARRAY_2D(NMAX_POINT, NVAR, double);
    vR    = ARRAY_2D(NMAX_POINT, NVAR, double);
    dvlim = ARRAY_2D(NMAX_POINT, NVAR, double);
    aa    = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    bb    = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    cc    = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    dd    = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    ee    = ARRAY_2D(DIMENSIONS, NMAX_POINT, double);
    PPM_COEFF (grid, aa, bb, cc, dd, ee);
  }

  vm = state->vm;
  vp = state->vp;
  vc = state->v;

  a = aa[g_dir];
  b = bb[g_dir];
  c = cc[g_dir];
  d = dd[g_dir];
  e = ee[g_dir];

/* -------------------------------------------------------
   1.         Compute undivided differences
   ------------------------------------------------------- */

  for (i = beg-2; i <= end+1; i++) {
  for (nv = 0; nv < NVAR; nv++) {
    dvF[i][nv] = vc[i + 1][nv] - vc[i][nv];
  }}

/* ---------------------------------------------------------
   2. compute limited slopes with monotonized central limiter
   --------------------------------------------------------- */

  for (i = beg-1; i <= end+1; i++){
    for (nv = 0; nv < NVAR; nv++) {
      dp = dvF[i][nv]; dm = dvF[i-1][nv];
      if (dp*dm > 0.0){
        d2 = 2.0*ABS_MIN(dp,dm);
        dc = d[i]*dp + e[i]*dm;
        dvlim[i][nv] = ABS_MIN(d2,dc);
      }else{
        dvlim[i][nv] = 0.0;
      }
    }
  }

/* ---------------------------------------------------------
   3. define right interface single value
   --------------------------------------------------------- */

  for (i = beg - 1; i <= end; i++) {
  for (nv = 0; nv < NVAR; nv++){
    vR[i][nv] = vc[i][nv] + a[i]*dvF[i][nv]  
                + b[i]*dvlim[i][nv] - c[i]*dvlim[i + 1][nv];
  }}

/* ---------------------------------------------------------
   4. compute right and left interpolated values using PPM
   --------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    dvp = dvF[i]; dvm = dvF[i-1];
    #if SHOCK_FLATTENING == MULTID    
     if (CheckZone (i, FLAG_MINMOD)) {
       for (nv = NVAR; nv--;    ) {
         dc = MINMOD(dvp[nv], dvm[nv]);
         vp[i][nv] = vc[i][nv] + 0.5*dc;
         vm[i][nv] = vc[i][nv] - 0.5*dc;
       }
     }else
    #endif
    for (nv = 0; nv < NVAR; nv++){
      dp = vR[i][nv]   - vc[i][nv];
      dm = vR[i-1][nv] - vc[i][nv];

      if (dp*dm >= 0.0) dp = dm = 0.0;
      else{
        if      (fabs(dp) >= 2.0*fabs(dm)) dp = - 2.0*dm;
        else if (fabs(dm) >= 2.0*fabs(dp)) dm =  - 2.0*dp;
      }
      vp[i][nv] = vc[i][nv] + dp; 
      vm[i][nv] = vc[i][nv] + dm;
    }

    #if PHYSICS == RHD || PHYSICS == RMHD
     VelocityLimiter (vc[i], vp[i], vm[i]);
    #endif
  }

/* --------------------------------------------------------
          1D shock flattening
   -------------------------------------------------------- */

  #if SHOCK_FLATTENING == YES
   Flatten (state, beg, end, grid);
  #endif

/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg-1; i <= end; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

  #if TIME_STEPPING == CHARACTERISTIC_TRACING
   CharTracingStep(state, beg, end, grid);
  #endif

/* -------------------------------------------
    compute states in conservative variables
   ------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);
}
#endif /* CHAR_LIMITING == NO */



#if CHAR_LIMITING == YES

/* ---------------------------------------------------------
    Set PARABOLIC_LIM to 0,1,2 to apply the parabolic
    limiter of PPM to either characteristic (0),
    primitive (1) or both variables (2). 
   --------------------------------------------------------- */

#define PARABOLIC_LIM  1
/* *************************************************************************** */
void States (const State_1D *state, int beg, int end, Grid *grid)
/*
 *
 * PURPOSE:
 * 
 *   Compute 1D left and right interface states using piecewise
 *   polynomial reconstruction (PPM or WENO) in the characteristic
 *   variables.
 *   This requires the eigenvector decomposition of the quasi-linear
 *   form of the equations.
 *   Interpolation is carried out at time level t^n by extrapolating
 *   the cell center value to the interfaces using appropriate limiting
 *   constraints
 *
 * LAST MODIFIED
 *
 *   June 20, 2012 by A. Mignone (mignone@ph.unito.it)
 *
 *************************************************************************** */
{
  int    i, j, k, nv, S=1;
  double dtdx, dx, dx2;
  double dp, dm, dvp[NVAR], dvm[NVAR];
  double dwp[NVAR], dwm[NVAR];
  double dwm2[NVAR], dwm1[NVAR], dwp1[NVAR], dwp2[NVAR];
  double Smm, Smp, Spp;
  double *vp, *vm, *vc, **v;
  double **L, **R, *lambda;
  double tau, a0, a1, w0, w1;
  const double one_sixth = 1.0/6.0;
  static double  **dvF;

/* --------------------------------------------
       local array memory allocation
   -------------------------------------------- */

  if (dvF == NULL){
    dvF = ARRAY_2D(NMAX_POINT, NVAR, double);
  } 
  v = state->v;

/* ---------------------------------------------
    define some useful quantities, compute
    source term and undivided differences
   --------------------------------------------- */

  #if INTERPOLATION == PARABOLIC
   S = 2;
  #endif
  SoundSpeed2 (v, state->a2, state->h, beg, end, CELL_CENTER, grid);

  for (i = beg-S; i <= end+S-1; i++){    
    for (nv = NVAR; nv--;   ) dvF[i][nv] = v[i+1][nv] - v[i][nv];
  }

/* --------------------------------------------------------------
                    main spatial loop
   -------------------------------------------------------------- */

  for (i = beg; i <= end; i++){    

    dx   = grid[g_dir].dx[beg];
    dx2  = dx*dx;
    dtdx = g_dt/dx;

    vc     = state->v[i]; 
    vp     = state->vp[i];
    vm     = state->vm[i];
    L      = state->Lp[i];
    R      = state->Rp[i];
    lambda = state->lambda[i];

    PrimEigenvectors(vc, state->a2[i], state->h[i], lambda, L, R);
    #if NVAR != NFLX
     for (k = NFLX; k < NVAR; k++) lambda[k] = vc[VXn]; 
    #endif

    #if SHOCK_FLATTENING == MULTID    
     if (CheckZone (i, FLAG_MINMOD)) {
       for (nv = 0; nv < NVAR; nv++){  
         dvp[nv] = 0.5*(MINMOD(dvF[i][nv], dvF[i-1][nv]));
         vp[nv] = vc[nv] + dvp[nv];
         vm[nv] = vc[nv] - dvp[nv];
       }
       continue;
     }
    #endif  /* SHOCK_FLATTENING == MULTID */

  /* ------------------------------------------------------------------
     1. Project undivided difference of primitive variables on
        along characteristics.
        Compute limited characteristic increments dwp and dwm
        by suitable reconstruction and such that 

          wp = w + dwp
          wm = w + dwm
     ------------------------------------------------------------------ */

    #if INTERPOLATION == WENO3

   /* -- compute undivided differences and 
         reconstruct characteristic fields -- */

     PrimToChar(L, dvF[i-1], dwm1);
     PrimToChar(L, dvF[i  ], dwp1);
     for (k = 0; k < NVAR; k++){
       tau = (dwp1[k] - dwm1[k]); 
       tau = tau*tau;

       a0 = 1.0 + tau/(dx2 + dwp1[k]*dwp1[k]);
       a1 = 1.0 + tau/(dx2 + dwm1[k]*dwm1[k]);

       dwp[k] =  (a0*dwp1[k] + 0.5*a1*dwm1[k])/(2.0*a0 + a1);
       dwm[k] = -(a1*dwm1[k] + 0.5*a0*dwp1[k])/(2.0*a1 + a0);
     }

    #endif     /* INTERPOLATION == WENO3 */

    #if INTERPOLATION == PARABOLIC

     PrimToChar(L, dvF[i-2], dwm2);
     PrimToChar(L, dvF[i-1], dwm1);
     PrimToChar(L, dvF[i  ], dwp1);
     PrimToChar(L, dvF[i+1], dwp2);

     for (k = 0; k < NVAR; k++){  
       Smm = MC(dwm2[k], dwm1[k]);
       Smp = MC(dwm1[k], dwp1[k]);
       Spp = MC(dwp1[k], dwp2[k]);
       dwp[k] =  0.5*dwp1[k] - (Spp - Smp)*one_sixth;
       dwm[k] = -0.5*dwm1[k] - (Smp - Smm)*one_sixth;

      /* -- parabolic limiter using characteristics -- */

       #if PARABOLIC_LIM == 0 || PARABOLIC_LIM == 2
        if (dwp[k]*dwm[k] >= 0.0) dwm[k] = dwp[k] = 0.0;
        else{
          if      (fabs(dwp[k]) >= 2.0*fabs(dwm[k])) dwp[k] = -2.0*dwm[k];
          else if (fabs(dwm[k]) >= 2.0*fabs(dwp[k])) dwm[k] = -2.0*dwp[k];
        }
       #endif
     }

    #endif   /* INTERPOLATION == PARABOLIC */

  /* -------------------------------------------------------------------
     2. Project limited differences in characteristic variable on right
        eigenvectors to obtain the corresponding primitive quantities,
 
        dv = \sum dw.R
     ------------------------------------------------------------------- */

    for (nv = 0; nv < NFLX; nv++) {
      dp = dm = 0.0;
      for (k = 0; k < NFLX; k++){
        dp += dwp[k]*R[nv][k];
        dm += dwm[k]*R[nv][k];
      }
      dvp[nv] = dp;
      dvm[nv] = dm;
    }
    #if NVAR != NFLX
     for (nv = NFLX; nv < NVAR; nv++){
       dvp[nv] = dwp[nv];
       dvm[nv] = dwm[nv];
     }
    #endif 

  /* --------------------------------------------------------------------
     3. Build L/R states in primitive variables and apply parabolic
        limiter to primitive variables if required.
     -------------------------------------------------------------------- */

    for (nv = 0; nv < NVAR; nv++) {
      #if PARABOLIC_LIM == 1 || PARABOLIC_LIM == 2
       if (dvp[nv]*dvm[nv] >= 0.0) dvp[nv] = dvm[nv] = 0.0;
       else {
         if      (fabs(dvp[nv]) >= 2.0*fabs(dvm[nv])) dvp[nv] = -2.0*dvm[nv];
         else if (fabs(dvm[nv]) >= 2.0*fabs(dvp[nv])) dvm[nv] = -2.0*dvp[nv];
       }
      #endif
      vp[nv] = vc[nv] + dvp[nv];
      vm[nv] = vc[nv] + dvm[nv];
    }

  /* ------------------------------------------------------------------
     4. Make sure that left and right states at time t^n are
        physically admissible. If not, use linear reconstruction on
        density and pressure.
     ------------------------------------------------------------------ */

    if (vp[RHO] < 0.0 || vm[RHO] < 0.0) {
      dvp[RHO] = 0.5*(MINMOD(dvF[i][RHO], dvF[i-1][RHO]));
      dvm[RHO] = - dvp[RHO]; 
      vp[RHO] = vc[RHO] + dvp[RHO];
      vm[RHO] = vc[RHO] + dvm[RHO];
    }
    #if EOS != ISOTHERMAL && EOS != BAROTROPIC
     if (vp[PRS] < 0.0 || vm[PRS] < 0.0) {
       dvp[PRS] = 0.5*(MINMOD(dvF[i][PRS], dvF[i-1][PRS]));
       dvm[PRS] = - dvp[PRS];       
       vp[PRS] = vc[PRS] + dvp[PRS];
       vm[PRS] = vc[PRS] + dvm[PRS];
     }
    #endif

    #if PHYSICS == RHD || PHYSICS == RMHD
     VelocityLimiter (vc, vp, vm);
    #endif

  }  /* -- end main loop on grid points -- */

/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg-1; i <= end; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

  #if TIME_STEPPING == CHARACTERISTIC_TRACING
   CharTracingStep(state, beg, end, grid);
  #endif

/* -------------------------------------------
    compute states in conservative variables
   ------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);
}
#undef PARABOLIC_LIM
#endif /* CHAR_LIMITING == YES */
