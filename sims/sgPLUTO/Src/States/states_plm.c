/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute states using piece-wise linear interpolation.

  Provide piece-wise linear reconstruction inside each cell using 
  primitive or characteristic variables.
  Note that vL and vR are left and right states with respect to the 
  \e interface, while vm and vp refer to the cell \e center:
  \verbatim
                      vL(i)-> <-vR(i)
       |----------*----------|----------*----------|
        <-vm(i)  (i)  vp(i)->         (i+1)
  \endverbatim    
 
  The default setting (LIMITER == DEFAULT) applies a different limiter
  to each variable: 
  - MC for density
  - VanLeer for velcoity and magnetic field
  - MinMod for pressure.
  
  Otherwise the same limiter can be imposed to all variables from 
  definitions.h.
  Slope limiters are function of left, right and centered slopes,
  limiter=f(dp, dc, dm).

  In non-cartesian geometries, function values are defined on the 
  geometrical centers thus giving an effective non-uniform grid spacing
  even if zones are equidistant.
  This is implemented, at the moment, only for cylindrical geometry.
  In this case, limiters use the \e divided difference and an extra 
  monotonicity constraint is added just before constructing the states.
  Assuming a monotonically increasing sequence of values,
  \f[
    \left\{\begin{array}{lcl}
     V_{i,+} &=& V_i + \Delta V_i \delta f_{i,R} \le V_{i+1} \\ \noalign{\medskip}
     V_{i,-} &=& V_i + \Delta V_i \delta f_{i,L} \ge V_{i-1}  
     \end{array}\right. 
     \quad\Longrightarrow\quad
     \Delta V_i < \min\left(\frac{V_{i+1} - V_i}{\delta f_{i,R}},\,
                            \frac{V_i - V_{i-1}}{-\delta f_{i,L}},\,\right)
  \f]
  In the general case, \f$\Delta V = {\rm MINMOD}\left[\Delta V,
  {\rm ABS\_MIN}\left(\frac{V_{i+1} - V_i}{\delta f_{i,R}},\,
   \frac{V_i - V_{i-1}}{-\delta f_{i,L}},\,\right)\right]\f$.

  A stencil of 3 zones is required for all limiters except for
  the FOURTH_ORDER_LIM which requires  5 zones. 

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 27, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if CHAR_LIMITING == NO
static void FourthOrderLinear(const State_1D *, int, int, Grid *);

/* ********************************************************************* */
void States (const State_1D *state, int beg, int end, Grid *grid)
/*! 
 * Compute states using piecewise linear interpolation.
 *
 * \param [in] state pointer to a State_1D structure
 * \param [in] beg   starting point where vp and vm must be computed
 * \param [in] end   final    point where vp and vm must be computed
 * \param [in] grid  pointer to array of Grid structures
 *
 * \return This function has no return value.
 *
 ************************************************************************ */
{
  int    nv, i;
  double **v, **vp, **vm;
  double dv[NVAR], dvc[NVAR];
  double dp, dm, d2, dc, hscale, dtdx;
  #if GEOMETRY != CARTESIAN
   double *xg, *xR;
   double *dfg, *df2g, *dfL, *dfR;
  #endif
  static double **dvF;

  #if LIMITER == FOURTH_ORDER_LIM
   FourthOrderLinear(state, beg, end, grid);
   return;
  #endif

/* -----------------------------------------------------------
         memory allocation and pointer shortcuts
   ----------------------------------------------------------- */

  if (dvF == NULL) {
    dvF = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  #if GEOMETRY != CARTESIAN
   dfg  = grid[g_dir].dfg;
   df2g = grid[g_dir].df2g;
   dfL  = grid[g_dir].dfL;
   dfR  = grid[g_dir].dfR;
  #endif

/* ----------------------------------------------------------
    scale factors for non-Cartesian geometries
   ---------------------------------------------------------- */
 
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
   hscale = 1.0;
  #elif GEOMETRY == POLAR
   if (g_dir == IDIR || g_dir == KDIR) hscale = 1.0;
   else hscale = grid[IDIR].x[*g_i];
  #elif GEOMETRY == SPHERICAL
   if      (g_dir == IDIR) hscale = 1.0;
   else if (g_dir == JDIR) hscale = grid[IDIR].x[*g_i];
   else if (g_dir == KDIR) {
     hscale = grid[IDIR].x[*g_i]*sin(grid[JDIR].x[*g_j]);
   }
  #endif

  v  = state->v;
  vp = state->vp;
  vm = state->vm;

/* -------------------------------------------
    compute undivided forward (F) differences
   ------------------------------------------- */

  for (i = beg-1; i <= end; i++){
    #if GEOMETRY == CARTESIAN
     for (nv = NVAR; nv--;   ) dvF[i][nv] = v[i+1][nv] - v[i][nv];
    #else
     for (nv = NVAR; nv--;   ) dvF[i][nv] = (v[i+1][nv] - v[i][nv])*dfg[i];
    #endif   
  }

/* -------------------------------------------
    2. Main spatial loop
   ------------------------------------------- */

  for (i = beg; i <= end; i++){

  /* ----------------------------------------
      2a. Compute centered slopes 
     ---------------------------------------- */
     
    #if GEOMETRY == CARTESIAN
     for (nv = NVAR; nv--;  ) dvc[nv] = 0.5*(dvF[i][nv] + dvF[i-1][nv]);
    #else
     for (nv = NVAR; nv--;  ) dvc[nv] = (v[i+1][nv] - v[i-1][nv])*df2g[i];
    #endif

  /* ----------------------------------------
      2b. if shock-flattening is enabled, 
          use minmod limiter
     ---------------------------------------- */
     
    #if SHOCK_FLATTENING == MULTID
     if (CheckZone(i,FLAG_MINMOD)) {
       for (nv = NVAR; nv--; ){
         dp = dvF[i][nv]; dm = dvF[i-1][nv];
         dv[nv] = MINMOD(dp,dm);
         #if GEOMETRY == CARTESIAN
          vp[i][nv] = v[i][nv] + 0.5*dv[nv];
          vm[i][nv] = v[i][nv] - 0.5*dv[nv];
         #else
          vp[i][nv] = v[i][nv] + dv[nv]*dfR[i];
          vm[i][nv] = v[i][nv] + dv[nv]*dfL[i];
         #endif
       }
       #if PHYSICS == RHD || PHYSICS == RMHD
        VelocityLimiter (v[i], vp[i], vm[i]);
       #endif
       continue;
     }else if (CheckZone(i,FLAG_FLAT)){
       for (nv = NVAR; nv--; ){
         vp[i][nv] = vm[i][nv] = v[i][nv];
       }
       continue;
     }
    #endif    

  /* ----------------------------------------
      2c. Compute limited slopes
     ---------------------------------------- */

    #if LIMITER == DEFAULT

   /* -- density: monotonized central difference -- */

     nv = DN;
     dp = dvF[i][nv]; dm = dvF[i-1][nv];
     if (dp*dm > 0.0){
       dc  = dvc[nv];
       d2  = 2.0*ABS_MIN(dp, dm);
       dv[nv] = ABS_MIN(d2, dc);
     }else{
       dv[nv] = 0.0;
     }
    
   /* -- velocity: vanleer (harmonic) limiter -- */

     for (nv = VX; nv < VX+COMPONENTS; nv++){
       dp = dvF[i][nv]; dm = dvF[i-1][nv];
       d2 = dp*dm;
       dv[nv] = (d2 > 0.0 ? 2.0*d2/(dp + dm):0.0);
     }

     /* -- magnetic field: vanleer (harmonic) limiter -- */

     #if PHYSICS == MHD || PHYSICS == RMHD
      for (nv = BX; nv < BX+COMPONENTS; nv++){
/*
    #ifdef STAGGERED_MHD
      if (nv == BXn) continue;
    #endif
*/
        dp = dvF[i][nv]; dm = dvF[i-1][nv];
        d2 = dp*dm;
        dv[nv] = (d2 > 0.0 ? 2.0*d2/(dp + dm):0.0);
      }
      #ifdef GLM_MHD /* -- mc limiter -- */
       nv = PSI_GLM;
       dp = dvF[i][nv]; dm = dvF[i-1][nv];
       if (dp*dm > 0.0){
         dc  = dvc[nv];      
         d2  = 2.0*ABS_MIN(dp, dm);
         dv[nv] = ABS_MIN(d2, dc);
       }else{
         dv[nv] = 0.0;
       }
      #endif
     #endif

   /* -- pressure: minmod limiter -- */
      
     #if EOS != ISOTHERMAL && EOS != BAROTROPIC
      nv = PR;
      dp = dvF[i][nv]; dm = dvF[i-1][nv];
      dv[nv] = MINMOD(dp, dm);
     #endif

     #if NFLX != NVAR   /* -- scalars: MC limiter -- */
      for (nv = NFLX; nv < NVAR; nv++){
        dp = dvF[i][nv]; dm = dvF[i-1][nv];
        if (dp*dm > 0.0){
          dc = dvc[nv];
          d2 = 2.0*ABS_MIN(dp, dm);
          dv[nv] = ABS_MIN(d2, dc);
        }else{
          dv[nv] = 0.0;
        }
      }
     #endif
      
    #else  /* -- use same limiter for all variables -- */

     for (nv = 0; nv < NVAR; nv++){
       dp = dvF[i][nv]; dm = dvF[i-1][nv];
       if (dp*dm > 0.0){
         #if LIMITER == FLAT_LIM
          dv[nv] = 0.0;
         #elif LIMITER == MINMOD_LIM
          dv[nv] = ABS_MIN(dp, dm);
         #elif LIMITER == VANALBADA_LIM
          #define EPS_VA 1.e-18
          double dp2  = dp*dp;
          double dm2  = dm*dm;
          dv[nv] = (dp*(dm2 + EPS_VA) + dm*(dp2 + EPS_VA))
                   /(dp2 + dm2 + EPS_VA);
          #undef EPS_VA
         #elif LIMITER == UMIST_LIM
          double ddp, ddm;
          ddp = 0.25*(dp + 3.0*dm);
          ddm = 0.25*(dm + 3.0*dp);
          d2  = 2.0*(fabs(dp) < fabs(dm) ? dp:dm);
          d2  = (fabs(d2) < fabs(ddp) ? d2:ddp);
          dv[nv] = (fabs(d2) < fabs(ddm) ? d2:ddm);
         #elif LIMITER == VANLEER_LIM
          dv[nv] = 2.0*dp*dm/(dp + dm);
         #elif LIMITER == MC_LIM
          dc = dvc[nv];
          d2     = 2.0*ABS_MIN(dp, dm);
          dv[nv] = ABS_MIN(d2, dc);
         #else
          print1 ("! Reconstruct: limiter not defined\n");
          QUIT_PLUTO(1);
         #endif
       }else{
         dv[nv] = 0.0;
       }
     }
    #endif

  /* ----------------------------------------
      2d. Compute + and - states 
     ---------------------------------------- */

    for (nv = NVAR; nv--;   ){
      #if GEOMETRY == CARTESIAN
       vp[i][nv] = v[i][nv] + 0.5*dv[nv];
       vm[i][nv] = v[i][nv] - 0.5*dv[nv];
      #else

    /* --------------------------------------------------------------
        In non-cartesian geometry we apply apply extra-monotonicity 
        constraints since interpolation is done on centroids:
       -------------------------------------------------------------- */
       
       #if GEOMETRY == CYLINDRICAL
        dp = v[i+1][nv] - v[i][nv];
        dm = v[i][nv] - v[i-1][nv];
        d2 = ABS_MIN(dp/dfR[i], dm/fabs(dfL[i]));
        dv[nv] = MINMOD(d2, dv[nv]);
       #endif

       vp[i][nv] = v[i][nv] + dv[nv]*dfR[i];
       vm[i][nv] = v[i][nv] + dv[nv]*dfL[i];
      #endif
    }

  /* --------------------------------------
      2e.  Relativistic Limiter
     -------------------------------------- */

    #if PHYSICS == RHD || PHYSICS == RMHD
     VelocityLimiter (v[i], vp[i], vm[i]);
    #endif
  } /* -- end loop on zones -- */

/* -------------------------------------------
    3. Shock flattening 
   -------------------------------------------  */

  #if SHOCK_FLATTENING == ONED
   Flatten (state, beg, end, grid);
  #endif

/* -------------------------------------------
    4.  Assign face-centered magnetic field
   -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg - 1; i <= end; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

/* --------------------------------------------------------
    5. evolve cell-center values by dt/2
   -------------------------------------------------------- */

  #if TIME_STEPPING == CHARACTERISTIC_TRACING
   CharTracingStep(state, beg, end, grid);
  #elif TIME_STEPPING == HANCOCK
   HancockStep(state, beg, end, grid);
  #endif

/* ---------------------------------------------
    6. compute states in conservative variables
   --------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);
}
/* ********************************************************************** */
void FourthOrderLinear(const State_1D *state, int beg, int end, Grid *grid)
/*
 * PURPOSE
 *
 *   Compute interface states using Colella's fourth-order slope limiter
 * 
 *  Ref:  Miller, G.H and P. COlella, 
 *        "A high-order Eulerian Godunov Method for 
 *         Elastic-Plastic Flow in Solids", JCP 167,131-176 (2001)
 *    
 *                             +
 *
 *        Saltzman, J, " An Unsplit 3D Upwind Method for 
 *                       Hyperbolic Conservation Laws", 
 *                       JCP 115, 153-168 (1994)
 *
 *********************************************************************** */
{
  int    i, nv;
  static double **s;
  static double **dv, **dvf, **dvc, **dvlim; 
  double scrh, dvp, dvm, dvl;
  double **v, **vp, **vm;

  v  = state->v;
  vp = state->vp;
  vm = state->vm;

  if (s == NULL){
    s     = ARRAY_2D(NMAX_POINT, NVAR, double);
    dv    = ARRAY_2D(NMAX_POINT, NVAR, double);
    dvf   = ARRAY_2D(NMAX_POINT, NVAR, double);
    dvc   = ARRAY_2D(NMAX_POINT, NVAR, double);
    dvlim = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  #if TIME_STEPPING == HANCOCK && PHYSICS != RMHD
   SoundSpeed2 (state->v, state->a2, state->h, beg, end, CELL_CENTER, grid);
  #endif
/*
  #if GEOMETRY != CARTESIAN
   print1 ("! FourthOrderLinear: only Cartesian geometry supported\n");
   QUIT_PLUTO(1);  
  #endif
*/
/* -----------------------------------------------------------
               compute undivided differences
   ----------------------------------------------------------- */

  for (i = beg-2; i <= end+1; i++){
    for (nv = 0; nv < NVAR; nv++) dv[i][nv] = v[i+1][nv] - v[i][nv];
  }

  for (i = beg - 1; i <= end + 1; i++){
  for (nv = 0; nv < NVAR; nv++){
    dvp = dv[i][nv]; dvm = dv[i-1][nv];
    dvc[i][nv] = 0.5*(dvp + dvm);
      s[i][nv] =  (dvp > 0.0 ? 0.5:-0.5) 
                + (dvm > 0.0 ? 0.5:-0.5);
    dvlim[i][nv] = 2.0*MIN(fabs(dvp), fabs(dvm));
    dvf[i][nv]   = MIN(fabs(dvc[i][nv]), dvlim[i][nv])*s[i][nv];
  }}

  for (i = beg; i <= end; i++){
    for (nv = 0; nv < NVAR; nv++){
      if (dv[i][nv]*dv[i-1][nv] > 0.0) {
        scrh = 4.0/3.0*dvc[i][nv] - (dvf[i+1][nv] + dvf[i-1][nv])/6.0; 
        dvlim[i][nv]  = MIN(fabs(scrh), dvlim[i][nv])*s[i][nv];
      }else{
        dvlim[i][nv] = 0.0;
      }
    }

    for (nv = NVAR; nv--;  ) {
      vp[i][nv] = v[i][nv] + 0.5*dvlim[i][nv];
      vm[i][nv] = v[i][nv] - 0.5*dvlim[i][nv];
    }
    #if PHYSICS == RHD || PHYSICS == RMHD
     VelocityLimiter(v[i], vp[i], vm[i]);
    #endif
  }

/*  -------------------------------------------
               Shock flattening
    -------------------------------------------  */

  #if SHOCK_FLATTENING == MULTID && CHAR_LIMITING == NO
   Flatten (state, beg, end, grid);
  #endif

/*  -------------------------------------------
        Shock flattening
    -------------------------------------------  */

  #if SHOCK_FLATTENING == ONED
   Flatten (state, beg, end, grid);
  #endif

/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg - 1; i <= end; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

/* --------------------------------------------------------
      evolve center values by dt/2
   -------------------------------------------------------- */


  #if TIME_STEPPING == HANCOCK
   HancockStep(state, beg, end, grid);
  #endif

/* -------------------------------------------
    compute states in conservative variables
   ------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);
}
#endif  /* CHAR_LIMITING == NO */



#if CHAR_LIMITING == YES
/* *********************************************************************** */
void States (const State_1D *state, int beg, int end,  Grid *grid)
/*
 *
 * PURPOSE:
 * 
 *   Compute 1D left and right interface states using piecewise
 *   linear reconstruction and the characteristic decomposition of the
 *   quasi-linear form of the equations.
 *
 *   This is done by first extrapolating the cell center value to the 
 *   interface using piecewise limited linear reconstruction
 *   on the characteristic variables.
 *
 *   Left and right states are then evolved for the half time step 
 *   using characteristic tracing if necessary.
 *
 * LAST MODIFIED
 *
 *   May 22, 2012 by A. Mignone (mignone@ph.unito.it)
 *
 ************************************************************************* */
{
  int    i, j, k, nv;
  double dx, dtdx, hscale;
  double d2w, dwc[NVAR], dw[NVAR], dwp[NVAR], dwm[NVAR];
  double *dvp, *dvm, dv[NVAR];
  double dp, dm, dc;
  double *vp, *vm, *vc, **v;
  double **L, **R, *lambda;
  double kstp[NVAR];
  #if GEOMETRY != CARTESIAN
   double *dfg, *df2g, *dfL, *dfR;
   double betaL[NVAR], betaR[NVAR];
  #endif
  static double **dvF;

/* --------------------------------------------
    allocate memory and set pointer shortcuts
   -------------------------------------------- */

  if (dvF == NULL){
    dvF = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  v = state->v; 
  #if GEOMETRY != CARTESIAN
   dfg  = grid[g_dir].dfg;
   df2g = grid[g_dir].df2g;
   dfL  = grid[g_dir].dfL;
   dfR  = grid[g_dir].dfR;
  #endif

/* ---------------------------------------------
    define some useful quantities, compute
    source term and undivided differences
   --------------------------------------------- */
 
  SoundSpeed2 (v, state->a2, state->h, beg, end, CELL_CENTER, grid);

  for (i = beg-1; i <= end; i++){
    #if GEOMETRY == CARTESIAN
     for (nv = NVAR; nv--;   ) dvF[i][nv] = v[i+1][nv] - v[i][nv];
    #else
     for (nv = NVAR; nv--;   ) dvF[i][nv] = (v[i+1][nv] - v[i][nv])*dfg[i];
    #endif
  }

/* --------------------------------------------------------------
    set the amount of steepening for each characteristic family.
    Default is 2, but nonlinear fields may be safely set to 1
    for strongly nonlinear problems.
   -------------------------------------------------------------- */

  for (k = NVAR; k--;  ) kstp[k] = 2.0;
  #if PHYSICS == HD || PHYSICS == RHD
   kstp[0] = kstp[1] = 1.0;
  #elif PHYSICS == MHD
   kstp[KFASTP] = kstp[KFASTM] = 1.0;
   #if COMPONENTS > 1
    kstp[KSLOWP] = kstp[KSLOWM] = 1.0; 
   #endif
  #endif

/* ----------------------------------------------------------
    scale factors for non-Cartesian geometries
   ---------------------------------------------------------- */
 
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
   hscale = 1.0;
  #elif GEOMETRY == POLAR
   if (g_dir == IDIR || g_dir == KDIR) hscale = 1.0;
   else hscale = grid[IDIR].x[*g_i];
  #elif GEOMETRY == SPHERICAL
   if      (g_dir == IDIR) hscale = 1.0;
   else if (g_dir == JDIR) hscale = grid[IDIR].x[*g_i];
   else if (g_dir == KDIR) {
     hscale = grid[IDIR].x[*g_i]*sin(grid[JDIR].x[*g_j]);
   }
  #endif

/* --------------------------------------------------------------
                    main spatial loop
   -------------------------------------------------------------- */

  for (i = beg; i <= end; i++){    

    vp     = state->vp[i];
    vm     = state->vm[i];
    vc     = state->v[i];
    L      = state->Lp[i];
    R      = state->Rp[i];
    lambda = state->lambda[i];

    dx   = hscale*grid[g_dir].dx[i];
    dtdx = g_dt/dx;

    PrimEigenvectors(vc, state->a2[i], state->h[i], lambda, L, R);

  /* ---------------------------------------------------------------
     1. Project forward, backward (and centered) undivided 
        differences of primitive variables along characteristics, 
        dw(k) = L(k).dv
     --------------------------------------------------------------- */

    dvp = dvF[i]; dvm = dvF[i-1];  /* -- pointer shorcuts -- */
    PrimToChar(L, dvm, dwm);
    PrimToChar(L, dvp, dwp);
    #if GEOMETRY == CARTESIAN
     for (k = NVAR; k--;    )  dwc[k] = 0.5*(dwm[k] + dwp[k]);
    #else
     for (nv = NVAR; nv--;   ) dv[nv] = (v[i+1][nv] - v[i-1][nv])*df2g[i];
     PrimToChar(L, dv, dwc);
    #endif

  /* ----------------------------------------------------------
     2. Apply slope limiter to characteristic differences for
        nv < NFLX, dw = Limiter(dwp, dwm). 
     ------------------------------------------------------- */

    #if SHOCK_FLATTENING == MULTID
     if (CheckZone (i, FLAG_MINMOD)) {
       for (k = NFLX; k--;    ) {
         dp = dwp[k]; dm = dwm[k];
         dw[k] = MINMOD(dp, dm);
       }
     } else
    #endif
    for (k = NFLX; k--;    ){
      #if LIMITER == DEFAULT

       if (dwp[k]*dwm[k] > 0.0) {
         dwc[k] = 0.5*(dwm[k] + dwp[k]);
         d2w    = kstp[k]*ABS_MIN(dwp[k],dwm[k]);
         dw[k]  = ABS_MIN(d2w, dwc[k]);
       }else dw[k] = 0.0;

      #else 

       dp = dwp[k]; dm = dwm[k];
       if (dp*dm > 0.0){
         #if LIMITER == FLAT_LIM
          dw[k] = 0.0;
         #elif LIMITER == MINMOD_LIM
          dw[k] = (fabs(dp) < fabs(dm) ? dp:dm);
         #elif LIMITER == VANALBADA_LIM
          #define EPS_VA 1.e-18
          double dp2  = dp*dp;
          double dm2  = dm*dm;
          dw[k] = (dp*(dm2 + EPS_VA) + dm*(dp2 + EPS_VA))
                    /(dp2 + dm2 + EPS_VA);
          #undef EPS_VA
         #elif LIMITER == VANLEER_LIM
          dw[k] = 2.0*dp*dm/(dp + dm);
         #elif LIMITER == MC_LIM
          dc = dwc[k];
          d2w   = 2.0*ABS_MIN(dp, dm);
          dw[k] = ABS_MIN(d2w, dc);
         #else
          print1 ("! CharSlopes: limiter not defined or available\n");
          QUIT_PLUTO(1);
         #endif
       }else{
         dw[k] = 0.0;
       }
      #endif
    }

  /* ------------------------------------------------------------------
     3. Project limited slopes in characteristic variables on right 
        eigenvectors to obtain primitive slopes: dv = \sum dw.R

        Also, enforce monotonicity in primitive variables as well.
     ------------------------------------------------------------------ */

    for (nv = NFLX; nv--;   ){
      dc = 0.0;
      for (k = 0; k < NFLX; k++) dc += dw[k]*R[nv][k];

      #if GEOMETRY != CYLINDRICAL 
       dp = dvp[nv]; dm = dvm[nv];
       if (dp*dm > 0.0){
         d2w    = 2.0*ABS_MIN(dp, dm);
         dv[nv] = MINMOD(d2w, dc);
       }else dv[nv] = 0.0;
      #else
       dp = v[i+1][nv] - v[i][nv]; 
       dm = v[i][nv] - v[i-1][nv];
       if (dp*dm > 0.0){
         d2w    = ABS_MIN(dp/dfR[i], -dm/dfL[i]);
         dv[nv] = MINMOD(d2w, dc);
       }else dv[nv] = 0.0;      
      #endif
    }

  /* -----------------------------------------------------------------
     4. Repeat construction for passive scalars (tracers).
        For a passive scalar, the primitive variable is the same as
        the characteristic one. We use the MC limiter always.
     ----------------------------------------------------------------- */
     
    #if NFLX != NVAR
     for (nv = NFLX; nv < NVAR; nv++ ){
       dp = dvp[nv]; dm = dvm[nv];
       if (dp*dm > 0.0){
         dc     = dwc[nv];
         d2w    = 1.3*ABS_MIN(dp, dm);
         dv[nv] = ABS_MIN(d2w, dc);
       }else{
         dv[nv] = 0.0;
       }
     }
    #endif

  /* --------------------------------------------------------------------
     5. Build L/R states at time level t^n
     -------------------------------------------------------------------- */

    for (nv = NVAR; nv--;   ) {
      #if GEOMETRY == CARTESIAN
       vp[nv] = vc[nv] + 0.5*dv[nv];
       vm[nv] = vc[nv] - 0.5*dv[nv];
      #else
       vp[nv] = vc[nv] + dv[nv]*dfR[i];
       vm[nv] = vc[nv] + dv[nv]*dfL[i];
      #endif
    }

    #if PHYSICS == RHD || PHYSICS == RMHD
     VelocityLimiter (vc, vp, vm);
    #endif

  }  /* -- end main loop on grid points -- */

/*  -------------------------------------------
        Shock flattening (only 1D)
    -------------------------------------------  */

  #if SHOCK_FLATTENING == ONED
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

/* --------------------------------------------------------
      evolve center values by dt/2
   -------------------------------------------------------- */

  #if TIME_STEPPING == CHARACTERISTIC_TRACING
   CharTracingStep(state, beg, end, grid);
  #elif TIME_STEPPING == HANCOCK
   HancockStep(state, beg, end, grid);
  #endif

/* -------------------------------------------
    compute states in conservative variables
   ------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);
}
#undef REF_STATE 
#endif  /* CHAR_LIMITING == YES */
