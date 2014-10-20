#include"pluto.h"

/* ********************************************************************* */
void HLL_Solver (const State_1D *state, int beg, int end, 
                 real *cmax, Grid *grid)
/*
 *
 * NAME
 *
 *   HLL_SOLVER
 *
 *
 * PURPOSE
 *
 *   Solve riemann problem for the Euler equations 
 *   using th HLLE solver;
 * 
 *   Reference:    --
 *   
 *
 * ARGUMENTS
 *
 *   vL (IN)       1-D array of left-edge primitive values at i+1/2
 *   vR (IN)       1-D array of right-edge primitive values at i+1/2
 *   flux (OUT)    1-D array of numerical fluxes, not including pressure
 *   grid (IN)     Array of grids
 *   cmax(OUT)     Array of maximum characteristic speeds in this direction
 *
 *
 *
 * LAST_MODIFIED
 *
 *   August, 17/2005, written by Andrea Mignone  (mignone@to.astro.it)
 *              
 *
 *
 **************************************************************************** */
{
  int    nv, i;
  real   scrh, bmin, bmax;
  double *uL, *uR, *SL, *SR;
  static double **fL, **fR;
  static double *pL, *pR, *a2L, *a2R, *hL, *hR;
  static double **Uhll;
  static double **VL, **VR, **UL, **UR;

  if (fL == NULL){
    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);
    pL = ARRAY_1D(NMAX_POINT, double);
    pR = ARRAY_1D(NMAX_POINT, double);
    SL = ARRAY_1D(NMAX_POINT, double);
    SR = ARRAY_1D(NMAX_POINT, double);

    a2L = ARRAY_1D(NMAX_POINT, double);
    a2R = ARRAY_1D(NMAX_POINT, double);

    hL = ARRAY_1D(NMAX_POINT, double);
    hR = ARRAY_1D(NMAX_POINT, double);

    Uhll = ARRAY_2D(NMAX_POINT, NVAR, double); /* NVAR and not NFLX since it needs  */
                                               /* to be converted in HLL_DIVB_SOUCE */
    #ifdef GLM_MHD
     VL = ARRAY_2D(NMAX_POINT, NVAR, double);
     VR = ARRAY_2D(NMAX_POINT, NVAR, double);
     UL = ARRAY_2D(NMAX_POINT, NVAR, double);
     UR = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
  }

  #ifdef GLM_MHD
   GLM_Solve (state, VL, VR, beg, end, grid);
   PrimToCons (VL, UL, beg, end);
   PrimToCons (VR, UR, beg, end);
  #else
   VL = state->vL; UL = state->uL;
   VR = state->vR; UR = state->uR;
  #endif

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (VL, a2L, hL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, hR, beg, end, FACE_CENTER, grid);

  Flux (UL, VL, hL, fL, pL, beg, end);
  Flux (UR, VR, hR, fR, pR, beg, end);

/* ----------------------------------------
     get max and min signal velocities
   ---------------------------------------- */

  SL = state->SL; SR = state->SR;
  HLL_Speed (VL, VR, a2L, a2R, hL, hR, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    uL = UL[i];  uR = UR[i];

/* ---- Compute fluxes ---- */
    
    bmin = MIN(0.0, SL[i]);
    bmax = MAX(0.0, SR[i]);
    scrh = 1.0/(bmax - bmin);
    for (nv = NFLX; nv--; ){
      state->flux[i][nv]  = bmin*bmax*(uR[nv] - uL[nv])
                         +  bmax*fL[i][nv] - bmin*fR[i][nv];
      state->flux[i][nv] *= scrh;
    }
    state->press[i] = (bmax*pL[i] - bmin*pR[i])*scrh;
  }

/* --------------------------------------------------------
              initialize source term
   -------------------------------------------------------- */
 
  #if MHD_FORMULATION == EIGHT_WAVES
/*
   POWELL_DIVB_SOURCE (state, beg, end, grid);
*/

  /* ----------------------------------------------------
       to avoid conversion problems in HLL_DIVB_SOURCE, 
       we use the HLL average provided by SR = -SL = 1 
     ---------------------------------------------------- */
/*
   for (i = beg; i <= end; i++) {
     uL = state->uL[i]; uR = state->uR[i];
     for (nv = 0; nv < NFLX; nv++) {
       Uhll[i][nv] = 0.5*(uR[nv] + uL[nv] + fL[i][nv] - fR[i][nv]);
     }
     Uhll[i][MXn] += (pL[i] - pR[i])*0.5;
     for (nv = NFLX; nv < NVAR; nv++) Uhll[i][nv] = 0.0;
   }
*/
   HLL_DIVB_SOURCE (state, Uhll, beg + 1, end, grid);
  #endif

}
