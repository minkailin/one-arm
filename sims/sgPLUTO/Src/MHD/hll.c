#include"pluto.h"

/* ************************************************************************* */
void HLL_Solver (const State_1D *state, int beg, int end, 
                 real *cmax, Grid *grid)
/* 
 *
 *
 * NAME
 *
 *   HLL_SOLVER
 *
 *
 * PURPOSE
 *
 *  - Solve riemann problem for the MHD equations using the 
 *    single state HLL Riemann solver
 * 
 *   - On input, it takes left and right primitive state
 *     vectors state->vL and state->vR at zone edge i+1/2;
 *     On output, return flux and pressure vectors at the
 *     same interface.
 *
 *   - Also, compute maximum wave propagation speed (cmax) 
 *     for  explicit time step computation
 *  
 *
 * LAST_MODIFIED
 *
 *   July 13th 2012, by Andrea Mignone  (mignone@ph.unito.it)
 *             
 *
 *************************************************************************** */
{
  int    nv, i, xdface;
  double scrh;
  double *uL, *uR, *SR, *SL;
  static double **fL, **fR, **Uhll;
  static double **VL, **VR, **UL, **UR;
  static double *pL, *pR, *a2L, *a2R;
  double **bgf;
    
  if (fL == NULL){
    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);

    pL  = ARRAY_1D(NMAX_POINT, double);
    pR  = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);
    a2R = ARRAY_1D(NMAX_POINT, double);

    Uhll = ARRAY_2D(NMAX_POINT, NFLX, double);
    #ifdef GLM_MHD
     VL = ARRAY_2D(NMAX_POINT, NVAR, double);
     VR = ARRAY_2D(NMAX_POINT, NVAR, double);
     UL = ARRAY_2D(NMAX_POINT, NVAR, double);
     UR = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
  }

  #if BACKGROUND_FIELD == YES
   bgf = GetBackgroundField (beg, end, FACE_CENTER, grid);
   #ifdef GLM_MHD 
    print1 ("! GLM MHD and BACKGROUND_FIELD incompatible\n");
    QUIT_PLUTO(1);
   #endif
  #endif

/* ------------------------------------------------
     Solve 2x2 Riemann problem with GLM Cleaning
   ------------------------------------------------ */

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

  SoundSpeed2 (VL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (UL, VL, a2L, bgf, fL, pL, beg, end);
  Flux (UR, VR, a2R, bgf, fR, pR, beg, end);

/* ----------------------------------------
     get max and min signal velocities
   ---------------------------------------- */
             
  SL = state->SL; SR = state->SR;
  HLL_Speed (VL, VR, a2L, a2R, bgf, SL, SR, beg, end);

/* ----------------------------------------
           compute HLL flux
   ---------------------------------------- */	     

  for (i = beg; i <= end; i++) {

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    if (SL[i] > 0.0){
    
      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = fL[i][nv];
      }
      state->press[i] = pL[i];
      
    }else if (SR[i] < 0.0){
    
      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = fR[i][nv];
      }
      state->press[i] = pR[i];
      
    }else{
    
      uL = UL[i]; uR = UR[i];

      scrh = 1.0 / (SR[i] - SL[i]);
    
      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                             SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        state->flux[i][nv] *= scrh;
      }
      state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
    }

  }


/* -----------------------------------------------------
               initialize source term
   ----------------------------------------------------- */

  #if MHD_FORMULATION == EIGHT_WAVES
/*
   ROE_DIVB_SOURCE (state, beg, end, grid);
*/
/*
   for (i = beg; i <= end; i++) {
     uR = state->uR[i]; uL = state->uL[i];
     scrh = 1.0 / (SR[i] - SL[i]);
     for (nv = 0; nv < NFLX; nv++) {
       Uhll[i][nv] = SR[i]*uR[nv] - SL[i]*uL[nv] +
                     fL[i][nv] - fR[i][nv];
       Uhll[i][nv] *= scrh;
     }
     Uhll[i][MXn] += (pL[i] - pR[i])*scrh;
   }
*/
   HLL_DIVB_SOURCE (state, Uhll, beg + 1, end, grid);

  #endif

}

