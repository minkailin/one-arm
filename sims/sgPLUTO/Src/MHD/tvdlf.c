#include "pluto.h"

/* ********************************************************************* */
void LF_Solver (const State_1D *state, int beg, int end, 
                double *cmax, Grid *grid)
/*
 *
 * NAME
 *
 *   LF_SOLVER
 *
 *
 * PURPOSE
 *
 *   - Solve Riemann problem using the Lax-Friedrichs 
 *     Rusanov solver:
 *
 *     F(i+1/2) = [FL + FR - c*(UR - UL)]*0.5
 *
 *     where c = max_speed[(UR + UL)*0.5]
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
 *   August 24, 2011 by Andrea Mignone  (mignone@ph.unito.it)
 *
 *
 **************************************************************************** */

{
  int     nv, i;
  double    cRL;
  double  *uR, *uL;
  static double **fL, **fR, **vRL;
  static double *pR, *pL, *cmin_RL, *cmax_RL, *a2L, *a2R;
  static double **VL, **VR, **UL, **UR;
  double **bgf;

  if (fR == NULL){
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    vRL = ARRAY_2D(NMAX_POINT, NFLX, double);
    cmin_RL = ARRAY_1D(NMAX_POINT, double);
    cmax_RL = ARRAY_1D(NMAX_POINT, double);
    pR      = ARRAY_1D(NMAX_POINT, double);
    pL      = ARRAY_1D(NMAX_POINT, double);

    a2R     = ARRAY_1D(NMAX_POINT, double);
    a2L     = ARRAY_1D(NMAX_POINT, double);
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

/* ------------------------------------------------------
            Compute max eigenvalue and fluxes
   ------------------------------------------------------ */

  for (i = beg; i <= end; i++) {
    for (nv = NFLX; nv--;      ) {
      vRL[i][nv] = 0.5*(VL[i][nv] + VR[i][nv]);
    }
    #if EOS == IDEAL
     g_maxMach = MAX(g_maxMach, fabs(vRL[i][VXn])/sqrt(g_gamma*vRL[i][PRS]/vRL[i][RHO])); 
    #elif EOS == BAROTROPIC
     cRL = sqrt(g_gamma*BAROTROPIC_PR(vRL[i][RHO])/vRL[i][RHO]);
     g_maxMach = MAX(g_maxMach, fabs(vRL[i][VXn])/cRL); 
    #elif EOS == ISOTHERMAL
     g_maxMach = MAX(g_maxMach, fabs(vRL[i][VXn])/g_isoSoundSpeed); 
    #endif
  }

  SoundSpeed2    (vRL, a2R, NULL, beg, end, FACE_CENTER, grid);
  MaxSignalSpeed (vRL, a2R, cmin_RL, cmax_RL, bgf, beg, end);

  for (i = beg; i <= end; i++) {
    cRL = MAX(fabs(cmin_RL[i]), fabs(cmax_RL[i]));
    state->SL[i] = -cRL, state->SR[i] = cRL;
    cmax[i] = cRL;
    uL = UL[i];
    uR = UR[i];
    for (nv = NFLX; nv--;    ) {
      state->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - cRL*(uR[nv] - uL[nv]));
    }
    state->press[i] = 0.5*(pL[i] + pR[i]);
  }

/* --------------------------------------------------------
              initialize source term
   -------------------------------------------------------- */

  #if MHD_FORMULATION == EIGHT_WAVES
   ROE_DIVB_SOURCE (state, beg + 1, end, grid);
  #endif

}

