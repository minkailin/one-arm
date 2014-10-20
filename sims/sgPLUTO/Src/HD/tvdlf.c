#include "pluto.h"

/* ********************************************************************* */
void LF_Solver (const State_1D *state, int beg, int end, 
                real *cmax, Grid *grid)
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
 *
 * LAST_MODIFIED
 *
 *   April 4th 2006, by Andrea Mignone  (mignone@to.astro.it)
 *
 *
 **************************************************************************** */
{
  int    nv, i;
  static real **fL, **fR, **vRL;
  static real *cRL_min, *cRL_max, *pL, *pR, *a2L, *a2R;
  double *uR, *uL;
  
  if (fR == NULL){
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    vRL = ARRAY_2D(NMAX_POINT, NFLX, double);
    cRL_min = ARRAY_1D(NMAX_POINT, double);
    cRL_max = ARRAY_1D(NMAX_POINT, double);
    pR  = ARRAY_1D(NMAX_POINT, double);
    pL  = ARRAY_1D(NMAX_POINT, double);
    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (state->uL, state->vL, a2L, fL, pL, beg, end);
  Flux (state->uR, state->vR, a2R, fR, pR, beg, end);

/* ---------------------------------------------------------------------
    Compute average state in order to get the local max characteristic
    velocity.
    In order to avoid underestimating this speed at strong reflecting
    boundaries (or when vxL = -vxR in general) , we average the
    absolute value of the normal velocity.
   ------------------------------------------------------------------- */

  for (i = beg; i <= end; i++)           {
    for (nv = NFLX; nv--;  ) {
      vRL[i][nv] = 0.5*(state->vL[i][nv] + state->vR[i][nv]);
    }
    vRL[i][VXn] = 0.5*(fabs(state->vL[i][VXn]) + fabs(state->vR[i][VXn]));
  }

/* ---------------------------------------------------------------------
           Compute max eigenvalue and flux
   --------------------------------------------------------------------- */

  SoundSpeed2    (vRL, a2R, NULL, beg, end, FACE_CENTER, grid);
  MaxSignalSpeed (vRL, a2R, cRL_min, cRL_max, beg, end);

  for (i = beg; i <= end; i++) {
    cRL_max[i] = MAX(fabs(cRL_max[i]), fabs(cRL_min[i]));
    cmax[i] = cRL_max[i];
    g_maxMach = MAX(g_maxMach, fabs(vRL[i][VXn])/sqrt(a2R[i]));
  }

  for (i = beg; i <= end; i++) {
    uR = state->uR[i];
    uL = state->uL[i];
    for (nv = NFLX; nv--;  ) {
      state->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] 
                                - cRL_max[i]*(uR[nv] - uL[nv]));
    }
    state->press[i] = 0.5*(pL[i] + pR[i]);
  }
}
