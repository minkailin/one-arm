#include"pluto.h"

/* ***************************************************************************** */
void HLL_Speed (double **vL, double **vR, double *a2L, double *a2R,
                double *SL, double *SR, int beg, int end)
/* 
 *
 *
 * NAME
 *
 *   HLL_SPEED
 *
 *
 * PURPOSE
 *
 *   Compute leftmost (SL) and rightmost (SR) speed for the Riemann fan
 * 
 *
 * ARGUMENTS
 *
 *   vL (IN)       1-D array of left-edge primitive values at i+1/2
 *   vR (IN)       1-D array of right-edge primitive values at i+1/2
 *   grid (IN)     Array of grids
 *   cmax(OUT)     Array of maximum characteristic speeds in this direction
 *
 *
 * LAST_MODIFIED
 *
 *   April 11/2008, Andrea Mignone  (mignone@to.astro.it)
 *              
 *
 ******************************************************************************** */
{
  int    i;
  static real *sl_min, *sl_max;
  static real *sr_min, *sr_max;

  if (sl_min == NULL){
    sl_min = ARRAY_1D(NMAX_POINT, double);
    sl_max = ARRAY_1D(NMAX_POINT, double);

    sr_min = ARRAY_1D(NMAX_POINT, double);
    sr_max = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------
    use Davis estimate for the signal velocities
   ---------------------------------------------- */

  MaxSignalSpeed (vL, a2L, sl_min, sl_max, beg, end);
  MaxSignalSpeed (vR, a2R, sr_min, sr_max, beg, end);
  for (i = beg; i <= end; i++) {
    SL[i] = MIN(sl_min[i], sr_min[i]);
    SR[i] = MAX(sl_max[i], sr_max[i]);
  }
}
