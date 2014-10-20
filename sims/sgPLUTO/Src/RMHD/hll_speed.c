#include"pluto.h"

/* *********************************************************************** */
void HLL_Speed (double **vL, double **vR, double *a2L, double *a2R,
                double *hL, double *hR, double *SL, double *SR, int beg, int end)
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
 *   vL (IN)       1-D array of left-edge  primitive values at i+1/2
 *   vR (IN)       1-D array of right-edge primitive values at i+1/2
 *   SL (OUT)      Array of left speeds
 *   SR (OUT)      Array of left speeds
 *   beg(IN)       initial grid index
 *   end(IN)       final grid index
 *
 *
 * LAST_MODIFIED
 *
 *   June 25, 2012 by Andrea Mignone  (mignone@ph.unito.it)
 *              
 *
 ************************************************************************** */
{
  int    i, err;
  static real *sl_min, *sl_max;
  static real *sr_min, *sr_max;

  if (sl_min == NULL){
    sl_min = ARRAY_1D(NMAX_POINT, double);
    sl_max = ARRAY_1D(NMAX_POINT, double);

    sr_min = ARRAY_1D(NMAX_POINT, double);
    sr_max = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------
              DAVIS Estimate  
   ---------------------------------------------- */

  MaxSignalSpeed (vL, a2L, hL, sl_min, sl_max, beg, end);
  MaxSignalSpeed (vR, a2R, hR, sr_min, sr_max, beg, end);
/*
  err = MAX_CH_SPEED (vL, sl_min, sl_max, beg, end);
  if (err != 0) return err;
  err = MAX_CH_SPEED (vR, sr_min, sr_max, beg, end);
  if (err != 0) return err;
*/
  for (i = beg; i <= end; i++) {
    SL[i] = MIN(sl_min[i], sr_min[i]);
    SR[i] = MAX(sl_max[i], sr_max[i]);
  }
/*  return 0; */
}

