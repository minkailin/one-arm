#include "pluto.h"

/* **************************************************************************  */
void SplitSource (const Data *d, double dt, Time_Step *Dts, Grid *grid)
/*
 *
 * PURPOSE
 *
 *   Main driver for handling source terms as a separate
 *   step using operator splitting.
 *
 *   Source terms may be:
 *
 *    - optically thin radiative losses (cooling)
 *    - Diffusion operators: 
 *       * resistivity 
 *       * Thermal conduction
 *       * Viscosity
 *
 *
 *
 ***************************************************************************** */
{
  int i, j, k, nv;
  static real **v;
  real t_save;

/*  ---- GLM source term treated in main ----  */

/*
  #ifdef GLM_MHD
   GLM_SOURCE (d->Vc, dt, grid);
  #endif
*/

/*  ---------------------------------------------
             Cooling/Heating losses
    ---------------------------------------------  */

  #if COOLING != NO
   #if COOLING == POWER_LAW  /* -- solve exactly -- */
    PowerLawCooling (d->Vc, dt, Dts, grid);
   #else
    CoolingSource (d, dt, Dts, grid);
   #endif
  #endif

/* ----------------------------------------------
     Parabolic terms using STS:

       - resistivity 
       - thermal conduction
       - viscosity 
   ---------------------------------------------- */

  #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
   STS (d, Dts, grid);
  #endif

  #if (PARABOLIC_FLUX & RK_CHEBYSHEV)
   RKC (d, Dts, grid);
  #endif
                                                                                                                                                                             
}

