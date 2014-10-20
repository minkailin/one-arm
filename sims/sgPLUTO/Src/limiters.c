/* *******************************************************************
 
   FILE: limiters.c

   PURPOSE: Compute slope limiters needed in the interpolation 
            routines; slopes are returned as

                           dV|
            dV_LIM[i]  =   --|  Delta x_i 
                           dx|i 
  
            where dV/dx depends on the choice of the limiter and the 
            grid type (uniform or not).


   ARGUMENTS:

     dwp (in)    : array w/forward  undivided differences 
     dwm (in)    : array w/backward undivided differences
     dw_lim(out) : array of limited slopes
     ibeg (in)   : starting point
     iend (in)   : ending point
     GG (in)     : grid structure in the current sweep direction 


   NOTE:   In developing new slopes limiters, keep in mind that 
           in order to satisfy the TVD condition you must have
           (see TORO, pg 477):

        min(U[i], U[i+1]) <= UR = U[i] + dV_LIM[i] <= max (U[i], U[i+1])

        min(U[i], U[i-1]) <= UL = U[i] - dV_LIM[i] <= max (U[i], U[i-1])
 
   ******************************************************************* */

#include"pluto.h"


/* ********************************************************************* */
inline double Lflat_lim (const double dp, const double dm)
/*! 
 * 
 *   Zero-slope; first order interpolation
 *
 *********************************************************************** */
{
  return 0.0;
}


/* ********************************************************************* */
inline double Lminmod_lim (const double dp, const double dm)
/*! 
 *
 * Evaluate slopes using the minmod limiter:
 *
 *
 *                   xm  when |xm| < |xp|  AND xp*xm > 0
 *  minmod(xp, xm) = xp  when |xp| < |xm|  AND xp*xm > 0
 *                   0   otherwise
 *
 *  Here xp and xm are approximations to forward and backward first 
 *  derivative.
 * 
 *********************************************************************** */
{
  return (dp*dm > 0.0 ? (fabs(dp) < fabs(dm) ? dp:dm):0.0);
}

/* ********************************************************************* */
inline double Lvanalbada_lim (const double dp, const double dm)
/*!
 *   va(x,y) = x*y*(x + y)/(x*x + y*y)
 *
 *
 * Notice although this limiter does not require the switch x*y > 0,
 * we nevertheless employ it.
 * 
 *********************************************************************** */
#define EPS_VA  1.e-18
{
  double dp2, dm2;

  if (dp*dm > 0.0){
    dp2 = dp*dp;
    dm2 = dm*dm;
    return (dp*(dm2 + EPS_VA) + dm*(dp2 + EPS_VA))/(dp2 + dm2 + EPS_VA);
  }else{
    return 0.0;
  }
}
#undef EPS_VA

      
/* ********************************************************************* */
inline double Lvanleer_lim (const double dp, const double dm)
/*! 
 *
 *    vl = 2*x*y/(x + y)   when x*y > 0
 *    vl = 0               otherwise
 *
 *********************************************************************** */
{
  double s;
  s = dp*dm;
  return (s > 0.0 ? 2.0*s/(dp + dm): 0.0);
}

/* ********************************************************************* */
inline double Lmc_lim (const double dp, const double dm, const double dc)
/*! 
 * Evaluate slopes using the monotonized central difference limiter:
 *
 *      mc_lim (xp, xm, xc) = MinMod [2*xp , 2*xm , xc]
 *
 *  Where xp, xm and xc are, respectively, approximations to 
 *  the forward, backward and central derivatives.
 *
 *********************************************************************** */
{
  double d2;

  if (dp*dm < 0.0) return 0.0;
  d2 = 2.0*(fabs(dp) < fabs(dm) ? dp:dm);
  return (fabs(d2) < fabs(dc) ? d2:dc);
}

/* ********************************************************************* */
inline double Lgminmod_lim (const double dp, const double dm, const double dc)
/*! 
 *
 * Evaluate slopes using the general mindmod limiter:
 *
 *      gminmod (xp, xm, xc) = MinMod [a*xp , a*xm , xc]
 *
 *  Where xp, xm and xc are, respectively, approximations to 
 *  the forward, backward and central derivatives, and 
 *  1 < a < 2.
 *
 * setting a = 1  yields  the minmod limiter
 * setting a = 2  yields  the mc     limiter
 *
 *********************************************************************** */
{
  int i;
  double scrh, a = 1.5;

  if (dp*dm > 0.0){
    scrh   = a*ABS_MIN(dp,dm);
    return ABS_MIN(scrh, dc);
  }
  return 0.0;
}

/* ********************************************************************* */
inline double Lumist_lim (const double dp, const double dm)
/*
 *
 *
 * Evaluate slopes using the Umist Limiter :
 *
 *   UMIST (x,y) = MinMod [2*xp , 2*xm , 1/4 (xp + 3xm) , 1/4 (xm + 3xp)]
 *
 *
 *********************************************************************** */
{
  double x1,x2, scrh;

  if (dp*dm > 0.0) {
      x1 = 0.25*(dp + 3.0*dm);
      x2 = 0.25*(dm + 3.0*dp);
      scrh = 2.0*ABS_MIN(dp, dm);
      scrh = ABS_MIN(scrh, x1);
      return ABS_MIN(scrh, x2);
    } else {
      return 0.0;
    }
  }

