#include "pluto.h"

/* ***************************************************************** */
void PowerLawCooling (Data_Arr VV, real dt, Time_Step *Dts, Grid *grid)
/*
 *
 * PURPOSE:  
 *
 *    Take a source step to account for bremmstrahlung cooling
 *    term.
 *
 *
 *     We integrate the following ODE:
 *
 *      dp_cgs/dt_cgs = -(gamma - 1) Lambda(rho_cgs, T_cgs)
 *
 *     where:   Lambda = a_br/(mu*mH)^2 rho_cgs^2 * sqrt(T_cgs)
 *              a_br = 2.e-27   c.g.s
 *
 *     Here the subscript 'cgs' means that the respective 
 *     quantity is given in c.g.s units. 
 *     We denote with mu the molecular weight, while mH is the
 *     hydrogen mass (also in c.g.s).
 * 
 *     The previous equation is solved analytically, 
 *     since the density does not change during this step.
 *     We solve the  non-dimensional form:
 *
 *               dp/dt = -cost * rho * (rho*p)^(1/2) 
 *
 *     [notice that since  p/rho = T/KELVIN this is equivalent to:
 *      dp/dt = -cost rho^2 (T/KELVIN)^(1/2) ]
 *
 *    The quantity cost is determined by transforming
 *    the dimensional equation into the non-dimensional one.
 *    If p, rho and t are in code (non-dimensional) units and
 *    if L_0, rho_0, and V_0 are the unit length, density 
 *    and velocity, then cost is found to be:
 *
 *                  a_br * (gamma - 1) * L_0 * rho_0
 *       cost = -------------------------------------------
 *               sqrt(kB * mu * mH) * kB * mu * mH * V_0^2
 *
 *     where a_br = 2.e-27 (in c.g.s), kB is the Boltmann constant
 *
 *
 *
 ******************************************************************* */
{
  int   i, j, k;
  real  cost, dE;
  real  rho, p, T, p_f, T_f;

  cost  = g_unitLength*g_unitDensity/(g_unitVelocity*g_unitVelocity);
  cost *= 2.e-27*(g_gamma-1.0)/ (0.5*CONST_mH*sqrt(0.5*CONST_mH*CONST_kB));

/*  -------------------------------------------------------------
                Integrate analytically
    -------------------------------------------------------------  */

  dE = 1.e-18;
  DOM_LOOP(k,j,i){

/*  ----  Find initial temperature in Kelvin  ----  */

    rho = VV[RHO][k][j][i];
    p   = VV[PRS][k][j][i];

    T   = (p/rho*KELVIN);

    if (T < g_minCoolingTemp) continue;

/*  ----  Find final energy  ----  */

    p_f = sqrt(p) - 0.5*cost*rho*sqrt(rho)*dt;
    p_f = MAX(p_f, 0.0);
    p_f = p_f*p_f;
    T_f = p_f/rho*KELVIN;
    T_f = MAX (T_f, g_minCoolingTemp);

/*  ----  Update Energy  ----  */

    p_f = T_f*rho/KELVIN;

    VV[PRS][k][j][i] = p_f;
    dE = fabs(1.0 - p_f/p) + 1.e-18;

    Dts->dt_cool = MIN(Dts->dt_cool, dt*g_maxCoolingRate/dE);
  }

}


real MeanMolecularWeight (real *v)
{
  return(0.5);
}
