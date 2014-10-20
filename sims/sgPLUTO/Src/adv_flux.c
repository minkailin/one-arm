/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute flux for passive scalars.                            

  Compute the interface upwind flux for passive scalars obeying 
  advection equations of the form:

  \f[ \partial_tq + v\cdot \nabla q = 0 
      \qquad\Longleftrightarrow\qquad
      \partial_t(\rho q) + \nabla\cdot(\rho\vec{v}q) = 0 \f]

  Fluxes are computed using an upwind selection rule based on the 
  density flux, already computed during a previous Riemann solver:
  \f[ 
    (\rho vq)_{i+\HALF} = \left\{\begin{array}{ll}
    (\rho v)_{i+\HALF}q_L & \;\textrm{if} 
           \quad (\rho v)_{i+\HALF} \ge 0 \\ \noalign{\medskip}
    (\rho v)_{i+\HALF}q_R & \; \textrm{otherwise}  
   \end{array}\right.
  \f]
  where \f$ (\rho v)_{i+\HALF}\f$ is the density flux computed with 
  the employed Riemann solver.

  \b Reference
  "The consistent multi-fluid advection method"
   Plewa and Muller, A&A (1999) 342, 179.
 
 \note
 The CMA (Consistent multi-fluid advection method), may be
 switched on (#define USE_CMA  YES) when the sum of several scalars 
 (e.g. different fluid phases) must be normalized to 1. 

  \author A. Mignone (mignone@ph.unito.it), O. Tesileanu
  \date   14 Aug, 2012
*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define USE_CMA NO

#ifndef C_IONS
 #define C_IONS 5
#endif

#ifndef N_IONS
 #define N_IONS 5
#endif

#ifndef O_IONS
 #define O_IONS 5
#endif

#ifndef Ne_IONS
 #define Ne_IONS 5
#endif

#ifndef S_IONS
 #define S_IONS 5
#endif

/* ********************************************************************* */
void AdvectFlux (const State_1D *state, int beg, int end, Grid *grid)
/*! 
 *
 * \param [in,out] state
 * \param [in]      beg    initial index of computation 
 * \param [in]      end    final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int    i, nv;
  real *ts, *flux, *vL, *vR;
  real s, rho;
  real phi;
  static double *sigma, **vi;
  
/* -- compute scalar's fluxes -- */

  for (i = beg; i <= end; i++){

    flux = state->flux[i];
    vL   = state->vL[i];
    vR   = state->vR[i];

    ts = flux[RHO] > 0.0 ? vL:vR;

    for (nv = NFLX; nv < (NFLX + NSCL); nv++){
      flux[nv] = flux[RHO]*ts[nv];
    }
    
    #if COOLING == MINEq

     /* -----   He   -----  */
                         
     phi = ts[HeI] + ts[HeII];
     for (nv = HeI; nv <= HeII; nv++) flux[nv] /= phi;

     /* -----   C   -----  */

     phi = 0.0;
     for (nv = CI; nv < CI+C_IONS; nv++) phi += ts[nv]; 
     for (nv = CI; nv < CI+C_IONS; nv++) flux[nv] /= phi;

     /* -----   N   -----  */

     phi = 0.0;
     for (nv = NI; nv < NI+N_IONS; nv++) phi += ts[nv]; 
     for (nv = NI; nv < NI+N_IONS; nv++) flux[nv] /= phi;

     /* -----   O   -----  */

     phi = 0.0;
     for (nv = OI; nv < OI+O_IONS; nv++) phi += ts[nv]; 
     for (nv = OI; nv < OI+O_IONS; nv++) flux[nv] /= phi;

     /* -----   Ne   -----  */
     phi = 0.0;
     for (nv = NeI; nv < NeI+Ne_IONS; nv++) phi += ts[nv];
     for (nv = NeI; nv < NeI+Ne_IONS; nv++) flux[nv] /= phi;

     /* -----   S   -----  */

     phi = 0.0;
     for (nv = SI; nv < SI+S_IONS; nv++) phi += ts[nv]; 
     for (nv = SI; nv < SI+S_IONS; nv++) flux[nv] /= phi;

    #endif

    #if USE_CMA == YES  /* -- only for tracers -- */
     phi = 0.0;
     for (nv = TR; nv < (TR+NTRACER); nv++) phi += ts[nv];
     for (nv = TR; nv < (TR+NTRACER); nv++) flux[nv] /= phi;
    #endif

    #if ENTROPY_SWITCH == YES
     if (flux[RHO] >= 0.0) flux[ENTR] = state->vL[i][ENTR]*flux[RHO];
     else                 flux[ENTR] = state->vR[i][ENTR]*flux[RHO];
    #endif
  }
}
