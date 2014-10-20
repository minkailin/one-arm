/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Flag zones where pressure is updated using entropy.

  Flag cells where pressure has to be updated from the conserved
  entropy density and compute entropy at the beginning of each
  integration stage.
  
  The flagging strategy is based on two switches designed to detect 
  the presence of compressive motion or shock waves in the fluid:
  - \f$ \nabla\cdot\vec{v} < 0 \f$
  - \f$ |\nabla p|/p > 0.1 \f$
  
  By default, if at least one of the switch is \e FALSE, we flag the 
  computational zone by turning the ::FLAG_ENTROPY bit of the 
  \c d->flag array on.
  In this way, the default strategy is to evolve the entropy equation
  everywhere, except at shocks (where both switches are \e TRUE).

  This flag will be checked later in the ConsToPrim() functions in 
  order to recover pressure from the entropy density rather than 
  from the total energy density.
  The update process is:
 
  - start with a vector of primitive variables  <tt> {V,s} </tt>
    where \c s is the entropy;
  - set boundary condition on \c {V}; 
    compute \c {s} from \c {V};
    flag zones where entropy may be used (flag = 1);
  - evolve equations for one time step;
  - convert <tt> {U,S} </tt> to primitive:
   
        where (flag == 0) {
          p = p(E);
        }else{  
          p = p(S);
          correct E using new p;
        }
                         
  \b Reference
     - "Maintaining Pressure Positivity in Magnetohydrodynamics Simulations"
       Balsara \& Spicer, JCP (1999) 148, 133
 

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#if (EOS != ISOTHERMAL) && (EOS != BAROTROPIC) && (ENTROPY_SWITCH == YES)

/* ********************************************************************* */
void EntropySwitch (const Data *d, Grid *grid)
/*!
 * Compute entropy and flag zones for updating pressure using the 
 * conserved entropy density.
 *
 * \param [in,out]  d   pointer to Data structure
 * \param [in]    grid  pointer to an array of Grid structures.
 *
 * \return This function has no return value.
 *
 *********************************************************************** */
{
  int i, j, k;
  int sw1, sw2;
  double dp, minpx, minpy, minpz, minp;
  double dv, scrh;
  double ***rho, ***pr, ***s;
  double ***vx, ***vy, ***vz;
  #if PHYSICS == MHD || PHYSICS == RMHD
   double ***bx, ***by, ***bz;
  #endif
  static double **v1d;

  if (v1d == NULL) v1d = ARRAY_2D(NMAX_POINT, NVAR, double);

  pr  = d->Vc[PRS];
  rho = d->Vc[RHO];
  EXPAND( vx = d->Vc[VX1];   ,
          vy = d->Vc[VX2];   , 
          vz = d->Vc[VX3]; )
  #if PHYSICS == MHD || PHYSICS == RMHD
   EXPAND( bx = d->Vc[BX1];   ,
           by = d->Vc[BX2];   , 
           bz = d->Vc[BX3]; )
  #endif
  s = d->Vc[ENTR];

/* -------------------------------------------------------------
                  Compute Entropy
   ------------------------------------------------------------- */

  KTOT_LOOP(k) {
  JTOT_LOOP(j) {
    ITOT_LOOP(i) {
      v1d[i][RHO] = rho[k][j][i];
      v1d[i][PRS] = pr[k][j][i];
    }
    Entropy(v1d, s[k][j], 0, NX1_TOT-1);
  }}

/* -------------------------------------------------------------
    Compute flags only at the beginning of the integration 
    process (g_intStage == 1).
   ------------------------------------------------------------- */

  if (g_intStage != 1) return;

  for (k = KOFFSET; k < NX3_TOT - KOFFSET; k++) { 
  for (j = JOFFSET; j < NX2_TOT - JOFFSET; j++) {
  for (i = IOFFSET; i < NX1_TOT - IOFFSET; i++) {

    D_EXPAND( dp  = fabs(pr[k][j][i + 1] - pr[k][j][i - 1]); ,
              dp += fabs(pr[k][j + 1][i] - pr[k][j - 1][i]); ,
              dp += fabs(pr[k + 1][j][i] - pr[k - 1][j][i]);  )

    D_EXPAND(minpx = MIN(pr[k][j][i - 1], pr[k][j][i + 1]);  ,
             minpy = MIN(pr[k][j - 1][i], pr[k][j + 1][i]);  ,
             minpz = MIN(pr[k - 1][j][i], pr[k + 1][j][i]);)

    D_EXPAND(minp = MIN(minpx, pr[k][j][i]);   ,
             minp = MIN(minp, minpy);          ,
             minp = MIN(minp, minpz);)

    dp /= 3.0;

/* -- flag regions with strong pressure gradients -- */

    sw1 = (dp/minp) > 0.1;

    D_EXPAND( dv  = vx[k][j][i + 1] - vx[k][j][i - 1]; ,
              dv += vy[k][j + 1][i] - vy[k][j - 1][i]; ,
              dv += vz[k + 1][j][i] - vz[k - 1][j][i];  )

/* -- flag regions with compressive motion -- */

    sw2 = dv < 0.0;

/* -- default: use entropy everywhere except at shocks -- */

    if (sw1 && sw2) {
      continue;
    }else{
      d->flag[k][j][i] |= FLAG_ENTROPY;
    }


/*    
    if (!(sw1 && sw2)){
      d->flag[k][j][i] |= FLAG_ENTROPY;
    }
*/
/*
    if (sw1 && sw2){
      int ip, jp, kp;
      for (jp = j-1; jp <= j+1; jp++){
      for (ip = i-1; ip <= i+1; ip++){
        d->flag[k][jp][ip] &= ~(FLAG_ENTROPY);
      }}
    }
*/

  }}}

}
#endif
