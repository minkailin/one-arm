/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advance equations using a directionally-split method.

  This is the main driver for dimensionally split integrations
  (DIMENSIONAL_SPLITTING = YES).
  Equations are advanced in time by taking contribution only from the 
  direction defined by the global variable ::g_dir.
  A full step requires as many calls as the number of DIMENSIONS.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)\n
           T. Matsakos
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ------------------------------------------------
     Need to define the extra array Vres ? 
   ------------------------------------------------ */

#if RESISTIVE_MHD == EXPLICIT
 #define BOUND_DIR ALL_DIR
#else
 #define BOUND_DIR (ALL_DIR)
#endif

/* ********************************************************************* */
int Sweep (const Data *d, Riemann_Solver *Riemann, 
           Time_Step *Dts, Grid *grid)
/*!
 * Advance the equations by incuding contribution from one direction
 * at a time.
 *
 * \param [in] d       pointer to Data structure
 * \param [in] Riemann pointer to a Riemann solver function
 * \param [in] Dts     pointer to a Time_Step structure
 * \param [in] grid    pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int  ii, jj, kk;
  int  *i, *j, *k;
  int  in, nv;
  Index indx;
  double dt, dl2, *inv_dl;
  static Data_Arr UU;
  static State_1D state;
  static double one_third = 1.0/3.0, **dcoeff;
  #if (PARABOLIC_FLUX & EXPLICIT)
   static Data_Arr Vres;
  #endif
static double **u;
  
/* -----------------------------------------------------------------
               Check algorithm compatibilities
   ----------------------------------------------------------------- */

  #ifdef STAGGERED_MHD
   print1 ("! CT works with unsplit integrators only\n");
   QUIT_PLUTO(1);
  #endif

/* -----------------------------------------------------------------
                      Allocate memory
   ----------------------------------------------------------------- */

  if (state.rhs == NULL){
    MakeState (&state);
u = ARRAY_2D(NMAX_POINT, NVAR, double);
     UU = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    #if (PARABOLIC_FLUX & EXPLICIT)
     Vres = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
     dcoeff = ARRAY_2D(NMAX_POINT, NVAR, double);  /* -- diffusion coefficient
                                                         array -- */ 
    #endif
  }

/* ----------------------------------------------------
                   STEP I   (PREDICTOR)
   ---------------------------------------------------- */

  g_intStage = 1;
  dt = g_dt;

  Boundary (d, BOUND_DIR, grid);
  #if SHOCK_FLATTENING == MULTID
   FindShock (d, grid);
  #endif

  #if (PARABOLIC_FLUX & EXPLICIT)
   for (nv = 0; nv < NVAR; nv++){
   TOT_LOOP(kk,jj,ii){
     Vres[nv][kk][jj][ii] = d->Vc[nv][kk][jj][ii];
   }}
  #endif

/* -----------------------------------------
    Convert primitive to conservative 
   ----------------------------------------- */
   
  KTOT_LOOP(kk) JTOT_LOOP(jj){
    ITOT_LOOP(ii){
      for (nv = NVAR; nv--;  ) state.v[ii][nv] = d->Vc[nv][kk][jj][ii];
    }
    PrimToCons(state.v, UU[kk][jj], 0, NX1_TOT-1);
  }

/* -------------------------------------------------
              Integration Loop
   ------------------------------------------------- */

  SetIndexes (&indx, grid);
  ResetState (d, &state, grid);
  TRANSVERSE_LOOP(indx,in,i,j,k){  

    inv_dl = GetInverse_dl(grid);

    for (in = 0; in < indx.ntot; in++) {
    for (nv = NVAR; nv--;  ) {
      state.v[in][nv] = d->Vc[nv][*k][*j][*i];
    }}

PrimToCons (state.v, u, 0, indx.ntot - 1);
    CheckNaN (state.v, 0, indx.ntot - 1, 0);
    
    #ifndef SINGLE_STEP
     for (in = indx.beg - 1; in <= indx.end + 1; in++) {
     for (nv = NVAR; nv--;  ) {
       UU[*k][*j][*i][nv] = u[in][nv]; 
     }}
    #endif

    States  (&state, indx.beg - 1, indx.end + 1, grid);
    Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);

    #if (PARABOLIC_FLUX & EXPLICIT)/* !! will be first order in time 
                                                 for single step algorithm   !! */
     ParabolicFlux (Vres, &state, dcoeff, indx.beg - 1, indx.end, grid);
    #endif
    #if UPDATE_VECTOR_POTENTIAL == YES
     VectorPotentialUpdate (d, NULL, &state, grid);
    #endif

    RightHandSide (&state, Dts, indx.beg, indx.end, dt, grid);

    for (in = indx.beg; in <= indx.end; in++) {
#if !GET_MAX_DT
      Dts->inv_dta = MAX(Dts->inv_dta, Dts->cmax[in]*inv_dl[in]);
#endif
      #if VISCOSITY == EXPLICIT
       dl2 = inv_dl[in]*inv_dl[in];
       Dts->inv_dtp = MAX(Dts->inv_dtp, dcoeff[in][MX1]*dl2);
      #endif
      #if RESISTIVE_MHD == EXPLICIT
       dl2 = inv_dl[in]*inv_dl[in];
       EXPAND(Dts->inv_dtp = MAX(Dts->inv_dtp, dcoeff[in][BX1]*dl2); ,
              Dts->inv_dtp = MAX(Dts->inv_dtp, dcoeff[in][BX2]*dl2); ,
              Dts->inv_dtp = MAX(Dts->inv_dtp, dcoeff[in][BX3]*dl2);)
      #endif
      #if THERMAL_CONDUCTION == EXPLICIT
       dl2 = inv_dl[in]*inv_dl[in];
       Dts->inv_dtp = MAX(Dts->inv_dtp, dcoeff[in][ENG]*dl2);
      #endif
      for (nv = NVAR; nv--;  ) {
        u[in][nv] += state.rhs[in][nv];
      }
    }

    ConsToPrim (u, state.v, indx.beg, indx.end, state.flag);
    for (in = indx.beg; in <= indx.end; in++) {
    for (nv = NVAR; nv--;  ) {
      d->Vc[nv][*k][*j][*i] = state.v[in][nv];
    }}
  }

/* ----------------------------------------------------
                   STEP II  (or CORRECTOR)
   ---------------------------------------------------- */

#if (TIME_STEPPING == RK2) || (TIME_STEPPING == RK3)

  g_intStage = 2;
  Boundary (d, BOUND_DIR, grid);
  #if (PARABOLIC_FLUX & EXPLICIT)
   for (nv = 0; nv < NVAR; nv++){
   TOT_LOOP(kk,jj,ii){
     Vres[nv][kk][jj][ii] = d->Vc[nv][kk][jj][ii];
   }}
  #endif

  ResetState (d, &state, grid);
  TRANSVERSE_LOOP(indx,in,i,j,k){  

    for (in = 0; in < indx.ntot; in++) {
    for (nv = NVAR; nv--;  ) {
      state.v[in][nv] = d->Vc[nv][*k][*j][*i];
    }}

    PrimToCons (state.v, u, 0, indx.ntot - 1);
    States  (&state, indx.beg - 1, indx.end + 1, grid);
    Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
    #if (PARABOLIC_FLUX & EXPLICIT)
     ParabolicFlux(Vres, &state, dcoeff, indx.beg - 1, indx.end, grid);
    #endif
    #if UPDATE_VECTOR_POTENTIAL == YES
     VectorPotentialUpdate (d, NULL, &state, grid);
    #endif

    RightHandSide (&state, Dts, indx.beg, indx.end, g_dt, grid);
    for (in = indx.beg; in <= indx.end; in++) {
    for (nv = NVAR; nv--;  ) {
      #if TIME_STEPPING == RK2
       u[in][nv] = 0.5*(UU[*k][*j][*i][nv] + u[in][nv] + state.rhs[in][nv]);
      #elif TIME_STEPPING == RK3    
       u[in][nv] = 0.25*(3.0*UU[*k][*j][*i][nv] + u[in][nv] + state.rhs[in][nv]);
      #endif
    }}

    ConsToPrim (u, state.v, indx.beg, indx.end, state.flag);

    for (in = indx.beg; in <= indx.end; in++) {
    for (nv = NVAR; nv--;  ) {
      d->Vc[nv][*k][*j][*i] = state.v[in][nv];
    }}
  }
#endif

/* ----------------------------------------------------
                   STEP III  (or CORRECTOR)
   ---------------------------------------------------- */

#if TIME_STEPPING == RK3 

  g_intStage = 3;
  Boundary (d, BOUND_DIR, grid);

  #if (PARABOLIC_FLUX & EXPLICIT)
   for (nv = 0; nv < NVAR; nv++){
   TOT_LOOP(kk,jj,ii){
      Vres[nv][kk][jj][ii] = d->Vc[nv][kk][jj][ii];
   }}
  #endif

  ResetState (d, &state, grid);
  TRANSVERSE_LOOP(indx,in,i,j,k){  

    for (in = 0; in < indx.ntot; in++) {
    for (nv = NVAR; nv--;  ) {
      state.v[in][nv] = d->Vc[nv][*k][*j][*i];
    }}
      
    PrimToCons (state.v, u, 0, indx.ntot - 1);
    States  (&state, indx.beg - 1, indx.end + 1, grid);
    Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
    #if (PARABOLIC_FLUX & EXPLICIT)
     ParabolicFlux (Vres, &state, dcoeff, indx.beg - 1, indx.end, grid);
    #endif
    #if UPDATE_VECTOR_POTENTIAL == YES
     VectorPotentialUpdate (d, NULL, &state, grid);
    #endif

    RightHandSide (&state, Dts, indx.beg, indx.end, g_dt, grid);

    for (in = indx.beg; in <= indx.end; in++) {
    for (nv = NVAR; nv--;  ) {
      u[in][nv] = one_third*(UU[*k][*j][*i][nv] + 
                        2.0*(u[in][nv] + state.rhs[in][nv]));
    }}    

    ConsToPrim (u, state.v, indx.beg, indx.end, state.flag);

    for (in = indx.beg; in <= indx.end; in++) {
    for (nv = NVAR; nv--;  ) {
       d->Vc[nv][*k][*j][*i] = state.v[in][nv];
    }}

  }
  
#endif

  return(0); /* -- step has been achieved, return success -- */
}

