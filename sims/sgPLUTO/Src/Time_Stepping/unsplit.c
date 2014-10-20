/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advance equations using directionally-unsplit integrators 
         (method of lines).

  Main driver for RK unsplit integrations (DIMENSIONAL_SPLITTING == NO) 
  and for finite difference methods (RK3).
  Time stepping include Euler, RK2 and RK3.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date   Oct 1, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
int Unsplit (const Data *d, Riemann_Solver *Riemann, 
             Time_Step *Dts, Grid *grid)
/*!
 * Advance the equations by a single time step using unsplit 
 * integrators based on the method of lines.
 *
 * \param [in,out]      d  pointer to Data structure
 * \param [in]    Riemann  pointer to a Riemann solver function
 * \param [in,out]    Dts  pointer to time step structure
 * \param [in]       grid  pointer to array of Grid structures
 *    
 *********************************************************************** */
{
  int  ii, jj, kk;
  int  *i, *j, *k;
  int  in, nv;
  
  static double  one_third = 1.0/3.0;
  static State_1D state;
  static Data_Arr UU, UU_1;
  double dt;
  double *inv_dl, dl2;
  static double ***C_dt[NVAR], **dcoeff;
  Index indx;

/* ----------------------------------------------------
                   Allocate memory 
   ---------------------------------------------------- */

  if (state.rhs == NULL){
    MakeState (&state);
    UU      = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    UU_1    = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);

    #if (PARABOLIC_FLUX & EXPLICIT)
     dcoeff = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
    
  /* -------------------------------------------------------
      C_dt is an array used to storee the inverse time step 
      for advection and diffusion.
      We use C_dt[RHO] for advection, 
             C_dt[MX1] for viscosity,
             C_dt[BX1...BX3] for resistivity and
             C_dt[ENG] for thermal conduction.
     ------------------------------------------------------- */
      
    C_dt[RHO] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #if VISCOSITY == EXPLICIT
     C_dt[MX1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
    #if RESISTIVE_MHD == EXPLICIT
     EXPAND(C_dt[BX1] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double); ,
            C_dt[BX2] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double); ,
            C_dt[BX3] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);)
    #endif
    #if THERMAL_CONDUCTION == EXPLICIT
     C_dt[ENG] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
   
  }

/* ----------------------------------------------------
    STEP 1: predictor (EULER, RK2, RK3)
   ---------------------------------------------------- */

  g_intStage = 1;  
  dt = g_dt;

  Boundary (d, ALL_DIR, grid);
  #ifdef FARGO
   FARGO_SubtractVelocity (d,grid); 
  #endif
  #if SHOCK_FLATTENING == MULTID
   FindShock (d, grid);
  #endif
  
/* ----------------------------------------------------
    Convert primitive to conservative and reset arrays
   ---------------------------------------------------- */
     
  KTOT_LOOP(kk) JTOT_LOOP(jj){
    ITOT_LOOP(ii){
      for (nv = NVAR; nv--;  ) state.v[ii][nv] = d->Vc[nv][kk][jj][ii];
    }
    PrimToCons(state.v, UU[kk][jj], 0, NX1_TOT-1);
    memcpy ((void *)UU_1[kk][jj][0], UU[kk][jj][0], NX1_TOT*NVAR*sizeof(double));

    memset ((void *)C_dt[0][kk][jj],'\0', NX1_TOT*sizeof(double));
    #if VISCOSITY == EXPLICIT
     memset ((void *)C_dt[MX1][kk][jj],'\0', NX1_TOT*sizeof(double));
     C_dt[MX1][kk][jj][ii] = 0.0;
    #endif
    #if RESISTIVE_MHD == EXPLICIT
     EXPAND(memset ((void *)C_dt[BX1][kk][jj],'\0', NX1_TOT*sizeof(double));  ,
            memset ((void *)C_dt[BX2][kk][jj],'\0', NX1_TOT*sizeof(double));  ,
            memset ((void *)C_dt[BX3][kk][jj],'\0', NX1_TOT*sizeof(double));)
    #endif
    #if THERMAL_CONDUCTION == EXPLICIT
     memset ((void *)C_dt[ENG][kk][jj],'\0', NX1_TOT*sizeof(double));  
    #endif
    
/*    
    ITOT_LOOP(ii){
      for (nv = NVAR; nv--;  ) UU_1[kk][jj][ii][nv] = UU[kk][jj][ii][nv];
      
      C_dt[0][kk][jj][ii] = 0.0;
      #if VISCOSITY == EXPLICIT
       C_dt[MX1][kk][jj][ii] = 0.0;
      #endif
      #if RESISTIVE_MHD == EXPLICIT
       EXPAND(C_dt[BX1][*k][*j][*i] = 0.0;  ,
              C_dt[BX2][*k][*j][*i] = 0.0;  ,
              C_dt[BX3][*k][*j][*i] = 0.0;)
      #endif
      #if THERMAL_CONDUCTION == EXPLICIT
       C_dt[ENG][*k][*j][*i] = 0.0;
      #endif
    }
*/
  }

/* ----------------------------------
         Loop on dimensions
   ---------------------------------- */

  for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){
  
    SetIndexes (&indx, grid);  /* -- set normal and transverse indices -- */
    ResetState (d, &state, grid);
    TRANSVERSE_LOOP(indx,in,i,j,k){  

      inv_dl = GetInverse_dl(grid);

      for (in = 0; in < indx.ntot; in++) {
        for (nv = NVAR; nv--;  ) state.v[in][nv] = d->Vc[nv][*k][*j][*i];
/*        TestCopy(d->Vc[0][*k][*j], state.v[in], *i);  */
        
        #ifdef STAGGERED_MHD
         state.bn[in] = d->Vs[g_dir][*k][*j][*i];
        #endif
      }
      CheckNaN (state.v, 0, indx.ntot-1,0);
      States  (&state, indx.beg - 1, indx.end + 1, grid); 
      Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
      #ifdef STAGGERED_MHD
       CT_StoreEMF (&state, indx.beg - 1, indx.end, grid);
      #endif
      #if (PARABOLIC_FLUX & EXPLICIT)
       ParabolicFlux(d->Vc, &state, dcoeff, indx.beg - 1, indx.end, grid);
      #endif
      #if UPDATE_VECTOR_POTENTIAL == YES
       VectorPotentialUpdate (d, NULL, &state, grid);
      #endif
      #ifdef SHEARINGBOX
       SB_SaveFluxes (&state, grid);
      #endif

      RightHandSide (&state, Dts, indx.beg, indx.end, dt, grid);
      for (in = indx.beg; in <= indx.end; in++) { 
        #if !GET_MAX_DT
         C_dt[0][*k][*j][*i] += 0.5*(Dts->cmax[in-1] + Dts->cmax[in])*inv_dl[in];
        #endif
        #if VISCOSITY == EXPLICIT
         dl2 = 0.5*inv_dl[in]*inv_dl[in];
         C_dt[MX1][*k][*j][*i] += (dcoeff[in][MX1]+dcoeff[in-1][MX1])*dl2;
        #endif
        #if RESISTIVE_MHD == EXPLICIT
         dl2 = 0.5*inv_dl[in]*inv_dl[in];
         EXPAND(C_dt[BX1][*k][*j][*i] += (dcoeff[in-1][BX1]+dcoeff[in][BX1])*dl2;  ,
                C_dt[BX2][*k][*j][*i] += (dcoeff[in-1][BX2]+dcoeff[in][BX2])*dl2;  ,
                C_dt[BX3][*k][*j][*i] += (dcoeff[in-1][BX3]+dcoeff[in][BX3])*dl2;)
        #endif
        #if THERMAL_CONDUCTION  == EXPLICIT
         dl2 = 0.5*inv_dl[in]*inv_dl[in];
         C_dt[ENG][*k][*j][*i] += (dcoeff[in-1][ENG] + dcoeff[in][ENG])*dl2;
        #endif
        for (nv = NVAR; nv--;  )  UU_1[*k][*j][*i][nv] += state.rhs[in][nv];
      }
    }
  }

  #ifdef SHEARINGBOX 
   SB_CorrectFluxes (UU_1, g_time, g_dt, grid);
  #endif

  #ifdef STAGGERED_MHD
   CT_Update(d, grid);
   CT_AverageMagneticField (d->Vs, UU_1, grid);
  #endif

/* ----------------------------------------------------
    STEP 2: Corrector (RK2, RK3)
   ---------------------------------------------------- */

  #if TIME_STEPPING != EULER

  /* -------------------------
       convert to primitive
     ------------------------- */
     
   g_dir = IDIR;
   SetIndexes (&indx, grid);
   for (kk = KBEG; kk <= KEND; kk++){ g_k = &kk;
   for (jj = JBEG; jj <= JEND; jj++){ g_j = &jj;
     ConsToPrim (UU_1[kk][jj], state.v, IBEG, IEND, state.flag);
     for (ii = IBEG; ii <= IEND; ii++){
     for (nv = 0; nv < NVAR; nv++) {
       d->Vc[nv][kk][jj][ii] = state.v[ii][nv];
     }} 
   }}

   g_intStage = 2;
   #ifdef FARGO
    FARGO_AddVelocity (d,grid); 
   #endif
   Boundary (d, ALL_DIR, grid);
   #ifdef FARGO
    FARGO_SubtractVelocity (d,grid); 
   #endif
   #if SHOCK_FLATTENING == MULTID
    FindShock (d, grid);
   #endif

   DOM_LOOP(kk, jj, ii){
   for (nv = NVAR; nv--;  ) {
     #if TIME_STEPPING == RK2 
      UU_1[kk][jj][ii][nv] = 0.5*(UU[kk][jj][ii][nv] + UU_1[kk][jj][ii][nv]);
     #elif TIME_STEPPING == RK3 
      UU_1[kk][jj][ii][nv] = 0.75*UU[kk][jj][ii][nv] + 0.25*UU_1[kk][jj][ii][nv];
     #endif
   }}

   #if TIME_STEPPING == RK2 
    dt = 0.5*g_dt;
   #elif TIME_STEPPING == RK3
    dt = 0.25*g_dt;
   #endif

/* ----------------------------------
         Loop on dimensions
   ---------------------------------- */

   for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

     SetIndexes (&indx, grid);    /* -- set normal and transverse indices -- */
     ResetState (d, &state, grid);
     TRANSVERSE_LOOP(indx,in,i,j,k){  

       for (in = 0; in < indx.ntot; in++) {
         for (nv = NVAR; nv--;  ) state.v[in][nv] = d->Vc[nv][*k][*j][*i];
         #ifdef STAGGERED_MHD
          state.bn[in] = d->Vs[g_dir][*k][*j][*i];
         #endif
       }

       States  (&state, indx.beg - 1, indx.end + 1, grid);     
       Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
       #ifdef STAGGERED_MHD
        CT_StoreEMF (&state, indx.beg - 1, indx.end, grid);
       #endif

       #if (PARABOLIC_FLUX & EXPLICIT)
        ParabolicFlux (d->Vc, &state, dcoeff,  indx.beg - 1, indx.end, grid);
       #endif

       #if UPDATE_VECTOR_POTENTIAL == YES
        VectorPotentialUpdate (d, NULL, &state, grid);
       #endif
       #ifdef SHEARINGBOX
        SB_SaveFluxes (&state, grid);
       #endif
       RightHandSide (&state, Dts, indx.beg, indx.end, dt, grid);
      
       for (in = indx.beg; in <= indx.end; in++) {
       for (nv = NVAR; nv--;  ) {
         UU_1[*k][*j][*i][nv] += state.rhs[in][nv];
       }}
     }
   }

   #ifdef SHEARINGBOX
    SB_CorrectFluxes (UU_1, g_time+g_dt, dt, grid);
   #endif

   #ifdef STAGGERED_MHD
    CT_Update(d, grid);
    CT_AverageMagneticField (d->Vs, UU_1, grid);
   #endif

  #endif  /* TIME_STEPPING != EULER */

/* ----------------------------------------------------
    STEP 3: final corrector for third-order algorithms
   ---------------------------------------------------- */

  #if TIME_STEPPING == RK3

  /* -------------------------
       convert to primitive
     ------------------------- */
     
   g_dir = IDIR;
   SetIndexes (&indx, grid);
   for (kk = KBEG; kk <= KEND; kk++){ g_k = &kk;
   for (jj = JBEG; jj <= JEND; jj++){ g_j = &jj;
     ConsToPrim (UU_1[kk][jj], state.v, IBEG, IEND, state.flag);
     for (ii = IBEG; ii <= IEND; ii++){
     for (nv = 0; nv < NVAR; nv++) {
       d->Vc[nv][kk][jj][ii] = state.v[ii][nv];
     }} 
   }}

   g_intStage = 3;
   #ifdef FARGO
    FARGO_AddVelocity (d,grid); 
   #endif
   Boundary (d, ALL_DIR, grid);
   #ifdef FARGO
    FARGO_SubtractVelocity (d,grid); 
   #endif
   #if SHOCK_FLATTENING == MULTID
    FindShock (d, grid);
   #endif

   DOM_LOOP(kk,jj,ii){
   for (nv = NVAR; nv--;  ) {
     UU_1[kk][jj][ii][nv] = one_third*(    UU  [kk][jj][ii][nv] +
                                       2.0*UU_1[kk][jj][ii][nv]);
   }}

   dt = 2.0/3.0*g_dt;

/* ----------------------------------
         Loop on dimensions
   ---------------------------------- */

   for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

     SetIndexes (&indx, grid);  /* -- set normal and transverse indices -- */
     ResetState (d, &state, grid);
     TRANSVERSE_LOOP(indx,in,i,j,k){  

       for (in = 0; in < indx.ntot; in++) {
         for (nv = NVAR; nv--;  ) state.v[in][nv] = d->Vc[nv][*k][*j][*i];
         #ifdef STAGGERED_MHD
          state.bn[in] = d->Vs[g_dir][*k][*j][*i];
         #endif
       }

       States  (&state, indx.beg - 1, indx.end + 1, grid);
       Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
       #ifdef STAGGERED_MHD
        CT_StoreEMF (&state, indx.beg - 1, indx.end, grid);
       #endif

       #if (PARABOLIC_FLUX & EXPLICIT) 
        ParabolicFlux (d->Vc, &state, dcoeff,  indx.beg - 1, indx.end, grid);
       #endif

       #if UPDATE_VECTOR_POTENTIAL == YES
        VectorPotentialUpdate (d, NULL, &state, grid);
       #endif
       #ifdef SHEARINGBOX
        SB_SaveFluxes (&state, grid);
       #endif
       RightHandSide (&state, Dts, indx.beg, indx.end, dt, grid);

       for (in = indx.beg; in <= indx.end; in++) {
       for (nv = NVAR; nv--;  ) {
         UU_1[*k][*j][*i][nv] += state.rhs[in][nv];
       }}
     }
   }

   #ifdef SHEARINGBOX
    SB_CorrectFluxes (UU_1, g_time+g_dt, dt, grid);
   #endif

   #ifdef STAGGERED_MHD
    CT_Update(d, grid);
    CT_AverageMagneticField (d->Vs, UU_1, grid);
   #endif

  #endif /* TIME_STEPPING == RK3 */

/* ----------------------------------------------
     Shift solution with FARGO scheme
   ---------------------------------------------- */

  #ifdef FARGO
   FARGO_ShiftSolution (UU_1, d->Vs, grid);
  #endif 

/* ----------------------------------------------
    convert to primitive and get maximum of 
    inverse dt coefficients
   ---------------------------------------------- */
   
  g_dir = IDIR;
  SetIndexes (&indx, grid);
  for (kk = KBEG; kk <= KEND; kk++){ g_k = &kk;
  for (jj = JBEG; jj <= JEND; jj++){ g_j = &jj;
    ConsToPrim (UU_1[kk][jj], state.v, IBEG, IEND, state.flag);
    for (ii = IBEG; ii <= IEND; ii++){
      for (nv = 0; nv < NVAR; nv++) d->Vc[nv][kk][jj][ii] = state.v[ii][nv];

      #if !GET_MAX_DT
       Dts->inv_dta = MAX(Dts->inv_dta, C_dt[0][kk][jj][ii]);
      #endif
      #if VISCOSITY == EXPLICIT
       Dts->inv_dtp = MAX(Dts->inv_dtp, C_dt[MX1][kk][jj][ii]);
      #endif
      #if RESISTIVE_MHD == EXPLICIT
       EXPAND(Dts->inv_dtp = MAX(Dts->inv_dtp, C_dt[BX1][kk][jj][ii]);  ,
              Dts->inv_dtp = MAX(Dts->inv_dtp, C_dt[BX2][kk][jj][ii]);  ,
              Dts->inv_dtp = MAX(Dts->inv_dtp, C_dt[BX3][kk][jj][ii]);)
      #endif
      #if THERMAL_CONDUCTION == EXPLICIT
       Dts->inv_dtp = MAX(Dts->inv_dtp, C_dt[ENG][kk][jj][ii]);
      #endif
    }
  }}
  #if !GET_MAX_DT
   Dts->inv_dta /= (double)DIMENSIONS;
  #endif
  #if (PARABOLIC_FLUX & EXPLICIT)
   Dts->inv_dtp /= (double)DIMENSIONS;
  #endif

  #ifdef FARGO
   FARGO_AddVelocity (d,grid); 
  #endif
  return(0); /* -- step has been achieved, return success -- */
}

/*
void TestCopy(double *V, double *v1d, int i)
{
  long int nsize;
  int nv;
  
  nsize = NX1_TOT*NX2_TOT*NX3_TOT;
  for (V += i, nv = 0; nv < NVAR; V += nsize, nv++) {
    v1d[nv] = *V;
  }
}
*/ 
 