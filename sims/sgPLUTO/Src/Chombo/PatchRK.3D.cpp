/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Implements a RK-midpoint time stepping routine.

  The RK-midpoing advance the system of conservation law by first
  evolving the solution at the half-time level and then doing a full
  step using midpoint quadrature rule:
  \f[ \left\{\begin{array}{lcl}
       U^{n+\HALF} &=& U^n + \frac{\Delta t}{2} R^n \\ \noalign{\medskip}
       U^{n+1}     &=& U^n + \Delta t R^{n+\HALF} 
      \end{array}\right.
  \f]
  
  \authors A. Mignone (mignone@ph.unito.it)\n
           C. Zanni   (zanni@oato.inaf.it)
  \date   Sep 20, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

/* ********************************************************************* */
void PatchPluto::updateSolution(FArrayBox&  a_U,
                                FArrayBox&  a_Utmp,
                                FArrayBox&  split_tags,
                                BaseFab<unsigned char>& a_Flags,
                                FluxBox&    a_F,
                                Time_Step   *Dts,
                                const Box&  UBox,
                                Grid *grid)
/*
 *
 *
 *
 *
 *********************************************************************** */
{
  CH_assert(isDefined());
  CH_assert(UBox == m_currentBox);

  int nv, in;
  int nxf, nyf, nzf, indf;
  int nxb, nyb, nzb;
  int *i, *j, *k;
  int ii, jj, kk;

  double ***UU[NVAR];
 #ifdef SKIP_SPLIT_CELLS
  double ***splitcells;
 #endif
  double *inv_dl, dl2;
  static Data d;
  static Data_Arr UU_1;
  static double ***C_dt[NVAR];
  static double **u;
  #if (PARABOLIC_FLUX & EXPLICIT)
   static double **dcoeff;
  #endif   
  Index indx;
  static State_1D state;

  Riemann_Solver *Riemann;
  Riemann = rsolver;

/* -----------------------------------------------------------------
               Check algorithm compatibilities
   ----------------------------------------------------------------- */

  #if !(GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL || GEOMETRY == SPHERICAL)
   print1 ("! RK only works in cartesian or cylindrical/spherical coordinates\n");
   QUIT_PLUTO(1);
  #endif     

  if (NX1_TOT > NMAX_POINT || NX2_TOT > NMAX_POINT || NX3_TOT > NMAX_POINT){
    print ("!updateSolution (RK_MID): need to re-allocate matrix\n");
    QUIT_PLUTO(1);
  }

/* -----------------------------------------------------------------
                          Allocate memory
   ----------------------------------------------------------------- */

  for (nv = 0; nv < NVAR; nv++){
    UU[nv] = ArrayMap(NX3_TOT, NX2_TOT, NX1_TOT, a_U.dataPtr(nv));
  }

  #ifdef SKIP_SPLIT_CELLS
   splitcells = ArrayBoxMap(KBEG, KEND, JBEG, JEND, IBEG, IEND, 
                                split_tags.dataPtr(0));
  #endif

/* -----------------------------------------------------------
         Allocate static memory areas
   -----------------------------------------------------------  */

  if (state.flux == NULL){

    MakeState (&state);

    nxf = nyf = nzf = 1;
    D_EXPAND(nxf = NMAX_POINT;  ,
             nyf = NMAX_POINT;  ,
             nzf = NMAX_POINT;)

    u      = ARRAY_2D(NMAX_POINT, NVAR, double);
    d.Vc   = ARRAY_4D(NVAR, nzf, nyf, nxf, double);
    UU_1   = ARRAY_4D(nzf, nyf, nxf, NVAR, double);
    d.flag = ARRAY_3D(nzf, nyf, nxf, unsigned char);
    
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

    C_dt[RHO] = ARRAY_3D(nzf, nyf, nxf, double);
    #if VISCOSITY == EXPLICIT
     C_dt[MX1] = ARRAY_3D(nzf, nyf, nxf, double);
    #endif
    #if RESISTIVE_MHD == EXPLICIT
     EXPAND(C_dt[BX1] = ARRAY_3D(nzf, nyf, nxf, double); ,
            C_dt[BX2] = ARRAY_3D(nzf, nyf, nxf, double); ,
            C_dt[BX3] = ARRAY_3D(nzf, nyf, nxf, double);)
    #endif
    #if THERMAL_CONDUCTION == EXPLICIT
     C_dt[ENG] = ARRAY_3D(nzf, nyf, nxf, double);
    #endif
  }

  g_intStage = 1;
  FlagReset (&d);
  getPrimitiveVars (UU, &d, grid);
  #if SHOCK_FLATTENING == MULTID
   FindShock (&d, grid);
  #endif

  #ifdef SKIP_SPLIT_CELLS
   DOM_LOOP(kk,jj,ii){
     if (splitcells[kk][jj][ii] < 0.5){
       d.flag[kk][jj][ii] |= FLAG_SPLIT_CELL;
     }
   }
  #endif

/* ----------------------------------------------------
     STEP 1: Predictor
   ---------------------------------------------------- */

  for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

    SetIndexes (&indx, grid);

    D_EXPAND(indx.beg    -= 2; indx.end    += 2; ,
             indx.t1_beg -= 2; indx.t1_end += 2;  ,
             indx.t2_beg -= 2; indx.t2_end += 2;)
             
    ResetState (&d, &state, grid);
    TRANSVERSE_LOOP(indx,in,i,j,k){  

      inv_dl = GetInverse_dl(grid);

      for (in = 0; in < indx.ntot; in++) {
      for (nv = 0; nv < NVAR; nv++) {
        state.v[in][nv] = d.Vc[nv][*k][*j][*i];
      }}

      CheckNaN (state.v, indx.beg-1, indx.end+1, 0);
      PrimToCons (state.v, u, 0, indx.ntot-1);

      for (in = 0; in < indx.ntot; in++){
      for (nv = NVAR; nv--;    ){
        state.vp[in][nv] = state.vm[in][nv] = state.v[in][nv];
      }}
      PrimToCons (state.vp, state.up, 0, indx.ntot-1);
      PrimToCons (state.vm, state.um, 0, indx.ntot-1);
      Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
      #if (PARABOLIC_FLUX & EXPLICIT)
       ParabolicFlux(d->Vc, &state, dcoeff, indx.beg - 1, indx.end, grid);
      #endif
      RightHandSide (&state, Dts, indx.beg, indx.end, 0.5*g_dt, grid);

      if (g_dir == IDIR){
        for (in = indx.beg; in <= indx.end; in++){
          #if !GET_MAX_DT
           C_dt[0][*k][*j][*i] = 0.5*(Dts->cmax[in-1] + Dts->cmax[in])*inv_dl[in];
          #endif
          #if VISCOSITY == EXPLICIT
           dl2 = 0.5*inv_dl[in]*inv_dl[in];
           C_dt[MX1][*k][*j][*i] = (dcoeff[in][MX1]+dcoeff[in-1][MX1])*dl2;
          #endif
          #if RESISTIVE_MHD == EXPLICIT
           dl2 = 0.5*inv_dl[in]*inv_dl[in];
           EXPAND(C_dt[BX1][*k][*j][*i] = (dcoeff[in-1][BX1]+dcoeff[in][BX1])*dl2;  ,
                  C_dt[BX2][*k][*j][*i] = (dcoeff[in-1][BX2]+dcoeff[in][BX2])*dl2;  ,
                  C_dt[BX3][*k][*j][*i] = (dcoeff[in-1][BX3]+dcoeff[in][BX3])*dl2;)
          #endif
          #if THERMAL_CONDUCTION  == EXPLICIT
           dl2 = 0.5*inv_dl[in]*inv_dl[in];
           C_dt[ENG][*k][*j][*i] = (dcoeff[in-1][ENG] + dcoeff[in][ENG])*dl2;
          #endif
            for (nv = NVAR; nv--;  ) UU_1[*k][*j][*i][nv] = u[in][nv] + state.rhs[in][nv];
        }
      }else{
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
          for (nv = NVAR; nv--;  ) UU_1[*k][*j][*i][nv] += state.rhs[in][nv];
        }
      }

    }
  }

/* ---- convert predictor state to primitive ---- */

  g_dir = IDIR;
  SetIndexes (&indx, grid);
  for (kk = KBEG - 2*KOFFSET; kk <= KEND + 2*KOFFSET; kk++){ g_k = &kk;
  for (jj = JBEG - 2*JOFFSET; jj <= JEND + 2*JOFFSET; jj++){ g_j = &jj;
    ConsToPrim (UU_1[kk][jj], state.v, IBEG-2, IEND+2, state.flag);
    for (ii = IBEG-2; ii <= IEND+2; ii++){
      for (nv = 0; nv < NVAR; nv++) d.Vc[nv][kk][jj][ii] = state.v[ii][nv];
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
  #ifdef GLM_MHD
   ch_max_loc = MAX(ch_max_loc, Dts->inv_dta*m_currentDlMin);
  #endif

/* ----------------------------------------------------
    STEP 2: Corrector
   ---------------------------------------------------- */

  int numFlux = numFluxes();
  a_F.resize(UBox,numFlux);
  a_F.setVal(0.0);

  g_intStage = 2;
  for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

    SetIndexes (&indx, grid);

    nxf = grid[IDIR].np_int + (g_dir == IDIR);
    nyf = grid[JDIR].np_int + (g_dir == JDIR);
    nzf = grid[KDIR].np_int + (g_dir == KDIR);

    nxb = grid[IDIR].lbeg - (g_dir == IDIR);
    nyb = grid[JDIR].lbeg - (g_dir == JDIR);
    nzb = grid[KDIR].lbeg - (g_dir == KDIR);

    ResetState (&d, &state, grid);
    TRANSVERSE_LOOP(indx,in,i,j,k){  

      for (in = 0; in < indx.ntot; in++) {
        for (nv = NVAR; nv--;  ) state.v[in][nv] = d.Vc[nv][*k][*j][*i];
      }
      States   (&state, indx.beg - 1, indx.end + 1, grid);
      Riemann  (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
      #if (PARABOLIC_FLUX & EXPLICIT)
       ParabolicFlux (d->Vc, &state, dcoeff,  indx.beg - 1, indx.end, grid);
      #endif
      RightHandSide (&state, Dts, indx.beg, indx.end, g_dt, grid);
      saveFluxes (&state, indx.beg-1, indx.end, grid);

      for (in = indx.beg; in <= indx.end; in++) {
      for (nv = NVAR; nv--;  ) {
        UU[nv][*k][*j][*i] += state.rhs[in][nv];
      }}

// Put fluxes in the FarrayBox a_F to be passed to Chombo

      for (in = indx.beg-1; in <= indx.end; in++) {
        #if AMR_EN_SWITCH == YES
         state.flux[in][ENG] = 0.0;
        #endif
        #if ENTROPY_SWITCH == YES
         state.flux[in][ENTR] = 0.0;
        #endif
        for (nv = 0; nv < NVAR; nv++) {
          indf = nv*nzf*nyf*nxf + (*k - nzb)*nyf*nxf 
                                + (*j - nyb)*nxf 
                                + (*i - nxb);
          a_F[g_dir].dataPtr(0)[indf] = state.flux[in][nv];
        }
      }
    }
  }

  #ifdef GLM_MHD
   double dtdx = g_dt/coeff_dl_min/m_dx;
   GLM_Source (UU, dtdx, grid);
  #endif

/* ----------------------------------------------
    Source terms included via operator splitting
   ---------------------------------------------- */

  #if COOLING != NO
   convertConsToPrim(UU, d.Vc, IBEG, JBEG, KBEG, 
                          IEND, JEND, KEND, grid);
   SplitSource (&d, g_dt, Dts, grid);
   convertPrimToCons(d.Vc, UU, IBEG, JBEG, KBEG, 
                         IEND, JEND, KEND, grid);
  #endif

/* ----------------------------------------------------------
    Convert total energy before returning to Chombo. 
    This should be done only inside the computational domain 
    and physical boundaries. For ease of implementation we
    carry it out everywhere (internal boundaries will be 
    overwritten later).
   ---------------------------------------------------------- */

  #if AMR_EN_SWITCH == YES
   totEnergySwitch (UU, IBEG-IOFFSET, IEND+IOFFSET, 
                        JBEG-JOFFSET, JEND+JOFFSET, 
                        KBEG-KOFFSET, KEND+KOFFSET, -1);
   #if ENTROPY_SWITCH == YES
    totEnergySwitch (UU, IBEG-IOFFSET, IEND+IOFFSET,
                         JBEG-JOFFSET, JEND+JOFFSET,
                         KBEG-KOFFSET, KEND+KOFFSET, 0);
   #endif
  #endif

/* ---------------------------------------------------------------
    In cylindrical geom. we pass U*r back to Chombo rather
    than U.
   --------------------------------------------------------------- */

  #if GEOMETRY == CYLINDRICAL
   for (nv = NVAR; nv--;   ){
   for (kk = KBEG-KOFFSET; kk <= KEND+KOFFSET; kk++){
   for (jj = JBEG-JOFFSET; jj <= JEND+JOFFSET; jj++){
   for (ii = IBEG-IOFFSET; ii <= IEND+IOFFSET; ii++){
     UU[nv][kk][jj][ii] *= grid[IDIR].x[ii];
   }}}}
  #endif
  #if GEOMETRY == SPHERICAL
   for (nv = NVAR; nv--;   ){
   for (kk = KBEG-KOFFSET; kk <= KEND+KOFFSET; kk++){
   for (jj = JBEG-JOFFSET; jj <= JEND+JOFFSET; jj++){
   for (ii = IBEG-IOFFSET; ii <= IEND+IOFFSET; ii++){
     UU[nv][kk][jj][ii] *= grid[IDIR].dV[ii]*grid[JDIR].dV[jj]/m_dx/m_dx;
   }}}}
  #endif

/* -------------------------------------------------
               Free memory 
   ------------------------------------------------- */

  for (nv = 0; nv < NVAR; nv++) FreeArrayMap(UU[nv]);

  #ifdef SKIP_SPLIT_CELLS
   FreeArrayBoxMap (splitcells, KBEG, KEND, JBEG, JEND, IBEG, IEND);
  #endif

}
