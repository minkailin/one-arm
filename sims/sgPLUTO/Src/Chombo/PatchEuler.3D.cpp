#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

/* ********************************************************************* */
void PatchPluto::updateSolution(FArrayBox&     a_U,
                                FArrayBox&     a_Utmp,
                                FArrayBox&     split_tags,
                                BaseFab<unsigned char>& a_Flags,
                                FluxBox&       a_F,
                                Time_Step      *Dts,
                                const Box&     UBox, 
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
  double *inv_dl, dl2;
  static Data d;
  static double ***C_dt[NVAR];
 #ifdef SKIP_SPLIT_CELLS
  double ***splitcells;
 #endif
 #if (PARABOLIC_FLUX & EXPLICIT)
  static double **dcoeff;
 #endif   
  Index indx;
  static State_1D state;

  Riemann_Solver *Riemann;
  Riemann = rsolver;

#if TIME_STEPPING == RK2 
  double wflux = 0.5;
#else
  double wflux = 1.;
#endif

/* -----------------------------------------------------------------
               Check algorithm compatibilities
   ----------------------------------------------------------------- */

  #if !(GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL || GEOMETRY == SPHERICAL)
   print1 ("! RK2/EULER only works in cartesian or cylindrical/spherical coordinates\n");
   QUIT_PLUTO(1);
  #endif     

  if (NX1_TOT > NMAX_POINT || NX2_TOT > NMAX_POINT || NX3_TOT > NMAX_POINT){
    print ("!updateSolution (Euler): need to re-allocate matrix\n");
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

  #if (TIME_STEPPING == RK2)
   d.flag = ArrayCharMap(NX3_TOT, NX2_TOT, NX1_TOT,a_Flags.dataPtr(0));
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

    d.Vc   = ARRAY_4D(NVAR, nzf, nyf, nxf, double);
    #if (TIME_STEPPING != RK2)
     d.flag = ARRAY_3D(nzf, nyf, nxf, unsigned char);
    #endif 
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

  FlagReset (&d);
  getPrimitiveVars (UU, &d, grid);
  #if SHOCK_FLATTENING == MULTID
   FindShock (&d, grid);
  #endif

  #ifdef SKIP_SPLIT_CELLS
   if (g_intStage == 1) {
    DOM_LOOP(kk,jj,ii){
     if (splitcells[kk][jj][ii] < 0.5){
       d.flag[kk][jj][ii] |= FLAG_SPLIT_CELL;
     }}}
  #endif

/* ----------------------------------------------------
    Reset arrays
   ---------------------------------------------------- */

  if (g_intStage == 1) {

   #if (TIME_STEPPING == RK2)
    a_Utmp.copy(a_U); //Temporary copy of old conserved variables
   #endif

   KTOT_LOOP(kk) JTOT_LOOP(jj){

     memset ((void *)C_dt[0][kk][jj],'\0', NX1_TOT*sizeof(double));
     #if VISCOSITY == EXPLICIT
      memset ((void *)C_dt[MX1][kk][jj],'\0', NX1_TOT*sizeof(double));
     #endif
     #if RESISTIVE_MHD == EXPLICIT
      EXPAND(memset ((void *)C_dt[BX1][kk][jj],'\0', NX1_TOT*sizeof(double));  ,
             memset ((void *)C_dt[BX2][kk][jj],'\0', NX1_TOT*sizeof(double));  ,
             memset ((void *)C_dt[BX3][kk][jj],'\0', NX1_TOT*sizeof(double));)
     #endif
     #if THERMAL_CONDUCTION == EXPLICIT
      memset ((void *)C_dt[ENG][kk][jj],'\0', NX1_TOT*sizeof(double));
     #endif
   }}

/* ----------------------------------------------------
   Loop on directions
   ---------------------------------------------------- */

  int numFlux = numFluxes();
  a_F.resize(UBox,numFlux);
  a_F.setVal(0.0);

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

    inv_dl = GetInverse_dl(grid);

      for (in = 0; in < indx.ntot; in++) {
       for (nv = 0; nv < NVAR; nv++) state.v[in][nv] = d.Vc[nv][*k][*j][*i];
      }

      CheckNaN (state.v, indx.beg-1, indx.end+1, 0);

      States  (&state, indx.beg - 1, indx.end + 1, grid);
      Riemann (&state, indx.beg - 1, indx.end, Dts->cmax, grid);
      #if (PARABOLIC_FLUX & EXPLICIT)
       ParabolicFlux(d->Vc, &state, dcoeff, indx.beg - 1, indx.end, grid);
      #endif
      RightHandSide (&state, Dts, indx.beg, indx.end, g_dt, grid);
      saveFluxes (&state, indx.beg-1, indx.end, grid);

      for (in = indx.beg; in <= indx.end; in++){
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
        for (nv = NVAR; nv--;  ) UU[nv][*k][*j][*i] += state.rhs[in][nv];
      }

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
          a_F[g_dir].dataPtr(0)[indf] = wflux*state.flux[in][nv];
        }
      }

    }
  }

// Compute advective/diffusive timestep (predictor only)

  if (g_intStage == 1) {
   DOM_LOOP(kk,jj,ii){
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

   #if !GET_MAX_DT
    Dts->inv_dta /= (double)DIMENSIONS;
   #endif
   #if (PARABOLIC_FLUX & EXPLICIT)
    Dts->inv_dtp /= (double)DIMENSIONS;
   #endif
   #ifdef GLM_MHD
    ch_max_loc = MAX(ch_max_loc, Dts->inv_dta*m_currentDlMin);
   #endif
  }

// Final update: average old and new conservative variables

  #if (TIME_STEPPING == RK2)
   if (g_intStage == 2) {
    a_U.plus(a_Utmp);
    a_U *= 0.5;
  #endif

/* ----------------------------------------------
    Source terms included via operator splitting
   ---------------------------------------------- */

  #ifdef GLM_MHD
    double dtdx = g_dt/coeff_dl_min/m_dx;
    GLM_Source (UU, dtdx, grid);
  #endif

  #if INCLUDE_COOLING != NO
   convertConsToPrim(UU, d.Vc, IBEG, JBEG, KBEG,
                          IEND, JEND, KEND, grid);
   SplitSource (&d, g_dt, Dts, grid);
   convertPrimToCons(d.Vc, UU, IBEG, JBEG, KBEG,
                         IEND, JEND, KEND, grid);
  #endif

 #if (TIME_STEPPING == RK2)
  }
 #endif

/* ----------------------------------------------------------
    Convert total energy into entropy before returning to Chombo
   ---------------------------------------------------------- */

  #if AMR_EN_SWITCH == YES
   totEnergySwitch (UU, IBEG, IEND,
                        JBEG, JEND,
                        KBEG, KEND, -1);
   #if ENTROPY_SWITCH == YES
    #if (TIME_STEPPING == RK2)
     if (g_intStage == 2)
    #endif
     totEnergySwitch (UU, IBEG, IEND,
                          JBEG, JEND,
                          KBEG, KEND, 0);
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
 
  #if (TIME_STEPPING == RK2)
   FreeArrayCharMap(d.flag);
  #endif

}
