/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute diffusion fluxes for explicit time stepping.

  Compute parabolic fluxes and corresponding source terms for explicit 
  time stepping only (::STS and ::RKC are treated separately) and add
  them to upwind fluxes.
  Note that source terms are only included for viscous terms, since 
  resistivity and thermal conduction do not need any.
  
  \note For explicit resistive MHD, the EMF is comprised of 2 terms:
      EMF = E(hyp) + E(par)
   The correct sequence of steps for building the EMF are:
   - E(hyp) is the flux computed with Riemann solver 
   - for STAGGERED_MHD E(hyp) is stored at appropriate 
     location for later re-use by calling ::CT_StoreEMF
   - Parabolic fluxes are stored in emf_res and stored
     in a different storage area by calling ::CT_StoreResistiveEMF
   - Parabolic fluxes are added to the hyperbolic flux 
     (useful only for cell-centered MHD).
    
   
  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    Oct 29, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

#if THERMAL_CONDUCTION != NO && (EOS == ISOTHERMAL || EOS == BAROTROPIC) 
 #error ! No Energy Equation: Thermal Conduction cannot be included
#endif

/* ********************************************************************* */
void ParabolicFlux (Data_Arr V, const State_1D *state,
                    double **dcoeff, int beg, int end, Grid *grid)
/*! 
 * Add the diffusion fluxes to the upwind fluxes for explicit time
 * integration.
 *
 * \param [in] V   pointer to the 3D array of cell-centered primitive
 *                 variables
 * \param [in,out] state pointer to a State_1D structure
 * \param [out]    dcoeff the diffusion coefficients 1D array
 * \param [in]      beg   initial index of computation
 * \param [in]      end   final   index of computation
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int i, j, k, nv;
  static double *par_Eflx;
  
/* ------------------------------------------------------------
    define par_Eflx as the sum of all parabolic contributions 
    (i.e., visc+tc+resistive) to the total energy flux.
   ------------------------------------------------------------ */

  if (par_Eflx == NULL) par_Eflx = ARRAY_1D(NMAX_POINT, double);

  for (i = beg; i <= end + 1; i++){    /* -- use memset instead -- */
    for (nv = NVAR; nv--;  ) {
      state->par_src[i][nv] = 0.0;
      state->par_flx[i][nv] = 0.0;
    }
    par_Eflx[i] = 0.0;
  }

/* -------------------------------------------------
     1. Viscosity 
   ------------------------------------------------- */
    
  #if VISCOSITY == EXPLICIT
   #ifdef FARGO
    print ("! ParabolicFlux: FARGO incompatible with explicit viscosity.\n");
    print ("                 Try STS or RKC instead\n");
    QUIT_PLUTO(1);
   #endif
   ViscousFlux (V, state->par_flx, state->par_src, dcoeff, beg, end, grid);
   for (i = beg; i <= end; i++){
     EXPAND(state->flux[i][MX1] += state->par_flx[i][MX1];  ,
            state->flux[i][MX2] += state->par_flx[i][MX2];  ,
            state->flux[i][MX3] += state->par_flx[i][MX3];  ) 
     #if EOS != ISOTHERMAL && EOS != BAROTROPIC
      par_Eflx[i]         += state->par_flx[i][ENG];
      state->flux[i][ENG] += state->par_flx[i][ENG];
     #endif

   }
  #endif 

/* -------------------------------------------------
      2. Thermal conduction
   ------------------------------------------------- */

  #if THERMAL_CONDUCTION == EXPLICIT
  {
    static int last_istep = -10;
    static double ***T;

    #ifdef CH_SPACEDIM 
     if (T == NULL) {
       int nxf = 1, nyf = 1, nzf = 1;
       D_EXPAND(nxf = NMAX_POINT;  , 
                nyf = NMAX_POINT;  ,
                nzf = NMAX_POINT;)

       T = ARRAY_3D(nzf, nyf, nxf, double); 
     }
    #else
     if (T == NULL) T = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double); 
    #endif

  /* ------------------------------------------------
      compute the temperature array just once at the 
      beginning of the next integration stage call.
     ------------------------------------------------ */

    if (g_intStage != last_istep){
      KTOT_LOOP(k) JTOT_LOOP(j) ITOT_LOOP(i){
        T[k][j][i] = V[PRS][k][j][i]/V[RHO][k][j][i];
      }
      last_istep = g_intStage;
    }

    TC_Flux (T, state, dcoeff, beg, end, grid);  
    for (i = beg; i <= end; i++) {
      par_Eflx[i]        -= state->par_flx[i][ENG];
      state->flux[i][ENG] -= state->par_flx[i][ENG];
    }
  }
  #endif 

/* -------------------------------------------------
      3. Resistivity 
   ------------------------------------------------- */

  #if RESISTIVE_MHD == EXPLICIT
   {
     ResistiveFlux (V, state->par_flx, dcoeff, beg, end, grid);
     for (i = beg; i <= end; i++){

     /* ------------------------------------------
         normal component of magnetic field does 
         not evolve during the current sweep. 
        ------------------------------------------ */

       state->par_flx[i][BXn] = 0.0;  

    /* ---------------------------------------------
        add the parabolic part of the EMF, although 
        this step is done after the hyperbolic part 
        has already been saved in CT_StoreEMF. 
        This is useful only for cell-centered MHD.
       --------------------------------------------- */

       EXPAND(state->flux[i][BX1] += state->par_flx[i][BX1];  ,
              state->flux[i][BX2] += state->par_flx[i][BX2];  ,
              state->flux[i][BX3] += state->par_flx[i][BX3]; )
       #if EOS != ISOTHERMAL && EOS != BAROTROPIC
        par_Eflx[i]        += state->par_flx[i][ENG];
        state->flux[i][ENG] += state->par_flx[i][ENG];
       #endif
     }

   /* ------------------------------------------
       now save the parabolic part of the EMF
      ------------------------------------------ */

     #ifdef STAGGERED_MHD
      CT_StoreResistiveEMF (state->par_flx, beg, end, grid);
     #endif
   }
  #endif 

/* --------------------------------------------
    at this point par_Eflx contains *ALL* 
    parabolic flux contributions. Copy its 
    content to state->par_flx for later re-use
   -------------------------------------------- */

  #if EOS != ISOTHERMAL && EOS != BAROTROPIC
   for (i = beg; i <= end; i++) state->par_flx[i][ENG] = par_Eflx[i];
  #endif

}
