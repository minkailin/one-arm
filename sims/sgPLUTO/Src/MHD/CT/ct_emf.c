/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Store or retrieve the Electromotive Force (EMF).              

  This file provides a database functionality for storing or 
  retrieving EMF components and related information at different 
  points and times in the code.

  The CT_StoreEMF() function is called immediately after a 1D Riemann 
  solver during the hydro sweeps in order to save Fluxes and 
  characteristic signal velocities into the emf structure for later 
  reuse. 
  The fluxes coming from different sweeps are the different components
  of the advective part (-v X B) part of the electric field.

  This CT_StoreResistiveEMF() function is used to save the 
  electric field components associated with resistive terms only.

  The function CT_GetEMF() is used to obtain the edge-centered electric 
  field by properly averaging the EMF components previously stored
  at the zone faces during the 1D sweeps.

  \author  A. Mignone (mignone@ph.unito.it)
  \date    Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

static EMF emf, emf_res;

#define eps_UCT_CONTACT   1.e-6
#define EX(k,j,i)  (vz[k][j][i]*By[k][j][i] - vy[k][j][i]*Bz[k][j][i])
#define EY(k,j,i)  (vx[k][j][i]*Bz[k][j][i] - vz[k][j][i]*Bx[k][j][i])
#define EZ(k,j,i)  (vy[k][j][i]*Bx[k][j][i] - vx[k][j][i]*By[k][j][i])

/* ********************************************************************* */
void CT_StoreEMF (const State_1D *state, int beg, int end, Grid *grid)
/*!
 * Store EMF components and related information available 
 * during 1D sweeps.
 *
 * \param [in]      state pointer to State_1D structure
 * \param [in]      beg    initial index of computation 
 * \param [in]      end    final   index of computation
 * \param [in]      grid  pointer to Grid structure;
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int i, j, k, s;

/* ----------------------------------------------------
     Allocate memory for EMF structure and 
     check for incompatible combinations of algorithms 
   ---------------------------------------------------- */

  if (emf.ez == NULL){

    emf.ibeg = emf.iend = 0;
    emf.jbeg = emf.jend = 0;
    emf.kbeg = emf.kend = 0;

  /* -- memory allocation -- */

    D_EXPAND(                                          ;  ,
      emf.ez = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
      emf.ex = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
      emf.ey = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    )

    #if CT_EMF_AVERAGE == UCT_CONTACT
     D_EXPAND(
       emf.svx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, signed char);  ,
       emf.svy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, signed char);  ,
       emf.svz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, signed char);
     )
    #endif

     D_EXPAND(                                           ;  ,
       emf.ezi = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double); 
       emf.ezj = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

       emf.exj = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.exk = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

       emf.eyi = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.eyk = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     )

    #if CT_EMF_AVERAGE == UCT_HLL

     D_EXPAND(
       emf.SxL = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.SxR = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

       emf.SyL = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.SyR = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

       emf.SzL = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
       emf.SzR = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     )
    #endif
  }

/* ------------------------------------------------------
     Store emf components or other necessary 1-D data
   ------------------------------------------------------ */

  if (g_dir == IDIR){

    emf.ibeg = beg; emf.iend = end;
    for (i = beg; i <= end; i++) { 

      D_EXPAND(emf.ezi[*g_k][*g_j][i] = -state->flux[i][BX2];  ,
                                                              ,
               emf.eyi[*g_k][*g_j][i] =  state->flux[i][BX3]; ) 

      #if CT_EMF_AVERAGE == UCT_CONTACT
       if      (state->flux[i][RHO] >  eps_UCT_CONTACT) s = 1;
       else if (state->flux[i][RHO] < -eps_UCT_CONTACT) s = -1;
       else s = 0;

       emf.svx[*g_k][*g_j][i] = s;
      #endif 

      #if CT_EMF_AVERAGE == UCT_HLL
       emf.SxL[*g_k][*g_j][i] = MAX(0.0, -state->SL[i]); 
       emf.SxR[*g_k][*g_j][i] = MAX(0.0,  state->SR[i]); 
      #endif

      #if CT_EMF_AVERAGE == RIEMANN_2D
       emf.ezi[*g_k][*g_j][i] = -2.0*(state->pnt_flx[i][BX2] + 
                                          state->dff_flx[i][BX2]); 
      #endif

    }

  }else if (g_dir == JDIR){

    emf.jbeg = beg; emf.jend = end;
    for (j = beg; j <= end; j++) {

      D_EXPAND(                                                ;   ,
               emf.ezj[*g_k][j][*g_i] =  state->flux[j][BX1];   ,
               emf.exj[*g_k][j][*g_i] = -state->flux[j][BX3]; )

      #if CT_EMF_AVERAGE == UCT_CONTACT
       if      (state->flux[j][RHO] >  eps_UCT_CONTACT) s = 1;
       else if (state->flux[j][RHO] < -eps_UCT_CONTACT) s = -1;
       else s = 0;
       emf.svy[*g_k][j][*g_i] = s;
      #endif

      #if CT_EMF_AVERAGE == UCT_HLL            
       emf.SyL[*g_k][j][*g_i] = MAX(0.0, -state->SL[j]); 
       emf.SyR[*g_k][j][*g_i] = MAX(0.0,  state->SR[j]); 
      #endif

      #if CT_EMF_AVERAGE == RIEMANN_2D 
//       emf.ezj[*g_k][j][*g_i] += state->dff_flx[j][BX1];  
       emf.ezj[*g_k][j][*g_i]   = 2.0*(state->pnt_flx[j][BX1] + 
                                           state->dff_flx[j][BX1]); 
      #endif
    }

  }else if (g_dir == KDIR){

    emf.kbeg = beg; emf.kend = end;
    for (k = beg; k <= end; k++) {

      emf.eyk[k][*g_j][*g_i] = -state->flux[k][BX1]; 
      emf.exk[k][*g_j][*g_i] =  state->flux[k][BX2]; 

      #if CT_EMF_AVERAGE == UCT_CONTACT
       if      (state->flux[k][RHO] >  eps_UCT_CONTACT) s = 1;
       else if (state->flux[k][RHO] < -eps_UCT_CONTACT) s = -1;
       else s = 0;
       emf.svz[k][*g_j][*g_i] = s;
      #endif

      #if CT_EMF_AVERAGE == UCT_HLL            
       emf.SzL[k][*g_j][*g_i] = MAX(0.0, -state->SL[k]); 
       emf.SzR[k][*g_j][*g_i] = MAX(0.0,  state->SR[k]); 
      #endif

    }
  }

/* ------------------------------------------------------
         Store velocity slopes if necessary 
   ------------------------------------------------------ */

  #if CT_EMF_AVERAGE == UCT_HLL
   #ifdef CTU 
    if (g_intStage == 2) return;    
   #endif

   /* -- "end+1" needed to save dvx_dx -- */

   CT_StoreVelSlopes (&emf, state, beg, end + 1); 

  #endif

}
/* ********************************************************************* */
void CT_StoreResistiveEMF (double **res_flx, int beg, int end, Grid *grid)
/*!
 * Store resistive part of EMF.
 *
 * \todo some memory could be saved here...
 *********************************************************************** */
{
  int i, j, k;
 
  if (emf_res.ezi == NULL){
    D_EXPAND(                                               ;  ,

      emf_res.ezi = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  
      emf_res.ezj = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,

      emf_res.exj = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
      emf_res.exk = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

      emf_res.eyi = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
      emf_res.eyk = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    )
  }

/* ------------------------------------------------------
     Store emf component or other necessary 1-D data
   ------------------------------------------------------ */

  if (g_dir == IDIR){

    for (i = beg; i <= end; i++) {
      D_EXPAND(emf_res.ezi[*g_k][*g_j][i] = -res_flx[i][BX2];  ,
                                                               ;  ,
               emf_res.eyi[*g_k][*g_j][i] =  res_flx[i][BX3]; ) 
     }

  }else if (g_dir == JDIR){

    for (j = beg; j <= end; j++) {
       D_EXPAND(                                                ;  ,
                emf_res.ezj[*g_k][j][*g_i] =  res_flx[j][BX1];  ,
                emf_res.exj[*g_k][j][*g_i] = -res_flx[j][BX3]; )
    }

  }else if (g_dir == KDIR){

    for (k = beg; k <= end; k++) {
       emf_res.eyk[k][*g_j][*g_i] = -res_flx[k][BX1]; 
       emf_res.exk[k][*g_j][*g_i] =  res_flx[k][BX2]; 
    }
  }
}

/* ********************************************************************* */
EMF *CT_GetEMF (const Data *d, Grid *grid)
/*!
 * Retrieve EMF by suitable average of 1D face-centered fluxes.
 * 
 * \param [in] d
 * \param [in] grid
 *
 * \return a pointer to an edge-centered EMF.
 *********************************************************************** */
{
  int    i, j, k;
  double ***vx, ***vy, ***vz;
  double ***Bx, ***By, ***Bz;

  #if !(CT_EMF_AVERAGE == ARITHMETIC  || \
        CT_EMF_AVERAGE == RIEMANN_2D  || \
        CT_EMF_AVERAGE == UCT_CONTACT || \
        CT_EMF_AVERAGE == UCT0        || \
        CT_EMF_AVERAGE == UCT_HLL)
    print1 ("! Unknown EMF average\n");
    QUIT_PLUTO(1);
  #endif

/* -------------------------------------
       set boundary conditions on 
       face-centered electric fields
   ------------------------------------- */
/*
  #ifdef CTU
   if (g_intStage == 2)
  #endif
  EMF_BOUNDARY (&emf, grid);
*/

/* ------------------------------------------------------
       Compute slopes of staggered magnetic fields 
   ------------------------------------------------------ */

  #if CT_EMF_AVERAGE == UCT_HLL
   #ifdef CTU
    if (g_intStage == 1)
   #endif
   CT_GetStagSlopes(d->Vs, &emf, grid);
  #endif

/* -----------------------------------------------------
                 Select average 
   ----------------------------------------------------- */
/*
printf ("%d %d   %d %d\n",emf.ibeg, emf.iend, emf.jbeg, emf.jend);
printf ("%d %d   %d %d\n",0, NX1_TOT-1, 0, NX2_TOT-1);
*/
  #if CT_EMF_AVERAGE == ARITHMETIC || CT_EMF_AVERAGE == RIEMANN_2D

   CT_EMF_ArithmeticAverage (&emf, 0.25);

  #elif CT_EMF_AVERAGE == UCT_CONTACT

   CT_EMF_ArithmeticAverage (&emf, 1.0);
   CT_EMF_IntegrateToCorner (d, &emf, grid);
   for (k = emf.kbeg; k <= emf.kend; k++){
   for (j = emf.jbeg; j <= emf.jend; j++){
   for (i = emf.ibeg; i <= emf.iend; i++){      
     #if DIMENSIONS == 3
      emf.ex[k][j][i] *= 0.25;
      emf.ey[k][j][i] *= 0.25;
     #endif
     emf.ez[k][j][i] *= 0.25;
   }}}

  #elif CT_EMF_AVERAGE == UCT_HLL
   #ifdef CTU
    if (g_intStage == 1) CT_EMF_CMUSCL_Average (d, &emf, grid);
    else   
   #endif
   CT_EMF_HLL_Solver (d, &emf, grid);

  #elif CT_EMF_AVERAGE == UCT0

/* -- Subtract cell-centered contribution -- */

   EXPAND(vx = d->Vc[VX1]; Bx = d->Vc[BX1];  ,
          vy = d->Vc[VX2]; By = d->Vc[BX2];  ,
          vz = d->Vc[VX3]; Bz = d->Vc[BX3];)

   for (k = emf.kbeg; k <= emf.kend + KOFFSET; k++){
   for (j = emf.jbeg; j <= emf.jend + JOFFSET; j++){
   for (i = emf.ibeg; i <= emf.iend + IOFFSET; i++){      
     #if DIMENSIONS == 3
      emf.exj[k][j][i] *= 2.0;
      emf.exk[k][j][i] *= 2.0;
      emf.eyi[k][j][i] *= 2.0;
      emf.eyk[k][j][i] *= 2.0;

      emf.exj[k][j][i] -= 0.5*(EX(k,j,i) + EX(k,j+1,i));
      emf.exk[k][j][i] -= 0.5*(EX(k,j,i) + EX(k+1,j,i));

      emf.eyi[k][j][i] -= 0.5*(EY(k,j,i) + EY(k,j,i+1));
      emf.eyk[k][j][i] -= 0.5*(EY(k,j,i) + EY(k+1,j,i));
     #endif
     emf.ezi[k][j][i] *= 2.0;
     emf.ezj[k][j][i] *= 2.0;
     emf.ezi[k][j][i] -= 0.5*(EZ(k,j,i) + EZ(k,j,i+1));
     emf.ezj[k][j][i] -= 0.5*(EZ(k,j,i) + EZ(k,j+1,i));
   }}}

   CT_EMF_ArithmeticAverage (&emf, 0.25);

  #endif

/* ------------------------------------------------------
    Contributions from resistive terms are accounted for
    using the standard arithmetic average.
    We'll keep this in absence of something better...
   ------------------------------------------------------ */

  #if RESISTIVE_MHD == EXPLICIT 
   #ifdef CTU
//    if (g_intStage == 2) return (&emf);  
   #endif 
   for (k = emf.kbeg; k <= emf.kend; k++){
   for (j = emf.jbeg; j <= emf.jend; j++){
   for (i = emf.ibeg; i <= emf.iend; i++){      
     #if DIMENSIONS == 3
      emf.ex[k][j][i] += 0.25*( emf_res.exk[k][j][i] + emf_res.exk[k][j + 1][i] 
                              + emf_res.exj[k][j][i] + emf_res.exj[k + 1][j][i]);
      emf.ey[k][j][i] += 0.25*( emf_res.eyi[k][j][i] + emf_res.eyi[k + 1][j][i] 
                              + emf_res.eyk[k][j][i] + emf_res.eyk[k][j][i + 1]);
     #endif 
     emf.ez[k][j][i] += 0.25*( emf_res.ezi[k][j][i] + emf_res.ezi[k][j + 1][i] 
                             + emf_res.ezj[k][j][i] + emf_res.ezj[k][j][i + 1]);
   }}}
  #endif
   
/* -------------------------------------------------------------
    Fine Tuning: EMF_USERDEF_BOUNDARY can be used to directly
    set the edge-centered electric field
   ------------------------------------------------------------- */
/*
  #ifdef CTU
   if (step == 2)
  #endif
  {
    int lside[3] = {X1_BEG, X2_BEG, X3_BEG};
    int rside[3] = {X1_END, X2_END, X3_END};
    int dir;

    for (dir = 0; dir < DIMENSIONS; dir++){
      if (grid[dir].lbound == USERDEF)
        EMF_USERDEF_BOUNDARY (&emf, lside[dir], EDGE_EMF, grid);  
      if (grid[dir].rbound == USERDEF)
        EMF_USERDEF_BOUNDARY (&emf, rside[dir], EDGE_EMF, grid);  
    }   
  }
*/

  return (&emf);
}
#undef EX
#undef EY
#undef EZ

#undef eps_UCT_CONTACT
