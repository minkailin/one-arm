/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the resistive MHD flux.

  Compute the resistive fluxes for the induction and energy equations.
  In the induction equation, fluxes are computed by explicitly writing
  the curl operator in components. In Cartesian components, for instance,
  one has
  \f[
     \pd{\vec{B}}{t} = - \nabla\times{\vec{E}_{\rm res}} =
      \pd{}{x}\left(\begin{array}{c}
       0           \\ \noalign{\medskip}
       \eta_z J_z  \\ \noalign{\medskip}
     - \eta_y J_y \end{array}\right)
       +
      \pd{}{y}\left(\begin{array}{c}
     - \eta_z J_z  \\ \noalign{\medskip}
           0       \\ \noalign{\medskip}
       \eta_x J_x \end{array}\right)
        +
      \pd{}{z}\left(\begin{array}{c}
        \eta_y J_y \\ \noalign{\medskip}
       -\eta_x J_x \\ \noalign{\medskip}
            0    \end{array}\right)      
       \,,\qquad
     \left(\vec{E}_{\rm res} = \tens{\eta}\cdot\vec{J}\right)
  \f]
  where \f$\tens{\eta}\f$ is the resistive diagonal tensor and 
  \f$J = \nabla\times\vec{B}\f$ is the current density.
  The corresponding contribution to the energy equation is
  \f[
    \pd{E}{t} = -\nabla\cdot\Big[(\tens{\eta}\cdot\vec{J})\times\vec{B}\Big]
              = - \pd{}{x}\left(\eta_yJ_yB_z - \eta_zJ_zB_y\right)
                - \pd{}{y}\left(\eta_zJ_zB_x - \eta_xJ_xB_z\right)
                - \pd{}{z}\left(\eta_xJ_xB_y - \eta_yJ_yB_x\right)
  \f] 
  
  \b References
     - "The PLUTO Code for Adaptive Mesh Computations in Astrophysical 
        Fluid Dynamics" \n
       Mignone et al, ApJS (2012) 198, 7M
       
  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date   Nov 5, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void ResistiveFlux (Data_Arr V, double **RF, double **dcoeff, 
                    int beg, int end, Grid *grid)
/*! 
 *  Compute resistive fluxes for the induction and energy
 *  equations. Also, compute the diffusion coefficient.
 *  Called by either ParabolicFlux or ParabolicRHS.
 *
 * \param [in]   V      3D data array of primitive variables
 * \param [out]  RF     the resistive fluxes
 * \param [out]  dcoeff the diffusion coefficients evaluated at 
 *                      cell interfaces
 * \param [in]   beg    initial index of computation
 * \param [in]   end    final   index of computation
 * \param [in]  grid    pointer to an array of Grid structures.
 *
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  double x1, x2, x3, scrh;
  double eta[3], vi[NVAR];
  static double **J;

  if (J == NULL) J = ARRAY_2D(NMAX_POINT, 3, double);
  
  GetCurrent (V, J, grid);
 
  D_EXPAND(x1 = grid[IDIR].x[*g_i];  ,
           x2 = grid[JDIR].x[*g_j];  ,
           x3 = grid[KDIR].x[*g_k]; )

/* ----------------------------------------------- 
     Compute resistive flux for the induction
     and energy equations
   ----------------------------------------------- */
  
  if (g_dir == IDIR){
  
    j = *g_j; k = *g_k;    
    for (i = beg; i <= end; i++){

    /* -- interface value -- */

      x1 = grid[IDIR].xr[i];
      for (nv = 0; nv < NVAR; nv++) {
        vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j][i + 1]);
      }
      ETA_Func (vi, x1, x2, x3, eta);

      EXPAND(                                       ,
             RF[i][BX2] = - eta[KDIR]*J[i][KDIR];   ,
             RF[i][BX3] =   eta[JDIR]*J[i][JDIR];)
      #if EOS != ISOTHERMAL
       RF[i][ENG] = EXPAND(0.0, + vi[BX2]*RF[i][BX2], + vi[BX3]*RF[i][BX3]);
      #endif
/*
      scrh = MAX(eta[0], eta[1]);
      scrh = MAX(scrh, eta[2]);
      eta_loc[i] = scrh;
*/
  /* ----------------------------------------------
      compute local inverse dt, dt^{-1} = eta/dl^2
     ---------------------------------------------- */

      EXPAND(dcoeff[i][BX1] = 0.0;        ,
             dcoeff[i][BX2] = eta[KDIR];  ,
             dcoeff[i][BX3] = eta[JDIR];)
    }

  }else if (g_dir == JDIR){
  
    i = *g_i; k = *g_k;    
    for (j = beg; j <= end; j++){

    /* -- interface value -- */

      x2 = grid[JDIR].xr[j];
      for (nv = 0; nv < NVAR; nv++) {
        vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j + 1][i]);
      }
      ETA_Func (vi, x1, x2, x3, eta);

      EXPAND(RF[j][BX1] =   eta[KDIR]*J[j][KDIR];   ,
                                                    ,
             RF[j][BX3] = - eta[IDIR]*J[j][IDIR];)
      #if EOS != ISOTHERMAL
       RF[j][ENG] = EXPAND(vi[BX1]*RF[j][BX1],     , + vi[BX3]*RF[j][BX3]);
      #endif
/*
      scrh = MAX(eta[0], eta[1]);
      scrh = MAX(scrh, eta[2]);
      eta_loc[j] = scrh;
*/
  /* ----------------------------------------------
      compute local inverse dt, dt^{-1} = eta/dl^2
     ---------------------------------------------- */

      EXPAND(dcoeff[j][BX1] = eta[KDIR];   ,
             dcoeff[j][BX2] = 0.0;         ,
             dcoeff[j][BX3] = eta[IDIR];)
    }

  }else if (g_dir == KDIR){
  
    i = *g_i; j = *g_j;    
    for (k = beg; k <= end; k++){

    /* -- interface value -- */

      x3 = grid[KDIR].xr[k];
      for (nv = 0; nv < NVAR; nv++) {
        vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k + 1][j][i]);
      }
      ETA_Func (vi, x1, x2, x3, eta);

      RF[k][BX1] = - eta[JDIR]*J[k][JDIR];
      RF[k][BX2] =   eta[IDIR]*J[k][IDIR];
      #if EOS != ISOTHERMAL
       RF[k][ENG] = vi[BX1]*RF[k][BX1] + vi[BX2]*RF[k][BX2];
      #endif
/*
      scrh = MAX(eta[0], eta[1]);
      scrh = MAX(scrh, eta[2]);
      eta_loc[k] = scrh;
*/ 
  /* ----------------------------------------------
      compute local inverse dt, dt^{-1} = eta/dl^2
     ---------------------------------------------- */

      EXPAND(dcoeff[k][BX1] = eta[JDIR];  ,
             dcoeff[k][BX2] = eta[IDIR];  ,
             dcoeff[k][BX3] = 0.0;)
    }
  }
}
