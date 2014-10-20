/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the right hand side of the conservative 
         HD/MHD equations.

  This function constructs the one-dimensional right hand side of 
  the conservative MHD or HD equations in the direction given by g_dir 
  in different geometries.
  The right hand side is computed as a two-point flux difference 
  term plus a source term:
  
  \f[ \mathrm{RHS}_i = 
      \frac{dt}{dV_i}\Big(A_{i+1/2}F_{i+1/2} - A_{i-1/2}F_{i-1/2}\Big)  
      + dt S_i  \f]
 
   where 
 
   - \f$ A_{i\pm 1/2} \f$ : interface areas
   - \f$ dV_i         \f$ : cell volume
   - \f$ F_{i\pm 1/2} \f$ : interface fluxes
   - \f$ dt           \f$ : time step
   - \f$ S_i          \f$ : source term including geometrical terms and 
                            body forces.
  
  See also \ref RHS_page
  The right hand side is assembled through the following steps:
 
   - If either one of FARGO, ROTATION or gravitational potential is 
      used, fluxes are combined to enforce conservation of total angular
      momentum and/or energy, see TotalFlux()
   - initialize rhs with flux differences
   - add geometrical source terms
   - enforce conservation of total angular  momentum and/or energy;
   - add gravity        

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void TotalFlux (const State_1D *, double *, int, int, Grid *);

#ifdef CH_SPACEDIM   /*  implies Chombo is being used  */
 #define USE_PR_GRADIENT  YES   
#else
 #ifdef FINITE_DIFFERENCE 
  #define USE_PR_GRADIENT  NO   /* -- default for Finite Difference schemes -- */
 #else
  #define USE_PR_GRADIENT  YES   /* -- default, do not change!! -- */
 #endif
#endif

#if defined FARGO && !defined SHEARINGBOX
 #define IF_FARGO(a)  a
#else 
 #define IF_FARGO(a)  
#endif

#if ROTATING_FRAME == YES
 #define IF_ROTATION(a)  a
#else 
 #define IF_ROTATION(a)  
#endif

#if EOS == IDEAL
 #define IDEAL_EOS 1
#else
 #define IDEAL_EOS 0
#endif

#if PHYSICS == MHD
#if BACKGROUND_FIELD == YES
 #define TotBB(v, b0, a, b) (v[a]*(v[b] + b0[b]) + b0[a]*v[b])
#else
 #define TotBB(v, b0, a, b) (v[a]*v[b])
#endif
#endif

/* *********************************************************************** */
void RightHandSide (const State_1D *state, Time_Step *Dts, 
              int beg, int end, double dt, Grid *grid)
/*! 
 *
 * \param [in,out]  state  pointer to State_1D structure
 * \param [in]      Dts    pointer to time step structure
 * \param [in]      beg    initial index of computation
 * \param [in]      end    final   index of computation
 * \param [in]      dt     time increment
 * \param [in]      grid  pointer to Grid structure
 *
 * \return This function has no return value.
 * \note    --
 * \todo    --
 ************************************************************************* */
{
  int  i, nv;
  double dtdx, dtdV, scrh;
  double *x1, *x1p, *dx1;
  double *x2, *x2p, *dx2;
  double *x3, *x3p, *dx3;
  double **vh, **vp, **vm;
  double **flux, **rhs, *p;
  double cl;
  double **Bg0, **wA, w, wp, vphi, gPhi_c;
  double g[3];
  static double **fA, *gPhi;

  #if GEOMETRY != CARTESIAN
   if (fA == NULL) fA = ARRAY_2D(NMAX_POINT, NVAR, double);
  #endif
  #ifdef FARGO
   wA = FARGO_GetVelocity();
  #endif

  if (gPhi == NULL) gPhi = ARRAY_1D(NMAX_POINT, double);

  #if (defined SHEARINGBOX) && (defined FARGO) && (EOS == IDEAL)
   print1 ("! ShearingBox+Fargo+Ideal EoS not properly implemented\n");
   QUIT_PLUTO(1);
  #endif

/* --------------------------------------------------
             Compute passive scalar fluxes
   -------------------------------------------------- */

  #if NSCL > 0
   AdvectFlux (state, beg - 1, end, grid);
  #endif

  #if (PHYSICS == MHD) && (BACKGROUND_FIELD == YES)
   Bg0 = GetBackgroundField (beg, end, CELL_CENTER, grid);
  #endif

/* --------------------------
      pointer shortcuts
   -------------------------- */

  rhs  = state->rhs;
  flux = state->flux;
  p    = state->press;
  vh   = state->vh;
  vp   = state->vp;
  vm   = state->vm;
  
  x1 = grid[IDIR].x; x1p = grid[IDIR].xr; dx1 = grid[IDIR].dx;
  x2 = grid[JDIR].x; x2p = grid[JDIR].xr; dx2 = grid[JDIR].dx;
  x3 = grid[KDIR].x; x3p = grid[KDIR].xr; dx3 = grid[KDIR].dx;

/* ------------------------------------------------
     Add pressure to normal component of 
     momentum flux if necessary.
   ------------------------------------------------ */

  #if USE_PR_GRADIENT == NO
   for (i = beg - 1; i <= end; i++) flux[i][MXn] += p[i];
  #endif

/* -------------------------------------------------
    compute gravitational potential and add 
    contribution to energy flux.
   ------------------------------------------------- */

  #if (defined FARGO && !defined SHEARINGBOX) ||\
      (ROTATING_FRAME == YES) || (BODY_FORCE & POTENTIAL)
   TotalFlux(state, gPhi, beg-1, end, grid);
  #endif

#if GEOMETRY == CARTESIAN
/* ***********************************************************
   
        Compute right-hand side of the MHD/HD 
        equations in CARTESIAN geometry.
    
   *********************************************************** */
{
  double x, y, z, *vc;

  if (g_dir == IDIR){

  /* ****************************************************
      Cartesian x-direction,

       - initialize rhs with flux differences (I1)
       - enforce conservation of total angular
         momentum and/or energy               (I3)
       - add gravity                          (I4)
     **************************************************** */

    y = x2[*g_j];
    z = x3[*g_k];
    for (i = beg; i <= end; i++) {
      x    = x1[i];
      dtdx = dt/dx1[i];

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      for (nv = NVAR; nv--;  ) {
        rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      }
      #if USE_PR_GRADIENT == YES
       rhs[i][MX1] -= dtdx*(p[i] - p[i-1]);
      #endif

    /* ----------------------------------------------------
       I3. modify rhs to achieve conservation
       ---------------------------------------------------- */

      #if (defined FARGO && !defined SHEARINGBOX) 
       w = wA[*g_k][i];
       rhs[i][MX2] -= w*rhs[i][RHO];
       #if IDEAL_EOS
        rhs[i][ENG] -= w*(rhs[i][MX2] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      vc = vh[i];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[*g_j], x3[*g_k]);
       rhs[i][MX1] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][MX1] -= dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if IDEAL_EOS
        gPhi_c      = BodyForcePotential(x, y, z);
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif

      #ifdef SHEARINGBOX
       rhs[i][MX1] += dt*2.0*vc[RHO]*vc[VX2]*sb_Omega;
      #endif

    }
  } else if (g_dir == JDIR){

  /* ****************************************************
      Cartesian y-direction,

       - initialize rhs with flux differences (J1)
       - add gravity                          (J4)
     **************************************************** */

    x = x1[*g_i];
    z = x3[*g_k];
    for (i = beg; i <= end; i++) {
      y    = x2[i];
      dtdx = dt/dx2[i];

    /* -----------------------------------------------
       J1. initialize rhs with flux difference
       ----------------------------------------------- */

      for (nv = NVAR; nv--;  ) {
        rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      }
      #if USE_PR_GRADIENT == YES
       rhs[i][MX2] -= dtdx*(p[i] - p[i-1]);
      #endif

    /* ----------------------------------------------------
       J4. Include gravity
       ---------------------------------------------------- */

      vc = vh[i];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[*g_i], x2[i], x3[*g_k]);
       rhs[i][MX2] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][MX2] -= dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if IDEAL_EOS
        gPhi_c      = BodyForcePotential(x, y, z);
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif

      #ifdef SHEARINGBOX
       rhs[i][MX2] -= dt*2.0*vc[RHO]*vc[VX1]*sb_Omega;
      #endif
    }

  }else if (g_dir == KDIR){

  /* ****************************************************
      Cartesian z-direction,

       - initialize rhs with flux differences (K1)
       - enforce conservation of total angular
         momentum and/or energy               (K3)
       - add gravity                          (K4)
     **************************************************** */

    x = x1[*g_i];
    y = x2[*g_j];
    for (i = beg; i <= end; i++) {
      z    = x3[i];
      dtdx = dt/dx3[i];

    /* -----------------------------------------------
       K1. initialize rhs with flux difference
       ----------------------------------------------- */

      for (nv = NVAR; nv--;  ) {
        rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      }
      #if USE_PR_GRADIENT == YES
       rhs[i][MX3] -= dtdx*(p[i] - p[i-1]);
      #endif

    /* ----------------------------------------------------
       K3. modify rhs to achieve conservation
       ---------------------------------------------------- */

      #if (defined FARGO && !defined SHEARINGBOX) 
       w = wA[i][*g_i];
       rhs[i][MX2] -= w*rhs[i][RHO];
       #if IDEAL_EOS
        rhs[i][ENG] -= w*(rhs[i][MX2] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       K4. Include gravity
       ---------------------------------------------------- */

      vc = vh[i];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[*g_i], x2[*g_j], x3[i]);
       rhs[i][MX3] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][MX3] -= dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if IDEAL_EOS
        gPhi_c      = BodyForcePotential(x, y, z);
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif
    }
  }
}
#elif GEOMETRY == CYLINDRICAL

/* ***********************************************************
   
        Compute right-hand side of the MHD/HD 
        equations in CYLINDRICAL geometry.
    
   *********************************************************** */
{
  double R, z, phi, R_1; 

  if (g_dir == IDIR) {  
    double vc[NVAR];

  /* ****************************************************
      Cylindrical radial direction:
      multiply fluxes times interface area
     **************************************************** */

    z   = x2[*g_j];
    phi = 0.0;
    for (i = beg - 1; i <= end; i++){ 
      R = grid[IDIR].A[i];

      fA[i][RHO] = flux[i][RHO]*R;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*R;     ,
             fA[i][iMZ]   = flux[i][iMZ]*R;     ,
             fA[i][iMPHI] = flux[i][iMPHI]*R*R;)       
      #if PHYSICS == MHD
       EXPAND(fA[i][iBR]   = flux[i][iBR]*R;  ,
              fA[i][iBZ]   = flux[i][iBZ]*R;  ,
              fA[i][iBPHI] = flux[i][iBPHI]*R;)
      #endif
      #if IDEAL_EOS
       fA[i][ENG] = flux[i][ENG]*R;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*R;
      #endif
      for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*R;
    }

  /* ****************************************************
      Cylindrical radial direction,

       - initialize rhs with flux differences (I1)
       - add source terms                     (I2)
       - add gravity                          (I4)
     **************************************************** */

    for (i = beg; i <= end; i++){ 
      R    = x1[i];
      dtdV = dt/grid[IDIR].dV[i];
      dtdx = dt/dx1[i];
      R_1  = 1.0/R;

    /* ---------------------------------------------------------------
       I1. initialize rhs with flux difference

           Note: when there's no explicit resistivity, we use the 
                 formulation with source terms since it seems to have 
                 better stability properties.
       ---------------------------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(rhs[i][iMR]   = - dtdV*(fA[i][iMR]   - fA[i-1][iMR]);  ,
             rhs[i][iMZ]   = - dtdV*(fA[i][iMZ]   - fA[i-1][iMZ]);  ,
             rhs[i][iMPHI] = - dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*fabs(R_1);)
      #if USE_PR_GRADIENT == YES
       rhs[i][iMR] -= dtdx*(p[i] - p[i-1]);  
      #endif
      #if PHYSICS == MHD

       #if (RESISTIVE_MHD == NO) || (RESISTIVE_MHD == SUPER_TIME_STEPPING)
        EXPAND(rhs[i][iBR]   = - dtdV*(fA[i][iBR]   - fA[i-1][iBR]);  ,
               rhs[i][iBZ]   = - dtdV*(fA[i][iBZ]   - fA[i-1][iBZ]);  ,
               rhs[i][iBPHI] = - dtdV*(fA[i][iBPHI] - fA[i-1][iBPHI]);)
       #else
        EXPAND(rhs[i][iBR]   = - dtdV*(fA[i][iBR]   - fA[i-1][iBR]);  ,
               rhs[i][iBZ]   = - dtdV*(fA[i][iBZ]   - fA[i-1][iBZ]);  ,
               rhs[i][iBPHI] = - dtdx*(flux[i][iBPHI] - flux[i-1][iBPHI]);)
       #endif

       #ifdef GLM_MHD
        rhs[i][iBR]     = - dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = - dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      #if IDEAL_EOS
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif

      for (nv = NFLX; nv < NVAR; nv++) {
        rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      }

    /* ----------------------------------------------------
       I2. Add source terms
       ---------------------------------------------------- */

      for (nv = NVAR; nv--;  ) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]); 
/* for (nv = NVAR; nv--;  ) vc[nv] = vh[i][nv]; */

      #if COMPONENTS == 3
       vphi = vc[iVPHI];
       #if ROTATING_FRAME == YES
        w     = g_OmegaZ*R;
        vphi += w;
       #endif

       #if PHYSICS == HD
        rhs[i][iMR] += dt*vc[RHO]*vphi*vphi*R_1;
       #elif PHYSICS == MHD
        rhs[i][iMR] += dt*(vc[RHO]*vphi*vphi - TotBB(vc, Bg0, iBPHI, iBPHI))*R_1;
        #if (RESISTIVE_MHD == NO) || (RESISTIVE_MHD == SUPER_TIME_STEPPING)
         rhs[i][iBPHI] -= dt*(vphi*vc[iBR] - vc[iBPHI]*vc[iVR])*R_1;
         #if BACKGROUND_FIELD == YES
          rhs[i][iBPHI] -= dt*(vphi*Bg0[i][iBR] - Bg0[i][iBPHI]*vc[iVR])*R_1;
         #endif
        #endif
       #endif /* PHYSICS == MHD */
      #endif  /* COMPONENTS == 3 */

    /* ----------------------------------------------------
       I3. modify rhs to achieve conservation
       ---------------------------------------------------- */

      #if ROTATING_FRAME == YES
       rhs[i][iMPHI] -= w*rhs[i][RHO];
       #if IDEAL_EOS
        rhs[i][ENG]  -= w*(rhs[i][iMPHI] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[*g_j], x3[*g_k]);
       rhs[i][iMR] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMR] -= dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if IDEAL_EOS
        gPhi_c      = BodyForcePotential(R, z, phi);
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif
    }
     
  } else if (g_dir == JDIR) { 
    double *vc;

  /* ****************************************************
      Cylindrical vertical direction:

       - initialize rhs with flux differences (J1)
       - add gravity                          (J4)
     **************************************************** */

    R   = x1[*g_i];
    phi = 0.0;
    for (i = beg; i <= end; i++){ 
      z    = x2[i];   
      dtdx = dt/dx2[i];

    /* -----------------------------------------------
       J1. initialize rhs with flux difference
       ----------------------------------------------- */

      for (nv = NVAR; nv--;  ) {
        rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      }
      #if USE_PR_GRADIENT == YES
       rhs[i][iMZ] += - dtdx*(p[i] - p[i-1]);
      #endif

    /* ----------------------------------------------------
       J4. Include gravity
       ---------------------------------------------------- */

      vc = vh[i];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[*g_i], x2[i], x3[*g_k]);
       rhs[i][iMZ] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMZ] += -dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if EOS == IDEAL
        gPhi_c      = BodyForcePotential(R, z, phi);
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif
    }
  }
}

#elif GEOMETRY == POLAR

/* ***********************************************************
   
        Compute right-hand side of the MHD/HD 
        equations in POLAR geometry.
    
   *********************************************************** */
{
  double R, phi, z; 
  double R_1;
   
  if (g_dir == IDIR) { 
    double vc[NVAR];

  /* ****************************************************
      Polar radial direction:
      multiply fluxes times interface area
     **************************************************** */

    phi = x2[*g_j];
    z   = x3[*g_k];
    for (i = beg - 1; i <= end; i++) { 
      R = grid[IDIR].A[i];

      fA[i][RHO] = flux[i][RHO]*R;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*R;      ,
             fA[i][iMPHI] = flux[i][iMPHI]*R*R;  ,
             fA[i][iMZ]   = flux[i][iMZ]*R;)       
      #if PHYSICS == MHD 
       EXPAND(fA[i][iBR]   = flux[i][iBR]*R;    ,
              fA[i][iBPHI] = flux[i][iBPHI];    ,
              fA[i][iBZ]   = flux[i][iBZ]*R;)
      #endif
      #if IDEAL_EOS
       fA[i][ENG] = flux[i][ENG]*R;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*R;
      #endif
      for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*R;
    }

  /* ****************************************************
      Polar radial direction,

       - initialize rhs with flux differences (I1)
       - add source terms                     (I2)
       - enforce conservation of total angular
         momentum and/or energy               (I3)
       - add gravity                          (I4)
     **************************************************** */

    for (i = beg; i <= end; i++) {
      R    = x1[i];
      dtdV = dt/grid[IDIR].dV[i];
      dtdx = dt/dx1[i];
      R_1  = grid[IDIR].r_1[i];

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(rhs[i][iMR]   = - dtdV*(fA[i][iMR]   - fA[i-1][iMR])
                             - dtdx*(p[i] - p[i-1]);                      ,      
             rhs[i][iMPHI] = - dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*R_1;  ,
             rhs[i][iMZ]   = - dtdV*(fA[i][iMZ]   - fA[i-1][iMZ]);)
      #if PHYSICS == MHD 
       EXPAND(rhs[i][iBR]   = -dtdV*(fA[i][iBR]   - fA[i-1][iBR]);    ,
              rhs[i][iBPHI] = -dtdx*(fA[i][iBPHI] - fA[i-1][iBPHI]);  ,
              rhs[i][iBZ]   = -dtdV*(fA[i][iBZ]   - fA[i-1][iBZ]);)
       #ifdef GLM_MHD
        rhs[i][iBR]     = -dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = -dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      #if IDEAL_EOS
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif
      for (nv = NFLX; nv < NVAR; nv++) {
        rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      }

    /* ----------------------------------------------------
       I2. Add source terms
       ---------------------------------------------------- */

      for (nv = NVAR; nv--;  ) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]); 
/* for (nv = NVAR; nv--;  ) vc[nv] = vh[i][nv]; */
      vphi = vc[iVPHI];
      #if (defined FARGO) || (ROTATING_FRAME == YES)
       w = 0.0;
       IF_FARGO   (w += wA[*g_k][i];)
       IF_ROTATION(w += g_OmegaZ*R;)
       vphi += w;
      #endif
      #if PHYSICS == HD
       rhs[i][iMR] += dt*vc[RHO]*vphi*vphi*R_1;
      #elif PHYSICS == MHD
       rhs[i][iMR] += dt*(  vc[RHO]*vphi*vphi 
                          - TotBB(vc, Bg0[i], iBPHI, iBPHI))*R_1;
      #endif

    /* ----------------------------------------------------
       I3. modify rhs to achieve conservation
       ---------------------------------------------------- */

      #if (defined FARGO) || (ROTATING_FRAME == YES)
       rhs[i][iMPHI] -= w*rhs[i][RHO];
       #if IDEAL_EOS
        rhs[i][ENG]   -= w*(rhs[i][iMPHI] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif
      
    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[*g_j], x3[*g_k]);
       rhs[i][iMR] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMR] -= dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if IDEAL_EOS
        gPhi_c      = BodyForcePotential(R, phi, z);
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif
    }
     
  } else if (g_dir == JDIR) {
    double *vc;

  /* ****************************************************
      Polar azimuthal direction:

       - initialize rhs with flux differences (J1)
       - add gravity                          (J4)
     **************************************************** */

    R = x1[*g_i];
    z = x3[*g_k];
    scrh = dt/R;
    for (i = beg; i <= end; i++){ 
      phi  = x2[i];
      dtdx = scrh/dx2[i];

    /* ------------------------------------------------
       J1. Compute equations rhs for phi-contributions
       ------------------------------------------------ */

      for (nv = NVAR; nv--;  ) {
        rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      }
      rhs[i][iMPHI] -= dtdx*(p[i] - p[i-1]);

    /* -------------------------------------------------------
       J4. Include gravity
       ------------------------------------------------------- */

      vc = vh[i];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[*g_i], x2[i], x3[*g_k]);
       rhs[i][iMPHI] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMPHI] -= dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if IDEAL_EOS
        gPhi_c      = BodyForcePotential(R, phi, z);
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif
    }

  } else if (g_dir == KDIR) { 
    double *vc;

  /* ****************************************************
      Polar vertical direction:

       - initialize rhs with flux differences (K1)
       - enforce conservation of total angular
         momentum and/or energy               (K3)
       - add gravity                          (K4)
     **************************************************** */

    R   = x1[*g_i];
    phi = x2[*g_j];
    for (i = beg; i <= end; i++){ 
      z    = x3[i];
      dtdx = dt/dx3[i];

    /* -----------------------------------------------
       K1. initialize rhs with flux difference
       ----------------------------------------------- */

      for (nv = NVAR; nv--;  ) {
        rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      }
      rhs[i][iMZ] -= dtdx*(p[i] - p[i-1]);

    /* ------------------------------------------------------
       K3. modify rhs to enforce conservation (FARGO only)
           (solid body rotations are not included the
            velocity depends on the cylindrical radius only)
       ------------------------------------------------------ */

      #ifdef FARGO 
       w = wA[i][*g_i];
       rhs[i][iMPHI] -= w*rhs[i][RHO];
       #if IDEAL_EOS
        rhs[i][ENG]  -= w*(rhs[i][iMPHI] + 0.5*w*rhs[i][RHO]);       
       #endif
      #endif

    /* ----------------------------------------------------
       K4. Include gravity
       ---------------------------------------------------- */

      vc = vh[i];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[*g_i], x2[*g_j], x3[i]); 
       rhs[i][iMZ] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMZ] += -dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if EOS == IDEAL
        gPhi_c      = BodyForcePotential(R, phi, z);
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif                                     
    }
  }
}
#elif GEOMETRY == SPHERICAL

/* ***********************************************************
   
        Compute right-hand side of the MHD/HD 
        equations in SPHERICAL geometry.
    
   *********************************************************** */
{
  double r, th, phi;
  double r2, r3, r_1;
  double s, s2, ct, s_1;

  if (g_dir == IDIR) { 
    double Sm, vc[NVAR];

  /* ****************************************************
      Spherical radial direction: 
      multiply fluxes by interface area 
     **************************************************** */

    th  = x2[*g_j]; s = sin(th);
    phi = x3[*g_k];
    for (i = beg - 1; i <= end; i++){
      r  = x1p[i];
      r2 = r*r; 
      r3 = r2*r;

      fA[i][RHO] = flux[i][RHO]*r2;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*r2;   ,
             fA[i][iMTH]  = flux[i][iMTH]*r2;  ,
             fA[i][iMPHI] = flux[i][iMPHI]*r3;)
      #if PHYSICS == MHD
       EXPAND(fA[i][iBR]   = flux[i][iBR]*r2;   ,
              fA[i][iBTH]  = flux[i][iBTH]*r;  ,
              fA[i][iBPHI] = flux[i][iBPHI]*r;)
      #endif
      #if IDEAL_EOS
       fA[i][ENG] = flux[i][ENG]*r2;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*r2;
      #endif
      for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*r2;
    } 

  /* ****************************************************
      Spherical radial direction:

       - initialize rhs with flux differences (I1)
       - add source terms                     (I2)
       - enforce conservation of total angular 
         momentum and/or energy               (I3)
       - add gravity                          (I4)
     **************************************************** */

    for (i = beg; i <= end; i++) { 
      r    = x1[i];
      dtdV = dt/grid[IDIR].dV[i];
      dtdx = dt/dx1[i];
      r_1  = grid[IDIR].r_1[i];

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(
        rhs[i][iMR]   = - dtdV*(fA[i][iMR] - fA[i-1][iMR])
                        - dtdx*(p[i] - p[i-1]);                    ,
        rhs[i][iMTH]  = -dtdV*(fA[i][iMTH]  - fA[i-1][iMTH]);      ,
        rhs[i][iMPHI] = -dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*r_1; 
      )
      #if PHYSICS == MHD
       EXPAND(                                                     
         rhs[i][iBR]   = -dtdV*(fA[i][iBR]   - fA[i-1][iBR]);       ,
         rhs[i][iBTH]  = -dtdx*(fA[i][iBTH]  - fA[i-1][iBTH])*r_1;  ,
         rhs[i][iBPHI] = -dtdx*(fA[i][iBPHI] - fA[i-1][iBPHI])*r_1;
       )
       #ifdef GLM_MHD
        rhs[i][iBR]     = -dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = -dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      #if IDEAL_EOS
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif

      for (nv = NFLX; nv < NVAR; nv++){
        rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      }

    /* ----------------------------------------------------
       I2. Add source terms 
       ---------------------------------------------------- */
  
      for (nv = NVAR; nv--;  ) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]);
/*  for (nv = NVAR; nv--;  ) vc[nv] = vh[i][nv]; */

      vphi = SELECT(0.0, 0.0, vc[iVPHI]);
      #if (defined FARGO) || (ROTATING_FRAME == YES)
       w = 0.0;
       IF_FARGO   (w += wA[*g_j][i];)
       IF_ROTATION(w += g_OmegaZ*r*s;)
       vphi += w;
      #endif
      Sm = vc[RHO]*(EXPAND(  0.0, + vc[iVTH]*vc[iVTH], + vphi*vphi));
      #if PHYSICS == MHD
       Sm += EXPAND(  0.0, - TotBB(vc, Bg0[i], iBTH, iBTH), 
                           - TotBB(vc, Bg0[i], iBPHI,iBPHI));
      #endif
      rhs[i][iMR] += dt*Sm*r_1;

    /* ----------------------------------------------------
       I3. modify rhs to enforce conservation
       ---------------------------------------------------- */

      #if (defined FARGO) || (ROTATING_FRAME == YES)
       rhs[i][iMPHI] -= w*rhs[i][RHO];
       #if IDEAL_EOS
        rhs[i][ENG]   -= w*(rhs[i][iMPHI] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[*g_j], x3[*g_k]); 
       rhs[i][iMR] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir]; 
       #endif
      #endif
      
      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMR] -= dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if IDEAL_EOS
        gPhi_c      = BodyForcePotential(r, th, phi); 
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif                                     
    }

  } else if (g_dir == JDIR) {

    double Sm, *vc;

  /* ****************************************************
      Spherical meridional direction:
      multiply fluxes by zone-interface area
     **************************************************** */

    r   = x1[*g_i];
    phi = x3[*g_k];
    for (i = beg - 1; i <= end; i++){ 
      s  = grid[JDIR].A[i];
      s2 = s*s;

      fA[i][RHO] = flux[i][RHO]*s;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*s;   ,
             fA[i][iMTH]  = flux[i][iMTH]*s;  ,
             fA[i][iMPHI] = flux[i][iMPHI]*s2;) 
      #if PHYSICS == MHD
       EXPAND(fA[i][iBR]   = flux[i][iBR]*s;   ,
              fA[i][iBTH]  = flux[i][iBTH]*s;  ,
              fA[i][iBPHI] = flux[i][iBPHI];)
      #endif
      #if IDEAL_EOS
       fA[i][ENG] = flux[i][ENG]*s;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*s;
      #endif
      for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*s;
    }

  /* ****************************************************
      Spherical meridional direction:

       - initialize rhs with flux differences (J1)
       - add source terms                     (J2)
       - enforce conservation of total angular
         momentum and/or energy               (J3)
       - add gravity                          (J4)
     **************************************************** */
    
    r_1 = grid[IDIR].r_1[*g_i];
    for (i = beg; i <= end; i++){
      th   = x2[i];
      dtdV = dt/grid[JDIR].dV[i]*r_1;
      dtdx = dt/dx2[i]*r_1;      
      s    = sin(th);
      s_1  = 1.0/s;   
      ct   = grid[JDIR].ct[i];         /* = cot(theta)  */

    /* -----------------------------------------------
       J1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(
        rhs[i][iMR]   = - dtdV*(fA[i][iMR] - fA[i-1][iMR]);  , 
        rhs[i][iMTH]  = - dtdV*(fA[i][iMTH] - fA[i-1][iMTH])
                        - dtdx*(p[i] - p[i-1]);              , 
        rhs[i][iMPHI] = - dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*fabs(s_1);
      )       
      #if PHYSICS == MHD
       EXPAND(                                                     
         rhs[i][iBR]   = -dtdV*(fA[i][iBR]   - fA[i-1][iBR]);  ,
         rhs[i][iBTH]  = -dtdV*(fA[i][iBTH]  - fA[i-1][iBTH]);  ,
         rhs[i][iBPHI] = -dtdx*(fA[i][iBPHI] - fA[i-1][iBPHI]);
       )
       #ifdef GLM_MHD
        rhs[i][iBTH]    = -dtdx*(flux[i][iBTH] - flux[i-1][iBTH]);
        rhs[i][PSI_GLM] = -dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      #if IDEAL_EOS
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif
      for (nv = NFLX; nv < NVAR; nv++){
        rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      }

    /* ----------------------------------------------------
       J2. Add source terms
       ---------------------------------------------------- */
       
      vc = vh[i];
      
      vphi = SELECT(0.0, 0.0, vc[iVPHI]);
      #if (defined FARGO) || (ROTATING_FRAME == YES)
       w = 0.0; 
       IF_FARGO   (w += wA[i][*g_i];)
       IF_ROTATION(w += g_OmegaZ*r*s;)
       vphi += w;
      #endif
      Sm = vc[RHO]*(EXPAND(  0.0, - vc[iVTH]*vc[iVR], + ct*vphi*vphi));
      #if PHYSICS == MHD
       Sm += EXPAND(0.0, +    TotBB(vc, Bg0[i], iBTH, iBR), 
                         - ct*TotBB(vc, Bg0[i], iBPHI, iBPHI));
      #endif
      rhs[i][iMTH] += dt*Sm*r_1;

    /* ----------------------------------------------------
       J3. modify rhs to enforce conservation
       ---------------------------------------------------- */

      #if (defined FARGO) || (ROTATING_FRAME == YES)
       rhs[i][iMPHI] -= w*rhs[i][RHO];
       #if IDEAL_EOS
        rhs[i][ENG]   -= w*(rhs[i][iMPHI] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       J4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[*g_i], x2[i], x3[*g_k]);
       rhs[i][iMTH] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMTH] -= dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if IDEAL_EOS
        gPhi_c      = BodyForcePotential(r, th, phi); 
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif
    }

  } else if (g_dir == KDIR) {

    double Sm, *vc;

  /* ****************************************************
      Spherical azimuthal direction:

       - initialize rhs with flux differences (K1)
       - add gravity                          (K4)
     **************************************************** */

    r  = x1[*g_i];
    th = x2[*g_j];
    s    = sin(th);
    scrh = dt/(r*s);
    for (i = beg; i <= end; i++) {
      phi  = x3[i];
      dtdx = scrh/dx3[i];

    /* ------------------------------------------------
       K1.  initialize rhs with flux difference
       ------------------------------------------------ */

      for (nv = NVAR; nv--;  ) {
        rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      }       
      rhs[i][iMPHI] -= dtdx*(p[i] - p[i-1]); 

    /* -------------------------------------------------------
       K4. Include gravity
       ------------------------------------------------------- */

      vc = vh[i];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[*g_i], x2[*g_j], x3[i]);
       rhs[i][iMPHI] += dt*vc[RHO]*g[g_dir];
       #if IDEAL_EOS
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMPHI] -= dtdx*vc[RHO]*(gPhi[i] - gPhi[i-1]);
       #if IDEAL_EOS
        gPhi_c      = BodyForcePotential(r, th, phi); 
        rhs[i][ENG] -= gPhi_c*rhs[i][RHO];
       #endif
      #endif
    }
  }
}
#endif  /* GEOMETRY == SPHERICAL */

/* ---------------------------------------------------------------
    Source terms coming from tensor discretazion of parabolic
    terms in curvilinear coordinates (only for viscosity)
  ---------------------------------------------------------------- */

  #if GEOMETRY != CARTESIAN
   #if VISCOSITY == EXPLICIT
    for (i = beg; i <= end; i++) {
      EXPAND(rhs[i][MX1] += dt*state->par_src[i][MX1];  ,
             rhs[i][MX2] += dt*state->par_src[i][MX2];  ,
             rhs[i][MX3] += dt*state->par_src[i][MX3];)
    }
   #endif
  #endif

/* --------------------------------------------------
              Powell's source terms
   -------------------------------------------------- */

  #if PHYSICS == MHD 
   #if MHD_FORMULATION == EIGHT_WAVES
    for (i = beg; i <= end; i++) {
      EXPAND(rhs[i][MX1] += dt*state->src[i][MX1];  ,
             rhs[i][MX2] += dt*state->src[i][MX2];  ,
             rhs[i][MX3] += dt*state->src[i][MX3];)

      EXPAND(rhs[i][BX1] += dt*state->src[i][BX1];  ,
             rhs[i][BX2] += dt*state->src[i][BX2];  ,
             rhs[i][BX3] += dt*state->src[i][BX3];)
      #if IDEAL_EOS
       rhs[i][ENG] += dt*state->src[i][ENG];
      #endif
    }
   #endif
  #endif

/* -------------------------------------------------
            Extended GLM source terms
   ------------------------------------------------- */

  #ifdef GLM_MHD
  #if   EGLM == YES
   EGLM_Source (state, dt, beg, end, grid);
  #endif
  #endif

/* -----------------------------------------------
      Entropy equation source terms 
   ----------------------------------------------- */

  #if RESISTIVE_MHD == EXPLICIT
   #if (ENTROPY_SWITCH == YES) && (EOS != ISOTHERMAL && EOS != BAROTROPIC)
    for (i = beg; i <= end; i++) {
      rhs[i][ENTR] += dt*state->src[i][ENTR];
    }
   #endif
  #endif

/* --------------------------------------------------
    Reset right hand side in internal boundary zones
   -------------------------------------------------- */
   
  #if INTERNAL_BOUNDARY == YES
   InternalBoundaryReset(state, Dts, beg, end, grid);
  #endif
  
/* --------------------------------------------------
           Time step determination
   -------------------------------------------------- */

#if !GET_MAX_DT
return;
#endif

  cl = 0.0;
  for (i = beg-1; i <= end; i++) {
    scrh = Dts->cmax[i]*grid[g_dir].inv_dxi[i];
    cl = MAX(cl, scrh);   
  }
  #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
   if (g_dir == JDIR) {
     cl /= fabs(grid[IDIR].xgc[*g_i]);
   }
   #if GEOMETRY == SPHERICAL
    if (g_dir == KDIR){
      cl /= fabs(grid[IDIR].xgc[*g_i])*sin(grid[JDIR].xgc[*g_j]);
    }
   #endif
  #endif
  Dts->inv_dta = MAX(cl, Dts->inv_dta);  
}

/* ********************************************************************* */
void TotalFlux (const State_1D *state, double *gPhi,
                int beg, int end, Grid *grid)
/*!
 *  Compute the total flux in order to enforce conservation of 
 *  angular momentum and energy in presence of FARGO source 
 *  terms, rotation or gravitational potential.
 *
 * \param [in]     state pointer to State_1D structure;
 * \param [in,out] gPhi  1D array defining the gravitational potential;
 * \param [in]     beg    initial index of computation; 
 * \param [in]     end    final   index of computation;
 * \param [in]     grid  pointer to Grid structure;
 *********************************************************************** */
#ifndef iMPHI
 #define iMPHI MY  /* -- for Cartesian coordinates -- */
#endif
{
  int i;
  double wp, R;
  double **flux, *vp;
  double *x1,  *x2,  *x3;
  double *x1p, *x2p, *x3p;
  #ifdef FARGO
   double **wA;
   wA = FARGO_GetVelocity();
  #endif

  flux = state->flux;
  x1  = grid[IDIR].x;  x1p = grid[IDIR].xr;
  x2  = grid[JDIR].x;  x2p = grid[JDIR].xr;
  x3  = grid[KDIR].x;  x3p = grid[KDIR].xr;

  if (g_dir == IDIR){ 
    for (i = beg; i <= end; i++){

    /* ----------------------------------------------------
        include flux contributions from FARGO or Rotation 
        Note: ShearingBox terms are not included here but
              only within the BodyForce function.
       ---------------------------------------------------- */

      #if (defined FARGO && !defined SHEARINGBOX) || (ROTATING_FRAME == YES)
       wp = 0.0;
       #if GEOMETRY == SPHERICAL
        IF_FARGO(wp = 0.5*(wA[*g_j][i] + wA[*g_j][i+1]);)
        R = x1p[i]*sin(x2[*g_j]);  /* -- cylindrical radius -- */
       #else
        IF_FARGO(wp = 0.5*(wA[*g_k][i] + wA[*g_k][i+1]);)
        R = x1p[i];                   /* -- cylindrical radius -- */
       #endif
       IF_ROTATION(wp += g_OmegaZ*R;)
       #if IDEAL_EOS
        flux[i][ENG] += wp*(0.5*wp*flux[i][RHO] + flux[i][iMPHI]);
       #endif
       flux[i][iMPHI] += wp*flux[i][RHO];
      #endif

    /* -- gravitational potential -- */

      #if (BODY_FORCE & POTENTIAL)
       gPhi[i] = BodyForcePotential(x1p[i], x2[*g_j], x3[*g_k]);
       #if IDEAL_EOS
        flux[i][ENG] += flux[i][RHO]*gPhi[i];                          
       #endif
      #endif
    }
  }else if (g_dir == JDIR){
    for (i = beg; i <= end; i++){ 

    /* ----------------------------------------------------
        include flux contributions from FARGO and Rotation 
       ---------------------------------------------------- */

      #if GEOMETRY == SPHERICAL
       #if defined FARGO || (ROTATING_FRAME == YES)
        wp = 0.0;
        R  = x1[*g_i]*sin(x2p[i]);
        IF_FARGO   (wp += 0.5*(wA[i][*g_i] + wA[i+1][*g_i]);)
        IF_ROTATION(wp += g_OmegaZ*R;)
        #if EOS != ISOTHERMAL && EOS != BAROTROPIC
         flux[i][ENG] += wp*(0.5*wp*flux[i][RHO] + flux[i][iMPHI]);
        #endif
        flux[i][iMPHI] += wp*flux[i][RHO];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       gPhi[i] = BodyForcePotential(x1[*g_i], x2p[i], x3[*g_k]);
       #if IDEAL_EOS
        flux[i][ENG] += flux[i][RHO]*gPhi[i];
       #endif
      #endif      

    }
  }else if (g_dir == KDIR){
    R = x1[*g_i];
    for (i = beg; i <= end; i++) {

    /* ----------------------------------------------------
        include flux contributions from FARGO
        (polar/carteisian geometries only)
       ---------------------------------------------------- */

      #if (GEOMETRY != SPHERICAL) && (defined FARGO) && (!defined SHEARINGBOX)
       wp = 0.5*(wA[i][*g_i] + wA[i+1][*g_i]);
       #if EOS == IDEAL
        flux[i][ENG] += wp*(0.5*wp*flux[i][RHO] + flux[i][iMPHI]);
       #endif
       flux[i][iMPHI] += wp*flux[i][RHO];
      #endif

      #if (BODY_FORCE & POTENTIAL)
       gPhi[i] = BodyForcePotential(x1[*g_i], x2[*g_j], x3p[i]); 
       #if IDEAL_EOS
        flux[i][ENG] += flux[i][RHO]*gPhi[i];                          
       #endif
      #endif
    }
  }
}
#undef TotBB
#undef IF_FARGO
#undef IF_ROTATION
#undef IDEAL_EOS
