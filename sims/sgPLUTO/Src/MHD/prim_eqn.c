/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the right hand side of the primitive MHD equations.
  
  Implements the right hand side of the quasi-linear form of the MHD 
  equations. 
  In 1D this may be written as
  \f[ 
      \partial_t{\mathbf{V}} = - A\cdot\partial_x\mathbf{V} + \mathbf{S}
  \f]
  where \f$ A \f$ is the matrix of the primitive form of the equations,
  \f$ S \f$ is the source term.

  \b Reference
 
  - "A solution adaptive upwind scheme for ideal MHD",
    Powell et al., JCP (1999) 154, 284.
 
  The function PrimRHS() implements the first term while PrimSource() 
  implements the source term part.
 
  \author A. Mignone (mignone@ph.unito.it)
  \date   Jul 31, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* *********************************************************************** */
void PrimRHS (double *v, double *dv, double cs2, double h, double *Adv)
/*!
 * Compute the matrix-vector multiplication \f$ A(\mathbf{v})\cdot 
 * d\mathbf{v} \f$ where A is the matrix of the quasi-linear form 
 * of the MHD equations.
 *
 *  \b References
 *
 *  - "A solution adaptive upwind scheme for ideal MHD",
 *    Powell et al., JCP (1999) 154, 284
 *
 *  - "An unsplit Godunov method for ideal MHD via constrained transport"
 *    Gardiner \& Stone, JCP (2005) 205, 509
 *
 * \param [in]  v    vector of primitive variables
 * \param [in]  dv   limited (linear) slopes
 * \param [in]  cs2  local sound speed
 * \param [in]  h    local enthalpy
 * \param [out] AdV  matrix-vector product
 *
 * \return This function has no return value.
 ************************************************************************* */
{
  int nv;
  real tau, scrh;
  double ch2;

  tau = 1.0/v[RHO];

  /*  ---------------------------------------------
           Adv[k]  Contains A[k][*]*dv[*]  
      ---------------------------------------------  */

  Adv[RHO] = v[VXn]*dv[RHO] + v[RHO]*dv[VXn];
  scrh = EXPAND(0.0, + v[BXt]*dv[BXt], + v[BXb]*dv[BXb]);

  #if EOS == IDEAL
   Adv[VXn] = v[VXn]*dv[VXn] + tau*(dv[PRS] + scrh);
  #elif EOS == BAROTROPIC
   Adv[VXn] = v[VXn]*dv[VXn] + tau*(cs2*dv[RHO] + scrh);
  #elif EOS == ISOTHERMAL
   Adv[VXn] = v[VXn]*dv[VXn] + tau*(cs2*dv[RHO] + scrh);
  #else 
   print ("! PRIM_RHS: not defined for this EoS\n");
   QUIT_PLUTO(1);
  #endif

  EXPAND(                                         ;    ,
         Adv[VXt] = v[VXn]*dv[VXt] - tau*v[BXn]*dv[BXt];    ,
         Adv[VXb] = v[VXn]*dv[VXb] - tau*v[BXn]*dv[BXb]; ) 

  #if MHD_FORMULATION == EIGHT_WAVES
   Adv[BXn] = v[VXn]*dv[BXn];
  #elif MHD_FORMULATION == DIV_CLEANING 
   ch2 = glm_ch*glm_ch;
   Adv[BXn]      = dv[PSI_GLM];             
   Adv[PSI_GLM] = dv[BXn]*ch2; 
  #else
   Adv[BXn] = 0.0;
  #endif
   
/* ------------------------------------------------------------------- */
/*! \note 
    In the 7-wave and 8-wave formulations we use the the same matrix 
    being decomposed into right and left eigenvectors during 
    the Characteristic Tracing step.
    Note, however, that it DOES NOT include two additional  terms 
    (-vy*dV[BX] for By, -vz*dv[BX] for Bz) that are needed in the 
    7-wave form and are added using source terms.
   --------------------------------------------------------------------- */

  EXPAND(                                                    ;   ,
         Adv[BXt] = v[BXt]*dv[VXn] - v[BXn]*dv[VXt] + v[VXn]*dv[BXt];   ,
         Adv[BXb] = v[BXb]*dv[VXn] - v[BXn]*dv[VXb] + v[VXn]*dv[BXb];)

  #if EOS == IDEAL
   Adv[PRS] = g_gamma*v[PRS]*dv[VXn] + v[VXn]*dv[PRS];
  #endif

/*  -------------------------------------------------------------
                         Now Define Tracers
    ------------------------------------------------------------- */

  #if NFLX != NVAR
   for (nv = NFLX; nv < (NFLX + NSCL); nv++) Adv[nv] = v[VXn]*dv[nv];
  #endif

}

/* ********************************************************************* */
void PrimSource (const State_1D *state, int beg, int end, double *a2, 
                  double *h, double **src, Grid *grid)
/*!
 * Compute source terms of the MHD equations in primitive variables.
 * These include:
 *
 *  - Geometrical sources;
 *  - Gravity;
 *  - terms related to divergence of B control (Powell eight wave and GLM);
 *  - FARGO source terms.
 *
 *  The rationale for choosing during which sweep a particular source 
 *  term has to be incorporated should match the same criterion used 
 *  during the conservative update. 
 *  For instance, in polar or cylindrical coordinates, curvilinear source
 *  terms are included during the radial sweep only.
 * 
 * \param [in]  state pointer to a State_1D structure
 * \param [in]  beg   initial index of computation
 * \param [in]  end   final   index of computation
 * \param [in]  a2    array of sound speed
 * \param [in]  h     array of enthalpies (not needed in MHD)
 * \param [out] src   array of source terms
 * \param [in]  grid  pointer to a Grid structure
 *
 * \note This function does not work in spherical coordinates yet. 
 *   For future implementations we annotate hereafter the induction 
 *   equation in spherical coordinates:
 *
 *  \f[ \partial_tB_r + \frac{1}{r}\partial_\theta E_\phi
 *    - \frac{1}{r\sin\theta}\partial_\phi E_\theta = -E_\phi\cot\theta/r \f]
 *  \f[ \partial_t B_\theta + \frac{1}{r\sin\theta}\partial_\phi E_r
 *    - \partial_rE_\phi =   E_\phi/r \f]
 *  \f[ \partial_t B_\phi + \partial_r E_\theta 
 *    - \frac{1}{r}\partial_\theta E_r = - E_\theta/r\f]
 *
 * where 
 *  \f[ E_\phi   = -(v \times B)_\phi   = - (v_r B_\theta - v_\theta B_r) 
 *     \,,\qquad
 *      E_\theta = -(v \times B)_\theta = - (v_\phi B_r    - v_r B_\phi) \f]
 *********************************************************************** */
{
  int nv, i;
  double tau, dA_dV;
  double hscale; /* scale factor */
  double *v, *vp,  *A, *dV, r_inv, ct;
  double *x1,  *x2,  *x3;
  double *x1p, *x2p, *x3p;
  double *dx1, *dx2, *dx3;
  static double *gPhi;
  double g[3], ch2, db, scrh;

#if ROTATING_FRAME == YES
 print1 ("! PrimSource: does not work with rotations\n");
 QUIT_PLUTO(1);
#endif

/* ----------------------------------------------------------
    memory allocation and pointer shortcuts 
   ---------------------------------------------------------- */

  if (gPhi == NULL) gPhi = ARRAY_1D(NMAX_POINT, double);

  #if GEOMETRY == CYLINDRICAL && (defined NEW_GEOM)
   x1 = grid[IDIR].xgc; x1p = grid[IDIR].xr; dx1 = grid[IDIR].dx;
   x2 = grid[JDIR].xgc; x2p = grid[JDIR].xr; dx2 = grid[JDIR].dx;
   x3 = grid[KDIR].xgc; x3p = grid[KDIR].xr; dx3 = grid[KDIR].dx;
  #else  
   x1 = grid[IDIR].x; x1p = grid[IDIR].xr; dx1 = grid[IDIR].dx;
   x2 = grid[JDIR].x; x2p = grid[JDIR].xr; dx2 = grid[JDIR].dx;
   x3 = grid[KDIR].x; x3p = grid[KDIR].xr; dx3 = grid[KDIR].dx;
  #endif
  
  #ifdef GLM_MHD
   ch2 = glm_ch*glm_ch;
  #endif

  A  = grid[g_dir].A;
  dV = grid[g_dir].dV;
  hscale  = 1.0;

/* ----------------------------------------------------------
     initialize all elements of src to zero 
   ---------------------------------------------------------- */

  for (i = beg; i <= end; i++) for (nv = NVAR; nv--;  ) src[i][nv] = 0.0;
  
/* ----------------------------------------------------------
          compute geometrical source terms
   ---------------------------------------------------------- */

#if (GEOMETRY == CARTESIAN) && (defined SHEARINGBOX)

/* -- include Coriolis terms -- */

  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) src[i][VX1] =  2.0*state->v[i][VX2]*sb_Omega;
  }else if (g_dir == JDIR){
    for (i = beg; i <= end; i++) src[i][VX2] = -2.0*state->v[i][VX1]*sb_Omega;
  }       

#elif GEOMETRY == CYLINDRICAL

  if (g_dir == IDIR) {
    for (i = beg; i <= end; i++){
      v = state->v[i]; 

      tau   = 1.0/v[RHO];
      dA_dV = 1.0/x1[i];

      src[i][RHO] = -v[RHO]*v[VXn]*dA_dV;
      EXPAND(                                                                  ,
             src[i][iBZ]   = (v[iBR]*v[iVZ] - v[iVR]*v[iBZ])*dA_dV;            ,
             src[i][iVR]   = (v[iVPHI]*v[iVPHI] - v[iBPHI]*v[iBPHI]*tau)*dA_dV; 
             src[i][iVPHI] = (-v[iVR]*v[iVPHI] + v[iBR]*v[iBPHI]*tau)*dA_dV;)
   
      #if EOS == IDEAL
       src[i][PRS] = a2[i]*src[i][RHO];
      #endif

      #ifdef GLM_MHD
       src[i][PSI_GLM] = -v[iBR]*dA_dV*ch2;
      #endif
    }
  }

#elif GEOMETRY == POLAR

  if (g_dir == IDIR) {
    for (i = beg; i <= end; i++){
      v = state->v[i]; 

      tau   = 1.0/v[RHO];
      dA_dV = 1.0/x1[i];
      src[i][RHO]  = -v[RHO]*v[VXn]*dA_dV;

      EXPAND(                                                                ,
         src[i][iVR]   = (v[iVPHI]*v[iVPHI] - v[iBPHI]*v[iBPHI]*tau)*dA_dV;  
         src[i][iVPHI] = (-v[iVR]*v[iVPHI]  + v[iBR]*v[iBPHI]*tau)*dA_dV;    ,
         src[i][iBZ]   = ( v[iBR]*v[iVZ]    - v[iVR]*v[iBZ])*dA_dV;)

      #if EOS == IDEAL
       src[i][PRS] = a2[i]*src[i][RHO];
      #endif
      #ifdef GLM_MHD
       src[i][PSI_GLM] = -v[iBR]*dA_dV*ch2;
      #endif

    }
  }

#elif GEOMETRY == SPHERICAL 

  print1 ("! PrimSource: not implemented in Spherical geometry\n");
  QUIT_PLUTO(1);
  
#endif

/* ----------------------------------------------------------
                   Add body forces 
   ---------------------------------------------------------- */

  #if (BODY_FORCE != NO)
   i = beg-1;
   if (g_dir == IDIR) {
     #if BODY_FORCE & POTENTIAL
      gPhi[i] = BodyForcePotential(x1p[i], x2[*g_j], x3[*g_k]);
     #endif
     for (i = beg; i <= end; i++){
       #if BODY_FORCE & VECTOR
        v = state->v[i];
        BodyForceVector(v, g, x1[i], x2[*g_j], x3[*g_k]);
        src[i][VX1] += g[g_dir];
       #endif
       #if BODY_FORCE & POTENTIAL
        gPhi[i]     = BodyForcePotential(x1p[i], x2[*g_j], x3[*g_k]); 
        src[i][VX1] -= (gPhi[i] - gPhi[i-1])/(hscale*dx1[i]);
       #endif

       #if DIMENSIONS == 1
        EXPAND(                        ,
               src[i][VX2] += g[JDIR];  ,
               src[i][VX3] += g[KDIR];)
       #endif
     }
   }else if (g_dir == JDIR){
     #if BODY_FORCE & POTENTIAL
      gPhi[i] = BodyForcePotential(x1[*g_i], x2p[i], x3[*g_k]);
     #endif
     for (i = beg; i <= end; i++){
       #if BODY_FORCE & VECTOR
        v = state->v[i];
        BodyForceVector(v, g, x1[*g_i], x2[i], x3[*g_k]);
        src[i][VX2] += g[g_dir];
       #endif
       #if BODY_FORCE & POTENTIAL
        gPhi[i]     = BodyForcePotential(x1[*g_i], x2p[i], x3[*g_k]);
        src[i][VX2] -= (gPhi[i] - gPhi[i-1])/(hscale*dx2[i]);
       #endif

       #if DIMENSIONS == 2 && COMPONENTS == 3
        src[i][VX3] += g[KDIR];
       #endif
     }
   }else if (g_dir == KDIR){
     #if BODY_FORCE & POTENTIAL
      gPhi[i] = BodyForcePotential(x1[*g_i], x2[*g_j], x3p[i]);
     #endif
     for (i = beg; i <= end; i++){
       #if BODY_FORCE & VECTOR
        v = state->v[i];
        BodyForceVector(v, g, x1[*g_i], x2[*g_j], x3[i]);
        src[i][VX3] += g[g_dir];
       #endif
       #if BODY_FORCE & POTENTIAL
        gPhi[i]     = BodyForcePotential(x1[*g_i], x2[*g_j], x3p[i]); 
        src[i][VX3] -=  (gPhi[i] - gPhi[i-1])/(hscale*dx3[i]);
       #endif
     }
   }
  #endif

/* -----------------------------------------------------------
      2 - MHD, div.B related source terms  
   ----------------------------------------------------------- */

  #if MHD_FORMULATION == DIV_CLEANING
   #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
    print1 ("! Error: Div. Cleaning does not work in this configuration.\n");
    print1 ("!        Try RK integrator instead\n");
    QUIT_PLUTO(1);
   #endif
   for (i = beg; i <= end; i++){
     v = state->v[i];

     tau = 1.0/v[RHO]; 
     db  = 0.5*(  A[i]  *(state->v[i+1][BXn] + state->v[i][BXn]) 
                - A[i-1]*(state->v[i-1][BXn] + state->v[i][BXn]))/dV[i];
     #if MHD_FORMULATION == DIV_CLEANING && EGLM == NO
      EXPAND(src[i][VXn] += v[BXn]*tau*db;  ,
             src[i][VXt] += v[BXt]*tau*db;  ,
             src[i][VXb] += v[BXb]*tau*db;)
     #endif
     EXPAND(                        ,
            src[i][BXt] += v[VXt]*db; ,
            src[i][BXb] += v[VXb]*db;)
     
     #if EOS == IDEAL
      scrh = EXPAND(v[VXn]*v[BXn], + v[VXt]*v[BXt], + v[VXb]*v[BXb]);
      src[i][PRS] += (1.0 - g_gamma)*scrh*db;
      #if MHD_FORMULATION == DIV_CLEANING && EGLM == NO
       scrh = 0.5*(state->v[i+1][PSI_GLM] - state->v[i-1][PSI_GLM])/grid[g_dir].dx[i];
       src[i][PRS] += (g_gamma - 1.0)*v[BXn]*scrh;
      #endif        
     #endif
   }
  #endif

/* -----------------------------------------------------------------
                        FARGO source terms

    (when SHEARINGBOX is also enabled, we do not include 
     these source terms since they're all provided by body_force)
   ----------------------------------------------------------------- */
  
  #ifdef FARGO
  #ifndef SHEARINGBOX
   #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
    print1 ("! Time Stepping works only in Cartesian or cylindrical coords\n");
    print1 ("! Use RK instead\n");
    QUIT_PLUTO(1);
   #endif

   double **wA, *dx, *dz;
   wA = FARGO_GetVelocity();
   if (g_dir == IDIR){
     dx = grid[IDIR].dx;
     for (i = beg; i <= end; i++){
       v = state->v[i];
       src[i][VX2] -= 0.5*v[VX1]*(wA[*g_k][i+1] - wA[*g_k][i-1])/dx[i];
     }
   }else if (g_dir == KDIR){
     dz = grid[KDIR].dx;
     for (i = beg; i <= end; i++){
       v = state->v[i];
       src[i][VX2] -= 0.5*v[VX3]*(wA[i+1][*g_i] - wA[i-1][*g_i])/dz[i];
     }
   }
  #endif
  #endif  
}
