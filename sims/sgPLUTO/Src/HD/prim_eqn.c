/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the right hand side of the primitive HD equations.     
  
  Implements the right hand side of the quasi-linear form of the hydro
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
void PrimRHS (double *w, double *dw, double cs2, double h, double *Adw)
/*!
 * Compute the matrix-vector multiplication \f$ A(\mathbf{v})\cdot 
 * d\mathbf{v} \f$ where A is the matrix of the quasi-linear form 
 * of the HD equations.
 *
 * \param [in]  w    vector of primitive variables;
 * \param [in]  dw   limited (linear) slopes;
 * \param [in]  cs2  local sound speed;
 * \param [in]  h    local enthalpy;
 * \param [out] Adw  matrix-vector product.
 * \return This function has no return value.
 ************************************************************************** */
{
  int  nv;
  real u;

  u   = w[VXn];
  
  Adw[RHO] = u*dw[RHO] + w[RHO]*dw[VXn];
  #if EOS == IDEAL
   EXPAND(Adw[VXn] = u*dw[VXn] + dw[PRS]/w[RHO];  ,
          Adw[VXt] = u*dw[VXt];                 ,
          Adw[VXb] = u*dw[VXb];)
   Adw[PRS] = cs2*w[RHO]*dw[VXn] + u*dw[PRS];
  #elif EOS == ISOTHERMAL
   EXPAND(Adw[VXn] = u*dw[VXn] + cs2*dw[RHO]/w[RHO];  ,
          Adw[VXt] = u*dw[VXt];                     ,
          Adw[VXb] = u*dw[VXb];)
  #endif

 /* ---- scalars  ---- */
     			    
  for (nv = NFLX; nv < NVAR; nv++) Adw[nv] = u*dw[nv];
}

/* ********************************************************************* */
void PrimSource (const State_1D *state, int beg, int end,
                 double *a2, double *h, double **src, Grid *grid)
/*!
 * Compute source terms of the HD equations in primitive variables.
 *
 *  - Geometrical sources;
 *  - Gravity;
 *  - Fargo source terms.
 *
 *  The rationale for choosing during which sweep a particular source 
 *  term has to be incorporated should match the same criterion used 
 *  during the conservative update. 
 *  For instance, in polar or cylindrical coordinates, curvilinear source
 *  terms are included during the radial sweep only.
 * 
 * \param [in]  state pointer to a State_1D structure;
 * \param [in]  beg   initial index of computation;
 * \param [in]  end   final   index of computation;
 * \param [in]  a2    array of sound speed; 
 * \param [in]  h     array of enthalpies (not needed in MHD);
 * \param [out] src   array of source terms;
 * \param [in]  grid  pointer to a Grid structure.
 * \return This function has no return value.
 *
 * \note This function does not work in spherical coordinates yet. 
 *********************************************************************** */
{
  int nv, i;
  double tau, dA_dV;
  double hscale; /* scale factor */
  double *v, *vp, *A, *dV, r_inv, ct;
  double *x1,  *x2,  *x3;
  double *x1p, *x2p, *x3p;
  double *dx1, *dx2, *dx3;
  static double *gPhi;
  double g[3], scrh;

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
      #if COMPONENTS == 3
       src[i][iVR]   =  v[iVPHI]*v[iVPHI]*dA_dV; 
       src[i][iVPHI] = -v[iVR]*v[iVPHI]*dA_dV;
      #endif
      #if EOS == IDEAL
       src[i][PRS] = a2[i]*src[i][RHO];
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

      #if COMPONENTS >= 2
       src[i][iVR] = v[iVPHI]*v[iVPHI]*dA_dV;  
      #endif

  /* -- in 1D, all sources should be included during this sweep -- */

      #if DIMENSIONS == 1 && COMPONENTS >= 2
       src[i][iVPHI] = -v[iVR]*v[iVPHI]*dA_dV; 
      #endif

      #if EOS == IDEAL
       src[i][PRS] = a2[i]*src[i][RHO];
      #endif
    }
  } else if (g_dir == JDIR) {
    dA_dV = 1.0/x1[*g_i];
    hscale = x1[*g_i];
    for (i = beg; i <= end; i++){
      v   = state->v[i]; 
      tau = 1.0/v[RHO];
      src[i][iVPHI] = -v[iVR]*v[iVPHI]*dA_dV;
    }
  }

#elif GEOMETRY == SPHERICAL 

  print1 ("! PrimSource: not implemented in Spherical geometry\n");
  QUIT_PLUTO(1);
  
#endif

/* ----------------------------------------------------------
                   add body forces 
   ---------------------------------------------------------- */

  #if (BODY_FORCE != NO)
   i = beg - 1;
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
