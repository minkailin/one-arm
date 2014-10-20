#include "pluto.h"

#if MHD_FORMULATION == EIGHT_WAVES
/* *************************************************************************** */
void ROE_DIVB_SOURCE (const State_1D *state, int is, int ie, Grid *grid)
/* 
 *
 * PURPOSE
 *
 *   Include Powell div.B source term to momentum, induction
 *   and energy equation for Roe and TVDLF solvers.
 *
 *********************************************************************** */
{
  int    i;
  real btx, bty, btz, bx, by, bz, vx, vy, vz;
  real r, s;
  real *Ar, *Ath;
  real *vm, **bgf;
  real *src, *v;
  static real *divB, *vp;
  Grid   *GG;

  if (divB == NULL){
    divB = ARRAY_1D(NMAX_POINT, double);
    vp   = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------- ---------------------
    compute magnetic field normal component 
    interface value by arithmetic averaging
   -------------------------------------------- */

  for (i = is - 1; i <= ie; i++) {
    vp[i] = 0.5*(state->vL[i][BXn] + state->vR[i][BXn]);
  }
  vm  = vp - 1;
  GG  = grid + g_dir;
  
/* --------------------------------------------
    Compute div.B contribution from the normal 
    direction (1) in different geometries 
   -------------------------------------------- */
  
  #if GEOMETRY == CARTESIAN

   for (i = is; i <= ie; i++) {
     divB[i] = (vp[i] - vm[i])/GG->dx[i];
   }

  #elif GEOMETRY == CYLINDRICAL

   if (g_dir == IDIR){   /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i]*Ar[i] - vm[i]*Ar[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- z -- */
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i] - vm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == POLAR

   if (g_dir == IDIR){  /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i]*Ar[i] - vm[i]*Ar[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- phi -- */
     r = grid[IDIR].x[*g_i];
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i] - vm[i])/(r*GG->dx[i]);
     }
   }else if (g_dir == KDIR){  /* -- z -- */
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i] - vm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == SPHERICAL

   if (g_dir == IDIR){  /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i]*Ar[i] - vm[i]*Ar[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- theta -- */
     Ath = grid[JDIR].A;
     r   = grid[IDIR].x[*g_i];
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i]*Ath[i] - vm[i]*Ath[i - 1]) /
                 (r*GG->dV[i]);
     }
   }else if (g_dir == KDIR){  /* -- phi -- */
     r = grid[IDIR].x[*g_i];
     s = sin(grid[JDIR].x[*g_j]);
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i] - vm[i])/(r*s*GG->dx[i]);
     }
   }

  #endif

  #if BACKGROUND_FIELD == YES
   bgf = GetBackgroundField (is - 1, ie, CELL_CENTER, grid);
  #endif

/* ---------------------
     Add source terms
   -------------------- */

  for (i = is; i <= ie; i++) {

    v   = state->vh[i];
    src = state->src[i];

    EXPAND (vx = v[VX1];  ,
            vy = v[VX2];  ,
            vz = v[VX3];)

    EXPAND (bx = btx = v[BX1];  ,
            by = bty = v[BX2];  ,
            bz = btz = v[BX3];)

    #if BACKGROUND_FIELD == YES
     btx += bgf[i][BX1];
     bty += bgf[i][BX2];
     btz += bgf[i][BX3];
    #endif

    src[RHO] = 0.0;
    EXPAND (src[MX1] = -divB[i]*btx;  ,
            src[MX2] = -divB[i]*bty;  ,
            src[MX3] = -divB[i]*btz;)

    #if EOS != ISOTHERMAL && EOS != BAROTROPIC
     src[ENG] = -divB[i]*(EXPAND(vx*bx, +vy*by, +vz*bz));
    #endif
    EXPAND (src[BX1] = -divB[i]*vx;   ,
            src[BX2] = -divB[i]*vy;   ,
            src[BX3] = -divB[i]*vz;)
  }
}

/* *********************************************************************  */
void HLL_DIVB_SOURCE (const State_1D *state, real **Uhll, int beg, int end,
                      Grid *grid)
/* 
 *
 * PURPOSE
 *
 *   Include div.B source term to momentum, induction
 *   and energy equation. Used in conjunction with 
 *   an HLL-type Riemann solver. 
 * 
 * LAST_MODIFIED
 *
 *   April 7th 2006, by Andrea Mignone  (mignone@to.astro.it)
 *
 *********************************************************************** */
{
  int i, nv;
  double vc[NVAR], *A, *src, *vm;
  double r, s, vB;
  static double *divB, *vp;
  Grid *GG;

  if (divB == NULL){
    divB = ARRAY_1D(NMAX_POINT, double);
    vp   = ARRAY_1D(NMAX_POINT, double);
  }

  GG = grid + g_dir;
  vm = vp - 1;
  A  = grid[g_dir].A;

/* --------------------------------------------
    Compute normal component of the field 
   -------------------------------------------- */

  for (i = beg - 1; i <= end; i++) {
    vp[i] = state->flux[i][RHO] < 0.0 ? state->vR[i][BXn]: state->vL[i][BXn];
  }

/* --------------------------------------------
    Compute div.B contribution from the normal 
    direction (1) in different geometries 
   -------------------------------------------- */

  
  #if GEOMETRY == CARTESIAN

   for (i = beg; i <= end; i++) {
     divB[i] = (vp[i] - vm[i])/GG->dx[i];
   }

  #elif GEOMETRY == CYLINDRICAL

   if (g_dir == IDIR){   /* -- r -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- z -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i] - vm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == POLAR

   if (g_dir == IDIR){        /* -- r -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- phi -- */
     r = grid[IDIR].x[*g_i];
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i] - vm[i])/(r*GG->dx[i]);
     }
   }else if (g_dir == KDIR){  /* -- z -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i] - vm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == SPHERICAL

   if (g_dir == IDIR){  /* -- r -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- theta -- */
     r   = grid[IDIR].x[*g_i];
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/(r*GG->dV[i]);
     }
   }else if (g_dir == KDIR){  /* -- phi -- */
     r = grid[IDIR].x[*g_i];
     s = sin(grid[JDIR].x[*g_j]);
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i] - vm[i])/(r*s*GG->dx[i]);
     }
   }

  #endif

  /* -----------------------------------------
          compute total source terms
     ----------------------------------------- */
     
  for (i = beg; i <= end; i++) {

    src = state->src[i];
    for (nv = NFLX; nv--;  ) vc[nv] = state->vh[i][nv];

    vB = EXPAND(vc[VX1]*vc[BX1], + vc[VX2]*vc[BX2], + vc[VX3]*vc[BX3]);
    src[RHO] = 0.0;
    EXPAND(src[MX1] = -vc[BX1]*divB[i];  ,
           src[MX2] = -vc[BX2]*divB[i];  ,
           src[MX3] = -vc[BX3]*divB[i];)

    EXPAND(src[BX1] = -vc[VX1]*divB[i];  ,
           src[BX2] = -vc[VX2]*divB[i];  ,
	       src[BX3] = -vc[VX3]*divB[i];)

    #if EOS != ISOTHERMAL  && EOS != BAROTROPIC
     src[ENG] = -vB*divB[i];
    #endif

  }  
}
#endif

