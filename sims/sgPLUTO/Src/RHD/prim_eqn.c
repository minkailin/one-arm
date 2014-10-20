#include "pluto.h"

/* *********************************************************************** */
void PrimRHS (real *w, real *dw, real cs2, real h, real *Adw)
/*
 *
 *  RHD matrix for the primitive equations
 *
 *
 * -----------------------------------------------------------
 *  Reference paper:
 *
 *   Mignone, Plewa & Bodo (2005) ApJ
 * -----------------------------------------------------------
 *
 *********************************************************************** */
{
  int nv;
  real rho, u1, u2, u3;
  real d_2, g2, scrh;
  real g, v1, v2, v3, gx2;

  #if USE_FOUR_VELOCITY == YES

   rho = w[RHO];
   EXPAND(u1  = w[VXn];  ,
          u2  = w[VXt];  ,
          u3  = w[VXb];)

   scrh = EXPAND(u1*u1, + u2*u2, + u3*u3);
   g2   = 1.0 + scrh;
   g    = sqrt(g2);
   d_2  = g/(g2 - scrh*cs2);

   /* get three-vel  */

   EXPAND(v1 = u1/g; ,
          v2 = u2/g; ,
          v3 = u3/g;)

   gx2 = 1.0/(1.0 - v1*v1);

   scrh = EXPAND(v1*dw[VXn], + v2*dw[VXt], + v3*dw[VXb]);

   Adw[PRS] =  d_2*(rho*h*cs2*(dw[VXn] - v1*scrh)
                       + u1*(1.0 - cs2)*dw[PRS]);

   Adw[RHO] = v1*dw[RHO] - (v1*dw[PRS] - Adw[PRS])/(h*cs2);
 
   scrh = 1.0/(g*rho*h);
   d_2  = u1*dw[PRS] - g*Adw[PRS];

   EXPAND(Adw[VXn] = v1*dw[VXn] + scrh*(dw[PRS] + u1*d_2);  ,
          Adw[VXt] = v1*dw[VXt] + scrh*u2*d_2;                ,
          Adw[VXb] = v1*dw[VXb] + scrh*u3*d_2;)

   for (nv = NFLX; nv < NVAR; nv++) Adw[nv] = v1*dw[nv];         

  #else

   rho = w[RHO];
   EXPAND(v1 = w[VXn];  ,
          v2 = w[VXt];  ,
          v3 = w[VXb];)

   g2  = EXPAND(v1*v1, + v2*v2, + v3*v3);
   d_2 = 1.0/(1.0 - g2*cs2);
   g2  = 1.0/(1.0 - g2);

   Adw[PRS] = d_2*(cs2*rho*h*dw[VXn] 
               + v1*(1.0 - cs2)*dw[PRS]);

   Adw[RHO] = v1*dw[RHO] - (v1*dw[PRS] - Adw[PRS])/(h*cs2);
  
   scrh = 1.0/(g2*rho*h);
   EXPAND(Adw[VXn] =  v1*dw[VXn] 
                   + scrh*(dw[PRS] - v1*Adw[PRS]); ,
          Adw[VXt] =  v1*dw[VXt] -  scrh*v2*Adw[PRS];  ,
          Adw[VXb] =  v1*dw[VXb] -  scrh*v3*Adw[PRS];)

   for (nv = NFLX; nv < NVAR; nv++) Adw[nv] = v1*dw[nv];
  #endif

}

/* ***********************************************#**************************  */
void PrimSource (const State_1D *state, int beg, int end, 
                 double *a2, double *h, double **src, Grid *grid)
/*
 *
 *
 *     Add source terms for primitive based time marching schemes;
 *     Source terms include:
 *
 *      - Geometrical sources
 *      - Gravity
 *
 *
 **************************************************************************** */
{
  int    nv, i;
  double r_1, scrh, alpha;
  double vel2, delta;
  double *q, *x1, *x2, *x3, g[3];

  #if GEOMETRY == CYLINDRICAL && (defined NEW_GEOM)
   x1 = grid[IDIR].xgc;
   x2 = grid[JDIR].xgc;
   x3 = grid[KDIR].xgc;
  #else  
   x1 = grid[IDIR].x; 
   x2 = grid[JDIR].x; 
   x3 = grid[KDIR].x; 
  #endif

  for (i = beg; i <= end; i++){
  for (nv = NVAR; nv--;  ){
    src[i][nv] = 0.0;
  }}

  #if GEOMETRY == CARTESIAN

   /* --  nothing to do here -- */
   
  #elif GEOMETRY == CYLINDRICAL 
  
   if (g_dir == IDIR){
     for (i = beg; i <= end; i++) {
 
       r_1  = 1.0/x1[i];
       q    = state->v[i];
       vel2 = EXPAND(q[VX1]*q[VX1], +q[VX2]*q[VX2], +q[VX3]*q[VX3]);

       #if USE_FOUR_VELOCITY == YES
        scrh    = sqrt(1.0 + vel2);
        alpha   = q[VXn]*r_1*scrh/(1.0 + vel2*(1.0 - a2[i]));
        scrh    = a2[i]*alpha;
        print1 ("! Primitive source terms not yet implemented\n");
        print1 ("! with 4-vel. Please try 3-vel\n");
        QUIT_PLUTO(1);
       #else
        alpha = q[VXn]*r_1/(1.0 - a2[i]*vel2);
        scrh  = a2[i]*(1.0 - vel2)*alpha;
       #endif

       src[i][RHO] = -q[RHO]*alpha;
       EXPAND (src[i][VX1] = scrh*q[VX1];  ,
               src[i][VX2] = scrh*q[VX2];  ,
               src[i][VX3] = scrh*q[VX3];)

       #if COMPONENTS == 3
        EXPAND(src[i][iVR]   +=  q[iVPHI]*q[iVPHI]*r_1;   ,
                                                          ,
               src[i][iVPHI] += -q[iVPHI]*q[iVR]*r_1;)
       #endif

       src[i][PRS] = -a2[i]*q[RHO]*h[i]*alpha;

     }
   }

  #elif GEOMETRY == SPHERICAL && DIMENSIONS == 1 && COMPONENTS == 1

   if (g_dir == IDIR){
     for (i = beg; i <= end; i++) {
 
       r_1  = 1.0/x1[i];
       q    = state->v[i];
       vel2 = EXPAND(q[VX1]*q[VX1], +q[VX2]*q[VX2], +q[VX3]*u[VX3]);

       #if USE_FOUR_VELOCITY == YES
        print1 ("! Primitive source terms not yet implemented\n");
        print1 ("! with 4-vel. Please try 3-vel\n");
        QUIT_PLUTO(1);
       #else
        delta = 1.0 - vel2*a2[i];
        scrh  = 2.0*q[iVR]/delta*r_1;
       #endif

       src[i][RHO] = -q[RHO]*scrh;
       EXPAND (src[i][VX1] = scrh*q[VX1]*a2[i]*(1.0 - vel2);  ,
               src[i][VX2] = scrh*q[VX2]*a2[i]*(1.0 - vel2);  ,
               src[i][VX3] = scrh*q[VX3]*a2[i]*(1.0 - vel2);)

       src[i][PRS] = src[i][RHO]*h[i]*a2[i];

     }
   }

  #else 

   print1 ("! Primitive source terms not available for this geometry\n");
   print1 ("! Please use RK integrators\n");
   QUIT_PLUTO(1);

  #endif   
  
/* -----------------------------------------------------------
               2 - Include Body Force
   ----------------------------------------------------------- */

  #if BODY_FORCE == EXPLICIT
  {
/*
    Force *f;
    f  = ACCELERATION (state->v, g_dir, beg, end, grid);
    for (i = beg; i <= end; i++){
      src[i][VXn] += f->a[i];
    }

    #if DIMENSIONS == 2 && COMPONENTS == 3
     if (g_dir == JDIR){
       f = ACCELERATION (state->v, KDIR, beg, end, grid);
       for (i = beg; i <= end; i++) {
         src[i][VX3] += f->a[i];
       }
     } 
    #endif
 */
    double scrh_phi, scrh_Omega[3];
    if (g_dir == IDIR) for (i = beg; i <= end; i++){
      v = state->v[i];
      BODY_FORCE (v, g, &scrh_phi, scrh_Omega, x1[i], x2[*g_j], x3[*g_k]);
      src[i][VX1] += g[IDIR];
    }
    if (g_dir == JDIR) for (i = beg; i <= end; i++){
      v = state->v[i];
      BODY_FORCE (v, g, &scrh_phi, scrh_Omega, x1[*g_i], x2[i], x3[*g_k]);
      src[i][VX2] += g[JDIR];
    }
    if (g_dir == KDIR) for (i = beg; i <= end; i++){
      v = state->v[i];
      BODY_FORCE (v, g, &scrh_phi, scrh_Omega, x1[*g_i], x2[*g_j], x3[i]);
      src[i][VX3] += g[KDIR];
    }

  }
  #endif
}


