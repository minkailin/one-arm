#include "pluto.h"

#define MAX_ITER   20
#define small_p    1.e-12
#define small_rho  1.e-12
#define INTERPOLATE_RAREFACTION   YES

#define accuracy   1.e-6

static real GLORENTZ (real *U, int n);
static void   SHOCK (real tau0, real u0, real p0, real g0, real V0, real h0,
                     real p1, real *u1, real *dudp, real *zeta, int istate);
static real RAREFACTION_SPEED (real *u, int side);

static real qglob_r[NFLX], qglob_l[NFLX], gmmr;

/* ********************************************************************  */
void TwoShock_Solver (const State_1D *state, int beg, int end, 
                      real *cmax, Grid *grid)
/*
 *
 * NAME
 *
 *   RIEMANN
 *
 *
 * PURPOSE
 *
 *   Solve riemann problem for the relativistic Euler equations 
 *   using the two-shock approximation; return numerical
 *   fluxes and source terms
 * 
 *   Reference:    "The Piecewise Parabolic Method for  
 *                  Multidimensional Relativistic Fluid Dynamics"
 *                  A. Mignone, T. Plewa and G. Bodo
 *   
 *  - On input, it takes left and right primitive state
 *     vectors state->vL and state->vR at zone edge i+1/2;
 *     On output, return flux and pressure vectors at the
 *     same interface.
 *
 *  - Also, compute maximum wave propagation speed (cmax) 
 *    for  explicit time step computation
 *  
 *
 * LAST_MODIFIED
 *
 *   April 4th 2006, by Andrea Mignone  (mignone@to.astro.it)
 *
 ******************************************************************************* */
{
  int     iter, i, nv;
  int     k, nfail = 0, izone_fail[1024];
  real  Ustar[NFLX];
  real  gL, gR, gL1, gR1;
  real  uxR, uxR1, pR, j2R, wR, VR, VR1;
  real  uxL, uxL1, pL, j2L, wL, VL, VL1;
  real  duR, duL, p1, u1, dp, a0, a1;
  real  tauR, tauL, am, ap;
  real  ql[NFLX], qr[NFLX], *qs;
  real  *vl, *vr;
  static double  **fl, **us, **ws;
  static double *hR, *hL, *a2L, *a2R, *cmax_loc, *cmin_loc;

  static double **fL, **fR, *prL, *prR;
  double bmin, bmax;
  double SL, SR;

  gmmr = g_gamma/(g_gamma - 1.0);

  if (fl == NULL){
    fl       = ARRAY_2D(NMAX_POINT, NFLX, double); 
    us       = ARRAY_2D(NMAX_POINT, NVAR, double);  
    ws       = ARRAY_2D(NMAX_POINT, NVAR, double);  
    cmax_loc = ARRAY_1D(NMAX_POINT, double);
    cmin_loc = ARRAY_1D(NMAX_POINT, double);

    a2R      = ARRAY_1D(NMAX_POINT, double);
    a2L      = ARRAY_1D(NMAX_POINT, double);
    hR       = ARRAY_1D(NMAX_POINT, double);
    hL       = ARRAY_1D(NMAX_POINT, double);

    fL  = ARRAY_2D(NMAX_POINT, NVAR, double);
    fR  = ARRAY_2D(NMAX_POINT, NVAR, double);
    prL = ARRAY_1D(NMAX_POINT, double);
    prR = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------------
        compute sound speed at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, hL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, hR, beg, end, FACE_CENTER, grid);
/*
  ENTHALPY (state->vR, hR, beg, end);
  ENTHALPY (state->vL, hL, beg, end);
*/
/*  ---------------------------------------------------------------
                       SOLVE RIEMANN PROBLEM
    ---------------------------------------------------------------   */

  for (i = beg; i <= end; i++) {

    #if SHOCK_FLATTENING == MULTID

      /* ---------------------------------------------
          Switch to HLL in proximity of strong shocks.
         --------------------------------------------- */

     if (CheckZone(i, FLAG_HLL) || CheckZone(i+1, FLAG_HLL)){

       HLL_Speed (state->vL, state->vR, a2L, a2R, &gL - i, &gR - i, i, i);
       Flux  (state->uL, state->vL, a2L, fL, prL, i, i);
       Flux  (state->uR, state->vR, a2R, fR, prR, i, i);

       a0 = MAX(fabs(gR), fabs(gL));
       cmax[i] = a0;

       gL = MIN(0.0, gL);
       gR = MAX(0.0, gR);
       a0 = 1.0/(gR - gL);
       for (nv = NFLX; nv--; ){
         state->flux[i][nv]  = gL*gR*(state->uR[i][nv] - state->uL[i][nv])
                              + gR*fL[i][nv] - gL*fR[i][nv];
         state->flux[i][nv] *= a0;
       }
       state->press[i] = (gR*prL[i] - gL*prR[i])*a0;
       continue;
     }
    #endif

    vl = state->vL[i];
    vr = state->vR[i];

    #if USE_FOUR_VELOCITY == YES  

    /*  ----  Transform four-velocity interface
              values to three-velocity           ----   */

     gL = 1.0 + EXPAND( vl[VXn]*vl[VXn],
                      + vl[VXt]*vl[VXt],
                      + vl[VXb]*vl[VXb]);

     gR = 1.0 + EXPAND( vr[VXn]*vr[VXn],
                      + vr[VXt]*vr[VXt],
                      + vr[VXb]*vr[VXb]);
     gL = sqrt(gL);
     gR = sqrt(gR);
    
     for (nv = 0; nv < NFLX; nv++) {
       qglob_l[nv] = ql[nv] = vl[nv];
       qglob_r[nv] = qr[nv] = vr[nv];
     }

   /*  ----  define left and right three-velocities  ----  */
    
     EXPAND(ql[VX1] /= gL;  ,
            ql[VX2] /= gL;  ,
            ql[VX3] /= gL;)

     EXPAND(qr[VX1] /= gR;  ,
            qr[VX2] /= gR;  ,
            qr[VX3] /= gR;)

    #else

     for (nv = 0; nv < NFLX; nv++) {
       ql[nv] = vl[nv];
       qr[nv] = vr[nv];
     }

     gR = GLORENTZ (qr,0);
     gL = GLORENTZ (ql,1);

    #endif 

    tauR = 1.0/qr[RHO];
    tauL = 1.0/ql[RHO];

    uxR = qr[VXn];
    uxL = ql[VXn];

    pR = qr[PRS];
    pL = ql[PRS];

    wR = pR*tauR;
    wL = pL*tauL;

    VR = tauR/gR;
    VL = tauL/gL;

/*  ---------------------------------------------
       First iteration is done out of the loop
    ---------------------------------------------  */

    #if EOS == IDEAL
     j2R = tauR*(hR[i]*(gmmr - 2.0) + 1.0) / (gmmr*pR);
     j2L = tauL*(hL[i]*(gmmr - 2.0) + 1.0) / (gmmr*pL);
    #endif
    #if EOS == TAUB
     a0  = (5.0*hR[i] - 8.0*wR)/(2.0*hR[i] - 5.0*wR);
     j2R = -tauR*(a0*wR + hR[i]*(1.0 - a0))/(a0*pR); 

     a0  = (5.0*hL[i] - 8.0*wL)/(2.0*hL[i] - 5.0*wL);
     j2L = -tauL*(a0*wL + hL[i]*(1.0 - a0))/(a0*pL);
    #endif

    gR1 = sqrt(VR*VR + (1.0 - uxR*uxR)*j2R);   /*    RIGHT    */
    wR  = VR*uxR + gR1;

    gL1 = -sqrt(VL*VL + (1.0 - uxL*uxL)*j2L);  /*    LEFT    */
    wL  = VL*uxL + gL1;

    duR = gR1/(hR[i]*gR);
    duL = gL1/(hL[i]*gL);

    p1  = ql[VXn] - qr[VXn] + duR*pR - duL*pL;
    p1 /= duR - duL;

    if (p1 < 0.0) {
      p1 = MIN(pR,pL);
    }

/*  -----------------------------------------------
             BEGIN iteration  loop   
    -----------------------------------------------  */

    for (iter = 1; iter < MAX_ITER; iter++)  {

      SHOCK (tauL, uxL, pL, gL, VL, hL[i], p1, &uxL1, &duL, &wL, -1);
      SHOCK (tauR, uxR, pR, gR, VR, hR[i], p1, &uxR1, &duR, &wR, 1);

  /*  ----  Find  next  approximate solution  ----  */

      dp  = (uxR1 - uxL1)/(duL - duR);
      if (-dp > p1) dp = -0.5*p1;
      p1 += dp;

      if (p1 < 0.0){
        p1 -= dp;
        p1 *= 0.5;
      }
      g_maxRiemannIter = MAX(g_maxRiemannIter,iter);
      if ( (fabs(dp) < accuracy*p1) ) break;
    }

/*  -----------------------------------------------
             END iteration  loop   
    -----------------------------------------------  */

    u1 = 0.5*(uxR1 + uxL1); 
    
  /*  ---- Check possible failures, and flag zones  ----  */

    if (u1 != u1 || iter == MAX_ITER) {

      u1 = 0.5*(uxL + uxR);
      p1 = 0.5*(pL  + pR);
      nfail++;
      izone_fail[nfail] = i;
      for (nv = NFLX; nv--; ){
        ws[i][nv] = vl[nv];
      }
      continue;
    } 

/*  This switch works to prevent shear smearing, by
    filtering noise in the velocity  */

/*
    u1 = (fabs(u1)<1.e-13 ? 0.0:u1);
*/

    Ustar[VXn] = u1;
    Ustar[PRS] = p1;

/*  ------------------------------------------------------------------ 
             Sample solution on  x/t = 0 axis       
    ------------------------------------------------------------------ */

    if (u1 >= 0.0) {          /*  ----  Solution is sampled to the LEFT of contact  ----  */

      EXPAND(       ,
             gL1   = gL*hL[i]/(gL*hL[i] + (p1 - pL)*(VL + wL*uxL));
             Ustar[VXt] = ql[VXt]*gL1;  ,
             Ustar[VXb] = ql[VXb]*gL1;)

      VL1        = VL - (u1 - uxL)*wL;
      gL1        = GLORENTZ (Ustar,2);
      Ustar[RHO] = MAX(small_rho,1.0/(VL1*gL1));

/*     Shock or rarefaction ?   */

      #if INTERPOLATE_RAREFACTION == YES
       if (p1 < ql[PRS]){
         am  = RAREFACTION_SPEED (ql,    -1);       /* -- rarefaction tail --  */
         ap  = RAREFACTION_SPEED (Ustar, -1);       /* -- rarefaction head --  */
         am  = MIN(am, ap);
       }else{        
         ap = am = VL/wL + uxL;              /*  --  SHOCK SPEED  --  */
       }
      #else
       am = ap = VL/wL + uxL;                /*  --  SHOCK SPEED  --  */
      #endif

      if (am >= 0.0) {                      /*  --  REGION L  --  */
        for (nv = NFLX; nv--;) ws[i][nv] = ql[nv];

      }else if (ap <= 0.0) {              /* REGION L1   */

        for (nv = NFLX; nv--;) ws[i][nv] = Ustar[nv];

      }else{                   /*  Solution is inside rarefaction fan, --> interpolate  */

        for (nv = NFLX; nv--;) {
          ws[i][nv] = (am*Ustar[nv] - ap*ql[nv])/(am - ap);
        }
      }
          
    } else {                /*  ----  Solution is sampled to the RIGHT of contact  ----  */

      EXPAND(        ,
             gR1 = gR*hR[i]/(gR*hR[i] + (p1 - pR)*(VR + wR*uxR));
             Ustar[VXt] = qr[VXt]*gR1;  ,
             Ustar[VXb] = qr[VXb]*gR1;)

      gR1        = GLORENTZ (Ustar,3);
      VR1        = VR - (u1 - uxR)*wR;
      Ustar[RHO] = MAX(small_rho,1.0/(VR1*gR1));

      #if INTERPOLATE_RAREFACTION == YES
       if (p1 < qr[PRS]){
         am  = RAREFACTION_SPEED (Ustar, 1);         /*  --  rarefaction tail  --  */
         ap  = RAREFACTION_SPEED (qr,    1);         /*  --  rarefaction head  --  */
         ap  = MAX(am, ap);
       }else{        
         am = ap = VR/wR + uxR;               /*  --  SHOCK SPEED  --  */
       }
      #else
       am = ap = VR/wR + uxR;                 /*  --  SHOCK SPEED  --  */
      #endif

      if (ap <= 0.0) {                        /*  --  REGION  R  --  */

        for (nv = NFLX; nv--;) ws[i][nv] = qr[nv];

      }else if (am >= 0.0){              /*   REGION R1   */

        for (nv = NFLX; nv--;) ws[i][nv] = Ustar[nv];

      }else{            /*  Solution is inside rarefaction fan, --> interpolate  */

        for (nv = NFLX; nv--;) {
          ws[i][nv] = (ap*Ustar[nv] - am*qr[nv])/(ap - am);
        }
      }
    }
 
/* ----------------------------------------------------------
                      Compute Fluxes               
   ---------------------------------------------------------- */

    #if USE_FOUR_VELOCITY == YES  /*  -- find 4-vel from the solution of the Riemann Problem --  */
     a0  = EXPAND( ws[i][VXn]*ws[i][VXn], 
                 + ws[i][VXt]*ws[i][VXt], 
                 + ws[i][VXb]*ws[i][VXb]);

     a0 = MIN(a0, 0.99999999);
     a0 = 1.0/sqrt(1.0 - a0);
     EXPAND(ws[i][VX1] *= a0;  ,
            ws[i][VX2] *= a0;  ,
            ws[i][VX3] *= a0;)
    #endif  
    PrimToCons (ws, us, i, i);
    SoundSpeed2 (ws, a2R, hR, i, i, FACE_CENTER, grid);
    Flux (us, ws, a2R, state->flux, state->press, i, i);
    MaxSignalSpeed (ws, a2R, cmin_loc, cmax_loc, i, i);
    a0 = MAX(fabs(cmax_loc[i]), fabs(cmin_loc[i]));
    cmax[i] = a0;
  }

/* ---------------------------------------------------------
    Compute HLL flux in those zones where the Riemann 
    Solver has failed (izone_fail[k] = 1, k > 0)
   --------------------------------------------------------- */

  for (k = 1; k <= nfail; k++){ 
    print1 ("! Failure in Riemann - substituting HLL flux: ");
    Where (izone_fail[k],NULL);

    HLL_Solver (state, izone_fail[k] - 2, izone_fail[k] + 3, 
                cmax, grid);

  }

  #if ARTIFICIAL_VISCOSITY == YES
   VISC_FLUX   (state, beg, end, grid);
  #endif
}

/* ############################################################# */
real GLORENTZ (real *U, int n)
/*
 #
 #               COMPUTE GAMMA LORENTZ FACTOR
 #
 ############################################################### */
{

   real wl2, beta_fix=0.9999, scrh;

   scrh = EXPAND(U[VX1]*U[VX1], + U[VX2]*U[VX2],  + U[VX3]*U[VX3]);

   if (scrh >= 1.0){

     print ("! u2 > 1 (%12.6e) in GLORENTZ\n", scrh);
     scrh = beta_fix/sqrt(scrh);
     EXPAND(U[VX1] *= scrh;  ,
            U[VX2] *= scrh;  ,
            U[VX3] *= scrh;)
     scrh = beta_fix*beta_fix;
     exit(1);
   }

   wl2 = 1.0/(1.0 - scrh);

   return(sqrt(wl2));
}

#undef  MAX_ITER


/* ############################################################### */
void SHOCK (real tau0, real u0, real p0, real g0, 
            real V0, real h0, real p1, real *u1,
            real *dudp, real *zeta, int istate)
/*
 #
 # Compute post shock quantities  u1, dudp, zeta, for a given
 # value of the post-shock pressure p1
 #
 #
 ############################################################### */
{
  real a,b,c, da,db,dc, dp;
  real tau1,h1,d_htau1,g1;
  real j2, dx;

/* ***********************************************************
        Use Taub Adiabat to find post-shock Enthalpy
        and mass flux; here j2 --> 1/(j)^2
   *********************************************************** */

  dp = p1 - p0;

  #if EOS == IDEAL
   a = 1.0 - dp/(gmmr*p1);
   b = 1.0 - a;
   c = -h0*(h0 + tau0*dp);
   
   h1 = 0.5/a*(-b + sqrt(b*b - 4.0*a*c));
   tau1 = (h1 - 1.0)/(gmmr*p1);
   g1   = 2.0*h1*gmmr/(2.0*h1 - 1.0);

   j2  = h0*gmmr*tau0 + (h1*tau1 + h0*tau0)*(1.0/(h0 + h1) - 1.0);
   j2 /= gmmr*p1;

   d_htau1  = (h1*tau1 + h0*tau0 - g1*h1*tau1);
   d_htau1 /= (g1*p1 - dp);
  #endif
  #if EOS == TAUB
   a = p0*(3.0*p1 + p0);
   b = -h0*h0*(3.0*p1 + 2.0*p0)
       -h0*tau0*(3.0*p1*p1 - 7.0*p1*p0 - 4.0*p0*p0) - dp;
   c = -h0*tau0*(h0*h0 + 2.0*h0*tau0*p1 + 2.0);
   
   dx = 2.0*c*dp/(-b + sqrt(b*b - 4.0*a*c*dp));/* = [h\tau] */

   g1 = 2.0*c/(-b + sqrt(b*b - 4.0*a*c*dp));/* = [h\tau]/[dp] */
   h1 = sqrt(h0*h0 + (dx + 2.0*h0*tau0)*dp);

   j2 = -g1;

   da = 3.0*p0;
   db = -h0*h0*3.0 - h0*tau0*(6.0*p1 - 7.0*p0) - 1.0;
   dc = c - h0*tau0*dp*2.0*h0*tau0;

   d_htau1 = -(da*dx*dx + db*dx + dc)/(2.0*a*dx + b);
  #endif

  g1    = ((real)istate)*sqrt(V0*V0 + (1.0 - u0*u0)*j2);
  *zeta = (V0*u0 + g1) / (1.0 - u0*u0);

/* ***********************************************************
                  Get post-shock velocity
   *********************************************************** */

  b   = 1.0/(h0*g0 + (p1 - p0)*((*zeta)*u0 + V0));
  *u1 = (h0*g0*u0 + (*zeta)*(p1 - p0)) * b;

/* ***********************************************************
                  Get du/dp for next iteration
   *********************************************************** */

  a     = -0.5*(d_htau1 + j2) / g1;
  *dudp = ((*zeta) + a - (*u1)*((*zeta)*u0 + V0 + a*u0)) * b;
}
/* ################################################################ */
real RAREFACTION_SPEED (real w[], int iside)
/*
 # 
 # PURPOSE
 #
 #   compute head or tail characteristic speeds enclosing
 #   the rarefaction fan
 #
 #
 # ARGUMENTS:
 #
 #  w    (IN)  :  a vector of primitive quantities containing 
 #                the three velocities.
 #  iside(IN)  :  an integer specifying a left (-1) or right (+1)
 #                rarefaction wave.
 #
 ################################################################## */
{
  int    nv;
  real   vx, vt2, vel2;
  real   sroot, delta2, cs2[1], h[1];
  static real **q;
  
  if (q == NULL){
    q = ARRAY_2D(10, NFLX, double);
  }
  
  for (nv = 0; nv < NFLX; nv++){
    q[0][nv] = w[nv];
  }
   
  SoundSpeed2(q, cs2, h, 0, 0, FACE_CENTER, NULL);
  
  vx   = w[VXn];
  vt2  = EXPAND(0.0, + w[VXt]*w[VXt], + w[VXb]*w[VXb]);
  vel2 = vx*vx + vt2;

  sroot = cs2[0]*(1.0 - vx*vx - vt2*cs2[0])*(1.0 - vel2);   /* this is eta/gamma */
    
  if (sroot < 0.0){
    print ("! sroot < 0 in RIemann \n");
    for (nv = 0; nv < NFLX; nv++){
      print ("%d  %12.6e  %12.6e\n",nv, qglob_l[nv], qglob_r[nv]);
    }
    QUIT_PLUTO(1);
  }
  sroot = sqrt(sroot);   /*this is eta/gamma */
  
  delta2   = 1.0 - vel2*cs2[0];
  return( (vx*(1.0 - cs2[0]) + (real)iside*sroot)/delta2);

}

#undef  MAX_ITER
#undef  small_p
#undef  small_rho
#undef  accuracy 
#undef  INTERPOLATE_RAREFACTION 

