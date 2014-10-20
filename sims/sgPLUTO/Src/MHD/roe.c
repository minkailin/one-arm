#include"pluto.h"

#define  sqrt_1_2  (0.70710678118654752440)
#define HLL_HYBRIDIZATION NO

/* **************************************************************************** */
void Roe_Solver (const State_1D *state, int beg, int end, 
                 real *cmax, Grid *grid)
/*
 *
 * PURPOSE
 *
 *   - Solve riemann problem for the adiabatic and isothermal
 *     MHD equations using the Roe linearized Riemann solver.
 *
 *     -----------------------------------------------------------
 *       Reference paper:
 *
 *       "Roe Matrices for Ideal MHD and Systematic Construction 
 *        of Roe Matrices for Systems of Conservation Laws"
 *
 *        P. Cargo, G. Gallice
 *        Journal of Computational Physics, 136, 446, (1997).
 *     -----------------------------------------------------------
 *
 *   - The isothermal version is recovered by taking the limit 
 *     of \bar{a}^2 for \gamma -> 1, which gives (page 451)
 *
 *         \bar{a}^2 -> g_isoSoundSpeed2 + X
 *  
 *     where X is defined as in the adiabatic case.
 *     Furthermore, the characteristic variables must 
 *     be modified by imposing zero jump across the entropy
 *     wave (first of Eq. 4.20, page 452), giving
 *
 *      dp = (\bar{a}^2 - X)*drho
 *
 *     and this condition must be used in the following 
 *     jumps, e.g., the term
 *
 *     (X*drho + dp)    becomes -->  (\bar{a}^2 * drho)
 *
 *
 * USEFUL SWITCH(ES):  
 *
 *   HLL_HYBRIDIZATION:  when set to YES, revert to the HLL solver
 *                       whenever an unphysical state appear
 *                       in the solution.
 *
 * LAST_MODIFIED
 *
 *   April 15th 2008, by Andrea Mignone  (mignone@to.astro.it)
 *
 *
 **************************************************************************** */
 {
  int  nv, i, j, k;
  int  ifail;
  real scrh0, scrh1, scrh2, scrh3;
  real rho, u, v, w, vel2, bx, by, bz, pr;
  real a2, a, ca2, cf2, cs2;
  real cs, ca, cf, b2;
  real S;
  real alpha_f, alpha_s, beta_y, beta_z;
  real dV[NFLX], dU[NFLX], *vL, *vR, *uL, *uR;
  real *SL, *SR;
  real Rc[NFLX][NFLX], eta[NFLX], lambda[NFLX];
  real alambda[NFLX], Uv[NFLX];

  real tau, sqrt_rho;
  real delta, delta_inv;
  
  real g1, sl, sr, H, HL, HR, Bx, By, Bz, X;
  real bt2, Btmag, sqr_rho_L, sqr_rho_R;

  static double *pR, *pL, *a2L, *a2R;
  static double **fL, **fR;
  static double **VL, **VR, **UL, **UR;
  real **bgf;
  real Us[NFLX];
  delta    = 1.e-6;

  if (fL == NULL){

    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);

    pR = ARRAY_1D(NMAX_POINT, double);
    pL = ARRAY_1D(NMAX_POINT, double);

    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);

    #ifdef GLM_MHD
     VL = ARRAY_2D(NMAX_POINT, NVAR, double);
     VR = ARRAY_2D(NMAX_POINT, NVAR, double);
     UL = ARRAY_2D(NMAX_POINT, NVAR, double);
     UR = ARRAY_2D(NMAX_POINT, NVAR, double);
    #endif
  }

  #if BACKGROUND_FIELD == YES
   print1 (" ! Background field not available for this solver\n");
   QUIT_PLUTO(1);
  #endif

  #ifdef GLM_MHD
   GLM_Solve (state, VL, VR, beg, end, grid);
   PrimToCons (VL, UL, beg, end);
   PrimToCons (VR, UR, beg, end);
  #else
   VL = state->vL; UL = state->uL;
   VR = state->vR; UR = state->uR;
  #endif

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (VL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux  (UL, VL, a2L, bgf, fL, pL, beg, end);
  Flux  (UR, VR, a2R, bgf, fR, pR, beg, end);

  SL = state->SL; SR = state->SR;

  #if EOS != ISOTHERMAL && EOS != BAROTROPIC
   g1 = g_gamma - 1.0;
  #endif

/* -------------------------------------------------
     Some eigenvectors components will always be 
     zero so set Rc = 0 initially  
   -------------------------------------------------- */
     
  for (k = NFLX; k--;  ) {
  for (j = NFLX; j--;  ) {
    Rc[k][j] = 0.0;
  }}

  for (i = beg; i <= end; i++) {

    vL = VL[i]; uL = UL[i];
    vR = VR[i]; uR = UR[i];

    #if SHOCK_FLATTENING == MULTID

  /* -- revert to HLL in proximity of strong shock -- */

     if (CheckZone(i, FLAG_HLL) || CheckZone(i+1,FLAG_HLL)){
       HLL_Speed (VL, VR, a2L, a2R, NULL, SL, SR, i, i);

       scrh0 = MAX(fabs(SL[i]), fabs(SR[i]));
       cmax[i] = scrh0;

       if (SL[i] > 0.0) {
         for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
         state->press[i] = pL[i];
       } else if (SR[i] < 0.0) {
         for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
         state->press[i] = pR[i];
       }else{
         scrh0 = 1.0/(SR[i] - SL[i]);
         for (nv = NFLX; nv--; ){
           state->flux[i][nv]  = SR[i]*SL[i]*(uR[nv] - uL[nv])
                              +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
           state->flux[i][nv] *= scrh0;
         }
         state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh0;
       }
       continue;
     }
    #endif

  /* -----------------------------------
       compute jumps in conservative 
       and primitive variables 
     ----------------------------------- */

    for (nv = 0; nv < NFLX; nv++) { 
      dV[nv] = vR[nv] - vL[nv];
      dU[nv] = uR[nv] - uL[nv];
    }

  /* ---------------------------------
         compute Roe averages 
     --------------------------------- */

    sqr_rho_L = sqrt(vL[RHO]);
    sqr_rho_R = sqrt(vR[RHO]);

    sl = sqr_rho_L/(sqr_rho_L + sqr_rho_R);
    sr = sqr_rho_R/(sqr_rho_L + sqr_rho_R);

/*      sl = sr = 0.5;    */
    
    rho = sr*vL[RHO] + sl*vR[RHO];

    tau      = 1.0/rho;
    sqrt_rho = sqrt(rho);

    EXPAND (u = sl*vL[VXn] + sr*vR[VXn];  ,
            v = sl*vL[VXt] + sr*vR[VXt];  ,
            w = sl*vL[VXb] + sr*vR[VXb];)

    EXPAND (Bx = sr*vL[BXn] + sl*vR[BXn];  ,
            By = sr*vL[BXt] + sl*vR[BXt];  ,
            Bz = sr*vL[BXb] + sl*vR[BXb];)
 
    S    = (Bx >= 0.0 ? 1.0 : -1.0);

    EXPAND(bx = Bx/sqrt_rho;  ,
           by = By/sqrt_rho;  ,
           bz = Bz/sqrt_rho; )
    
    bt2   = EXPAND(0.0  , + by*by, + bz*bz);
    b2    = bx*bx + bt2;
    Btmag = sqrt(bt2*rho);

    X  = EXPAND(dV[BXn]*dV[BXn], + dV[BXt]*dV[BXt], + dV[BXb]*dV[BXb]);
    X /= (sqr_rho_L + sqr_rho_R)*(sqr_rho_L + sqr_rho_R)*2.0;   

    scrh0 = EXPAND(u*dU[MXn],  + v*dU[MXt],  + w*dU[MXb]);
    scrh1 = EXPAND(Bx*dU[BXn], + By*dU[BXt], + Bz*dU[BXb]);

  /* -------------------------------
       compute enthalpy
     ------------------------------- */

    #if EOS == ISOTHERMAL 
     a2 = 0.5*(a2L[i] + a2R[i]) + X;  /* in most cases a2L = a2R
                                                       for isothermal MHD */
    #elif EOS == BAROTROPIC
     print ("! ROE_SOLVER: not implemented for barotropic EOS\n");
     QUIT_PLUTO(1);
    #elif EOS == IDEAL
     vel2  = EXPAND(u*u, + v*v, + w*w);
     HL = (uL[ENG] + pL[i])/vL[RHO];
     HR = (uR[ENG] + pR[i])/vR[RHO];
     H  = sl*HL + sr*HR;
     dV[PRS] = g1*((0.5*vel2 - X)*dV[RHO]
                 - scrh0 + dU[ENG] - scrh1);
        
     a2 = (2.0 - g_gamma)*X + g1*(H - 0.5*vel2 - b2);
    
     if (a2 < 0.0) {
      printf ("! Roe: a2 < 0.0 !! \n");
      Show(VL,i);
      Show(VR,i);
      QUIT_PLUTO(1);
     }      
    #endif
    
/* ------------------------------------------------------------
    Compute fast and slow magnetosonic speeds.

    The following expression appearing in the definitions
    of the fast magnetosonic speed 
    
     (a^2 - b^2)^2 + 4*a^2*bt^2 = (a^2 + b^2)^2 - 4*a^2*bx^2

    is always positive and avoids round-off errors.
   ------------------------------------------------------------ */
        
    scrh0 = a2 - b2;
    ca2   = bx*bx;
    scrh0 = scrh0*scrh0 + 4.0*bt2*a2;    
    scrh0 = sqrt(scrh0);    

    cf2 = 0.5*(a2 + b2 + scrh0); 
    cs2 = a2*ca2/cf2;   /* -- same as 0.5*(a2 + b2 - scrh0) -- */
    
    cf = sqrt(cf2);
    cs = sqrt(cs2);
    ca = sqrt(ca2);
    a  = sqrt(a2); 
    
    if (cf == cs) {
      alpha_f = 1.0;
      alpha_s = 0.0;
    }else if (a <= cs) {
      alpha_f = 0.0;
      alpha_s = 1.0;
    }else if (cf <= a){
      alpha_f = 1.0;
      alpha_s = 0.0;
    }else{ 
      scrh0   = 1.0/(cf2 - cs2);
      alpha_f = (a2  - cs2)*scrh0;
      alpha_s = (cf2 -  a2)*scrh0;
      alpha_f = MAX(0.0, alpha_f);
      alpha_s = MAX(0.0, alpha_s);
      alpha_f = sqrt(alpha_f);
      alpha_s = sqrt(alpha_s);
    }

    if (Btmag > 1.e-9) {
      SELECT(                     , 
             beta_y = DSIGN(By);  ,
             beta_y = By/Btmag; 
             beta_z = Bz/Btmag;)
    } else {
      SELECT(                       , 
             beta_y = 1.0;          ,
             beta_z = beta_y = sqrt_1_2;)
    }

  /* --------------------------------------------------------
      Compute non-zero entries of conservative
      eigenvectors (Rc), wave strength L*dU (=eta) for all 
      8 (or 7) waves.
      The expressions are given by eq. (4.18)--(4.21)    
      Fast and slow right (left) eigenvectors are divided by 
      rho (a^2). 
      
      NOTE: the expression on the paper has a typo in the 
            very last term of the energy component: 
            it should be + and not - !
     -------------------------------------------------------- */

  /* -----------------------
      FAST WAVE  (u - c_f) 
     ----------------------- */

    k = KFASTM;
    lambda[k] = u - cf;

    scrh0 = alpha_s*cs*S;
    scrh1 = EXPAND(0.0, + beta_y*dV[VXt], + beta_z*dV[VXb]);
    scrh2 = EXPAND(0.0, + beta_y*dV[BXt], + beta_z*dV[BXb]);
    scrh3 = EXPAND(0.0, + v*beta_y, + w*beta_z);

    Rc[RHO][k] = alpha_f;
    EXPAND( Rc[MXn][k] = alpha_f*lambda[k];         ,
            Rc[MXt][k] = alpha_f*v + scrh0*beta_y;  ,
            Rc[MXb][k] = alpha_f*w + scrh0*beta_z; ) 
    EXPAND(                                         ,                                
            Rc[BXt][k] = alpha_s*a*beta_y/sqrt_rho;  ,
            Rc[BXb][k] = alpha_s*a*beta_z/sqrt_rho; )

    #if EOS == IDEAL
     Rc[ENG][k] =   alpha_f*(H - b2 - u*cf) + scrh0*scrh3 
                 + alpha_s*a*Btmag/sqrt_rho;

     eta[k] = alpha_f*(X*dV[RHO] + dV[PRS]) 
              + rho*scrh0*scrh1
              - rho*alpha_f*cf*dV[VXn] + sqrt_rho*alpha_s*a*scrh2;
    #elif EOS == ISOTHERMAL
     eta[k] = alpha_f*(0.0*X + a2)*dV[RHO] 
              + rho*scrh0*scrh1
              - rho*alpha_f*cf*dV[VXn] + sqrt_rho*alpha_s*a*scrh2;
    #endif
    
    eta[k] *= 0.5/a2;

  /* -----------------------
      FAST WAVE  (u + c_f) 
     ----------------------- */

    k = KFASTP;
    lambda[k] = u + cf;

    Rc[RHO][k] = alpha_f;
    EXPAND( Rc[MXn][k] = alpha_f*lambda[k];         ,
            Rc[MXt][k] = alpha_f*v - scrh0*beta_y;  ,
            Rc[MXb][k] = alpha_f*w - scrh0*beta_z; ) 
    EXPAND(                              ,                                
            Rc[BXt][k] = Rc[BXt][KFASTM];  ,
            Rc[BXb][k] = Rc[BXb][KFASTM]; )

    #if EOS == IDEAL
     Rc[ENG][k] =   alpha_f*(H - b2 + u*cf) - scrh0*scrh3 
                 + alpha_s*a*Btmag/sqrt_rho;

     eta[k] = alpha_f*(X*dV[RHO] + dV[PRS]) 
              - rho*scrh0*scrh1
              + rho*alpha_f*cf*dV[VXn] + sqrt_rho*alpha_s*a*scrh2;
    #elif EOS == ISOTHERMAL
     eta[k] = alpha_f*(0.*X + a2)*dV[RHO] 
              - rho*scrh0*scrh1
              + rho*alpha_f*cf*dV[VXn] + sqrt_rho*alpha_s*a*scrh2;
    #endif

    eta[k] *= 0.5/a2;

  /* -----------------------
      Entropy wave  (u) 
     ----------------------- */

    #if EOS == IDEAL
     k = KENTRP;
     lambda[k] = u;

     Rc[RHO][k] = 1.0;
     EXPAND( Rc[MXn][k] = u; ,
             Rc[MXt][k] = v; ,
             Rc[MXb][k] = w; )
     Rc[ENG][k] = 0.5*vel2 + (g_gamma - 2.0)/g1*X;

     eta[k] = ((a2 - X)*dV[RHO] - dV[PRS])/a2;
    #endif

  /* --------------------------------------
         div.B wave  (u) 

      This wave exists when: 

       1) 8 wave formulation
       2) CT, since we
          always have 8 components, but it 
          carries zero jump.

     With Divergence_Cleaning, KDIVB is
     replaced by KPSI_GLMM, KPSI_GLMP.
     These two waves, however, should not
     enter in the Riemann solver since 
     the 2x2 linear system formed by (B,psi)
     has already been solved.
       
     -------------------------------------- */

    #ifdef GLM_MHD
     lambda[KPSI_GLMP] =  glm_ch;
     lambda[KPSI_GLMM] = -glm_ch;
     eta[KPSI_GLMP] = eta[KPSI_GLMM] = 0.0;
    #else
     k = KDIVB;
     lambda[k] = u;
     #if MHD_FORMULATION == EIGHT_WAVES
      Rc[BXn][k] = 1.0;
      eta[k]    = dU[BXn];
     #else
      Rc[BXn][k] = eta[k] = 0.0;
     #endif
    #endif
    
    #if COMPONENTS > 1    

   /* -----------------------
       SLOW WAVE  (u - c_s) 
      ----------------------- */

     k = KSLOWM;
     lambda[k] = u - cs;

     Rc[RHO][k] = alpha_s;
     EXPAND( Rc[MXn][k] = alpha_s*lambda[k];                ,
             Rc[MXt][k] = alpha_s*v - alpha_f*cf*beta_y*S;  ,
             Rc[MXb][k] = alpha_s*w - alpha_f*cf*beta_z*S; ) 
     EXPAND(                                           ,                                
             Rc[BXt][k] = - alpha_f*a*beta_y/sqrt_rho;  ,
             Rc[BXb][k] = - alpha_f*a*beta_z/sqrt_rho; )

     #if EOS == IDEAL
      Rc[ENG][k] =   alpha_s*(H - b2 - u*cs) - alpha_f*cf*S*scrh3
                  - alpha_f*a*Btmag/sqrt_rho; 

      eta[k] = alpha_s*(X*dV[RHO] + dV[PRS]) 
                - rho*alpha_f*cf*S*scrh1
                - rho*alpha_s*cs*dV[VXn] - sqrt_rho*alpha_f*a*scrh2;
     #elif EOS == ISOTHERMAL
      eta[k] = alpha_s*(0.*X + a2)*dV[RHO] 
               - rho*alpha_f*cf*S*scrh1
               - rho*alpha_s*cs*dV[VXn] - sqrt_rho*alpha_f*a*scrh2;
     #endif

     eta[k] *= 0.5/a2;

   /* -----------------------
       SLOW WAVE  (u + c_s) 
      ----------------------- */

     k = KSLOWP;
     lambda[k] = u + cs; 

     Rc[RHO][k] = alpha_s;
     EXPAND( Rc[MXn][k] = alpha_s*lambda[k];                ,
             Rc[MXt][k] = alpha_s*v + alpha_f*cf*beta_y*S;  ,
             Rc[MXb][k] = alpha_s*w + alpha_f*cf*beta_z*S; ) 
     EXPAND(                              ,                                
             Rc[BXt][k] = Rc[BXt][KSLOWM];  ,
             Rc[BXb][k] = Rc[BXb][KSLOWM]; )

     #if EOS == IDEAL
      Rc[ENG][k] =   alpha_s*(H - b2 + u*cs) + alpha_f*cf*S*scrh3
                  - alpha_f*a*Btmag/sqrt_rho;

      eta[k] = alpha_s*(X*dV[RHO] + dV[PRS]) 
               + rho*alpha_f*cf*S*scrh1
               + rho*alpha_s*cs*dV[VXn] - sqrt_rho*alpha_f*a*scrh2; 
     #elif EOS == ISOTHERMAL
      eta[k] = alpha_s*(0.*X + a2)*dV[RHO] 
               + rho*alpha_f*cf*S*scrh1
               + rho*alpha_s*cs*dV[VXn] - sqrt_rho*alpha_f*a*scrh2; 
     #endif

     eta[k] *= 0.5/a2;

    #endif

    #if COMPONENTS == 3

   /* ------------------------
       Alfven WAVE  (u - c_a) 
      ------------------------ */

     k = KALFVM;
     lambda[k] = u - ca;

     Rc[MXt][k] = - rho*beta_z;  
     Rc[MXb][k] = + rho*beta_y;
     Rc[BXt][k] = - S*sqrt_rho*beta_z;   
     Rc[BXb][k] =   S*sqrt_rho*beta_y;
     #if EOS == IDEAL
      Rc[ENG][k] = - rho*(v*beta_z - w*beta_y);
     #endif

     eta[k] = + beta_y*dV[VXb] - beta_z*dV[VXt] 
              + S/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

     eta[k] *= 0.5;

   /* -----------------------
       Alfven WAVE  (u + c_a) 
      ----------------------- */

     k = KALFVP;
     lambda[k] = u + ca;

     Rc[MXt][k] = - Rc[MXt][KALFVM];  
     Rc[MXb][k] = - Rc[MXb][KALFVM];
     Rc[BXt][k] =   Rc[BXt][KALFVM];   
     Rc[BXb][k] =   Rc[BXb][KALFVM];
     #if EOS == IDEAL
      Rc[ENG][k] = - Rc[ENG][KALFVM];
     #endif

     eta[k] = - beta_y*dV[VXb] + beta_z*dV[VXt] 
              + S/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

     eta[k] *= 0.5;
    #endif

   /* -----------------------------------------
          Compute maximum signal velocity
      ----------------------------------------- */

    cmax[i] = fabs (u) + cf;
    g_maxMach = MAX (fabs (u / a), g_maxMach);
    for (k = 0; k < NFLX; k++) alambda[k] = fabs(lambda[k]);

   /* --------------------------------
              Entropy Fix 
      -------------------------------- */
      
    if (alambda[KFASTM] < 0.5*delta) {
      alambda[KFASTM] = lambda[KFASTM]*lambda[KFASTM]/delta + 0.25*delta;
    }
    if (alambda[KFASTP] < 0.5*delta) {
      alambda[KFASTP] = lambda[KFASTP]*lambda[KFASTP]/delta + 0.25*delta;
    }
    #if COMPONENTS > 1
     if (alambda[KSLOWM] < 0.5*delta) {
       alambda[KSLOWM] = lambda[KSLOWM]*lambda[KSLOWM]/delta + 0.25*delta;
     }
     if (alambda[KSLOWP] < 0.5*delta) {
       alambda[KSLOWP] = lambda[KSLOWP]*lambda[KSLOWP]/delta + 0.25*delta; 
     }

    #endif
   
  /*  ---------------------------------
         Compute Roe numerical flux 
      --------------------------------- */

    for (nv = 0; nv < NFLX; nv++) {
      scrh0 = 0.0;
      for (k = 0; k < NFLX; k++) {
        scrh0 += alambda[k]*eta[k]*Rc[nv][k];
      }
      state->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - scrh0);

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

#if MHD_FORMULATION == CONSTRAINED_TRANSPORT && CT_EMF_AVERAGE == RIEMANN_2D
 state->pnt_flx[i][nv] =  0.5*(vL[nv]*vL[VXn] + vR[nv]*vR[VXn]); 
 state->dff_flx[i][nv] = -0.5*scrh0;
#endif

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
      
    }
    state->press[i] = 0.5*(pL[i] + pR[i]);

    #if CHECK_ROE_MATRIX == YES

     /* -----------------------------------
         check the Roe matrix condition,
 
             FR - FL = A*(UR - UL)

         where A*(UR - UL) = R*lambda*eta
        ----------------------------------- */

     for (nv = 0; nv < NFLX; nv++){
       dV[nv] = fR[i][nv] - fL[i][nv]; 
       if (nv == MXn) dV[MXn] += pR[i] - pL[i];
       for (k = 0; k < NFLX; k++){
         dV[nv] -= Rc[nv][k]*eta[k]*lambda[k];
       }
       if (fabs(dV[nv]) > 1.e-4){
         printf (" ! Roe matrix condition not satisfied, var = %d\n", nv);
         printf (" ! Err = %12.6e\n",dV[nv]); 
         Show(VL, i);
         Show(VR, i);
         exit(1);
       }
     } 
    #endif

  /* -------------------------------------
        Save max and min Riemann fan 
        speeds for EMF computation
     ------------------------------------- */

    SL[i] = lambda[KFASTM];
    SR[i] = lambda[KFASTP];

    #if HLL_HYBRIDIZATION == YES

  /* ------------------------------------------------------
      Hybridize with HLL solver: replace occurences 
      of unphysical states (p < 0, rho < 0) with 
      HLL Flux. Reference:
      
      "A Positive Conservative Method for MHD based
      based on HLL and Roe methods"
      P. Janhunen, JCP (2000), 160, 649

     ------------------------------------------------------ */

     if (SL[i] < 0.0 && SR[i] > 0.0){

       ifail = 0;    

      /* -----------------------
           check left state 
         ----------------------- */

       #if EOS == ISOTHERMAL
        Uv[RHO] = uL[RHO] + (state->flux[i][RHO] - fL[i][RHO])/SL[i];        
        ifail  = (Uv[RHO] < 0.0);
       #else
        for (nv = NFLX; nv--; ){
          Uv[nv] = uL[nv] + (state->flux[i][nv] - fL[i][nv])/SL[i];        
        }
        Uv[MXn] += (state->press[i] - pL[i])/SL[i];    
 
        scrh0 = EXPAND(Uv[MX1]*Uv[MX1], + Uv[MX2]*Uv[MX2], + Uv[MX3]*Uv[MX3]);
        scrh1 = EXPAND(Uv[BX1]*Uv[BX1], + Uv[BX2]*Uv[BX2], + Uv[BX3]*Uv[BX3]);    
        scrh2 = Uv[ENG] - 0.5*scrh0/Uv[RHO] - 0.5*scrh1;
        ifail = (scrh2 < 0.0) || (Uv[RHO] < 0.0);
       #endif

      /* -----------------------
           check right state 
         ----------------------- */

       #if EOS == ISOTHERMAL
        Uv[RHO] = uR[RHO] + (state->flux[i][RHO] - fR[i][RHO])/SR[i];
        ifail  = (Uv[RHO] < 0.0);
       #else
        for (nv = NFLX; nv--;  ){
          Uv[nv] = uR[nv] + (state->flux[i][nv] - fR[i][nv])/SR[i];
        }
        Uv[MXn] += (state->press[i] - pR[i])/SR[i];

        scrh0 = EXPAND(Uv[MX1]*Uv[MX1], + Uv[MX2]*Uv[MX2], + Uv[MX3]*Uv[MX3]);
        scrh1 = EXPAND(Uv[BX1]*Uv[BX1], + Uv[BX2]*Uv[BX2], + Uv[BX3]*Uv[BX3]);    
        scrh2 = Uv[ENG] - 0.5*scrh0/Uv[RHO] - 0.5*scrh1;
        ifail = (scrh2 < 0.0) || (Uv[RHO] < 0.0);
       #endif

       #if DIMENSIONS > 1
   
     /* ---------------------------------------------
          use the HLL flux function if the interface 
          lies within a strong shock.
          The effect of this switch is visible
          in the Mach reflection test.
        --------------------------------------------- */

       #if EOS == ISOTHERMAL
        scrh0  = fabs(vL[RHO] - vR[RHO]);
        scrh0 /= MIN(vL[RHO], vR[RHO]);
       #else       
        scrh0  = fabs(vL[PRS] - vR[PRS]);
        scrh0 /= MIN(vL[PRS], vR[PRS]);
       #endif
       if (scrh0 > 1.0 && (vR[VXn] < vL[VXn])) ifail = 1;

       #endif
      
       if (ifail){
         scrh0 = 1.0/(SR[i] - SL[i]);
    
         for (nv = 0; nv < NFLX; nv++) {
           state->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                                SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
           state->flux[i][nv] *= scrh0;
         }
         state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh0;
       }
     }
    #endif

  }

/* --------------------------------------------------------
              initialize source term
   -------------------------------------------------------- */
  
  #if MHD_FORMULATION == EIGHT_WAVES
   ROE_DIVB_SOURCE (state, beg + 1, end, grid);
  #endif

}
#undef sqrt_1_2
#undef HLL_HYBRIDIZATION
