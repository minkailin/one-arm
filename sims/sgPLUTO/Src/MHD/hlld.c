#include"pluto.h"

#define VERIFY_CONSISTENCY_CONDITION NO

#if EOS == IDEAL
/* ***************************************************************************** */
void HLLD_Solver (const State_1D *state, int beg, int end, 
                  real *cmax, Grid *grid)
/* 
 *
 *
 * NAME
 *
 *   HLLD_SOLVER 
 *
 *
 * PURPOSE
 *
 *   Solve riemann problem for the MHD equations using the 
 *   four-state HLLD 
 * 
 *   Reference:    "A muLti-state HLL approximate Riemann Solver for 
 *                  ideal MHD", Miyoshi, T., Kusano, K., JCP 2005 
 *
 * ARGUMENTS
 *
 *   vL (IN)       1-D array of left-edge primitive values at i+1/2
 *   vR (IN)       1-D array of right-edge primitive values at i+1/2
 *   flux (OUT)    1-D array of numerical fluxes, not including pressuRe
 *   grid (IN)     Array of grids
 *   cmax(OUT)     Array of maximum characteristic speeds in this direction
 *
 * SWITCHES
 *
 *
 * LAST_MODIFIED
 *
 *   June 8, 2007 by Andrea Mignone  (mignone@to.astro.it)
 *              
 *
 ******************************************************************************** */
{
  int    nv, i;
  int    revert_to_hllc;
  real   scrh, Uhll[NFLX];

  real usL[NFLX], ussl[NFLX];
  real usR[NFLX], ussr[NFLX];
  
  real vsL, wsL, scrhL, S1L, sqrL, duL;
  real vsR, wsR, scrhR, S1R, sqrR, duR;
  real Bx, SM, sBx, pts;
  real vss, wss;
  real *vL, *vR, *uL, *uR, *SL, *SR;
  static real *ptL, *ptR;
  static real **fL, **fR;
  static double **VL, **VR, **UL, **UR;
  static double *a2L, *a2R;

  if (fL == NULL){
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);

    ptR = ARRAY_1D(NMAX_POINT, double);
    ptL = ARRAY_1D(NMAX_POINT, double);

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
   print ("! Background field splitting not allowed with HLLD solver\n");
   QUIT_PLUTO(1);
  #endif
  #if MHD_FORMULATION == EIGHT_WAVES
   print ("! hlld Riemann solver does not work with Powell's 8-wave\n");
   QUIT_PLUTO(1);
  #endif
 
  #ifdef GLM_MHD
/*
double **uL0, **uR0;
static double *BxL, *BxR, *psiL, *psiR;
if (BxL == NULL){
 BxL = ARRAY_1D(NMAX_POINT, double);
 BxR = ARRAY_1D(NMAX_POINT, double);
 psiL = ARRAY_1D(NMAX_POINT, double);
 psiR = ARRAY_1D(NMAX_POINT, double);
}
uL0 = state->uL;  
uR0 = state->uR;
for (i = beg; i <= end; i++){
 BxL[i]  = state->vL[i][BXn];
 BxR[i]  = state->vR[i][BXn];
 psiL[i] = state->vL[i][PSI_GLM];
 psiR[i] = state->vL[i][PSI_GLM];
}
state->uL = UL;
state->uR = UR;

   GLM_Solve (state, beg, end, grid);
*/

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

  Flux (UL, VL, a2L, NULL, fL, ptL, beg, end);
  Flux (UR, VR, a2R, NULL, fR, ptR, beg, end);

/* ----------------------------------------
     get max and min signal velocities
   ---------------------------------------- */
             
  SL = state->SL; SR = state->SR;
  HLL_Speed (VL, VR, a2L, a2R, NULL, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {
    
  /* ----------------------------------------
      get max propagation speed for dt comp.
     ---------------------------------------- */             

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    vL = VL[i]; uL = UL[i];
    vR = VR[i]; uR = UR[i];

/* ---------------------------------------------------------- 
                COMPUTE FLUXES and STATES
   ---------------------------------------------------------- */

    if (SL[i] >= 0.0){                     /*  ----  Region L  ---- */

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
      state->press[i] = ptL[i];

    }else if (SR[i] <= 0.0) {              /*  ----  Region R  ---- */

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
      state->press[i] = ptR[i];
 
    } else {

      #if SHOCK_FLATTENING == MULTID

      /* -- revert to HLL in proximity of strong shocks -- */

       if (CheckZone(i, FLAG_HLL) || CheckZone(i+1, FLAG_HLL)){
         scrh = 1.0/(SR[i] - SL[i]);
         for (nv = NFLX; nv--; ){
           state->flux[i][nv]  = SR[i]*SL[i]*(uR[nv] - uL[nv])
                              +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
           state->flux[i][nv] *= scrh;
         }
         state->press[i] = (SR[i]*ptL[i] - SL[i]*ptR[i])*scrh;
         continue;
       }
      #endif

    /* ---------------------------
              Compute U*  
       --------------------------- */

      scrh = 1.0/(SR[i] - SL[i]);
      Bx   = (SR[i]*vR[BXn] - SL[i]*vL[BXn])*scrh; 
      sBx  = (Bx > 0.0 ? 1.0 : -1.0);

      duL  = SL[i] - vL[VXn];
      duR  = SR[i] - vR[VXn];

      scrh = 1.0/(duR*uR[RHO] - duL*uL[RHO]);
      SM   = (duR*uR[MXn] - duL*uL[MXn] - ptR[i] + ptL[i])*scrh;

      pts  = duR*uR[RHO]*ptL[i] - duL*uL[RHO]*ptR[i] + 
             vL[RHO]*vR[RHO]*duR*duL*(vR[VXn]- vL[VXn]);
      pts *= scrh;

      usL[RHO] = uL[RHO]*duL/(SL[i] - SM);
      usR[RHO] = uR[RHO]*duR/(SR[i] - SM);

      sqrL = sqrt(usL[RHO]);
      sqrR = sqrt(usR[RHO]);

      S1L = SM - fabs(Bx)/sqrL;
      S1R = SM + fabs(Bx)/sqrR;

    /* ---------------------------------------------
        When S1L -> SL or S1R -> SR a degeneracy 
        occurs. Although Miyoshi & Kusano say that 
        no jump exists, we don't think this is 
        actually true. Indeed, vy*, vz*, By*, Bz* 
        cannote be solved independently. 
        In this case we revert to the HLLC solver 
        of Li (2005), except for the term v.B in the
        * region, which we compute in our own way.
        Note, that by comparing the expressions of 
        Li (2005) and Miyoshi & Kusano (2005), the 
        only change involves a re-definition
        of By* and Bz* in terms of By(HLL), Bz(HLL).
       --------------------------------------------- */

      revert_to_hllc = 0;

      if ( (S1L - SL[i]) <  1.e-4*(SM - SL[i]) ) revert_to_hllc = 1;
      if ( (S1R - SR[i]) > -1.e-4*(SR[i] - SM) ) revert_to_hllc = 1;

      if (revert_to_hllc){

        scrh = 1.0/(SR[i] - SL[i]);
        for (nv = NFLX; nv--; ){  
          Uhll[nv]  = SR[i]*uR[nv] - SL[i]*uL[nv] + fL[i][nv] - fR[i][nv];
          Uhll[nv] *= scrh;
        }

        EXPAND(usL[BXn] = usR[BXn] = Uhll[BXn];   ,
               usL[BXt] = usR[BXt] = Uhll[BXt];   ,
               usL[BXb] = usR[BXb] = Uhll[BXb];)
 
        S1L = S1R = SM; /* region ** should never be computed since */ 
                        /* fluxes are given in terms of UL* and UR* */

      }else{

        scrhL = (uL[RHO]*duL*duL - Bx*Bx)/(uL[RHO]*duL*(SL[i] - SM) - Bx*Bx);
        scrhR = (uR[RHO]*duR*duR - Bx*Bx)/(uR[RHO]*duR*(SR[i] - SM) - Bx*Bx);
 
        EXPAND(usL[BXn]  = Bx;            ,
               usL[BXt]  = uL[BXt]*scrhL;  ,
               usL[BXb]  = uL[BXb]*scrhL;)           

        EXPAND(usR[BXn]  = Bx;            ,
               usR[BXt]  = uR[BXt]*scrhR;  ,
               usR[BXb]  = uR[BXb]*scrhR;)           
      }

      scrhL = Bx/(uL[RHO]*duL);
      scrhR = Bx/(uR[RHO]*duR);

      EXPAND(                                       ;  ,
             vsL = vL[VXt] - scrhL*(usL[BXt] - uL[BXt]);
             vsR = vR[VXt] - scrhR*(usR[BXt] - uR[BXt]);  ,

             wsL = vL[VXb] - scrhL*(usL[BXb] - uL[BXb]);
             wsR = vR[VXb] - scrhR*(usR[BXb] - uR[BXb]); )
         
      EXPAND(usL[MXn] = usL[RHO]*SM; 
             usR[MXn] = usR[RHO]*SM;   ,
    
             usL[MXt] = usL[RHO]*vsL;
             usR[MXt] = usR[RHO]*vsR;  ,

             usL[MXb] = usL[RHO]*wsL;
             usR[MXb] = usR[RHO]*wsR;)

      scrhL  = EXPAND(vL[VXn]*Bx, + vL[VXt]*uL[BXt] , + vL[VXb]*uL[BXb]);
      scrhL -= EXPAND(    SM*Bx, +    vsL*usL[BXt], +    wsL*usL[BXb]);
     
      usL[ENG]  = duL*uL[ENG] - ptL[i]*vL[VXn] + pts*SM + Bx*scrhL;
      usL[ENG] /= SL[i] - SM;

      scrhR  = EXPAND(vR[VXn]*Bx, + vR[VXt]*uR[BXt] , + vR[VXb]*uR[BXb]);
      scrhR -= EXPAND(    SM*Bx, +    vsR*usR[BXt], +    wsR*usR[BXb]);
     
      usR[ENG]  = duR*uR[ENG] - ptR[i]*vR[VXn] + pts*SM + Bx*scrhR;
      usR[ENG] /= SR[i] - SM;

      #ifdef GLM_MHD
       usL[PSI_GLM] = usR[PSI_GLM] = vL[PSI_GLM];
      #endif

  /* ------------------------------
         compute HLLD flux 
     ------------------------------ */

      if (S1L >= 0.0){       /*  ----  Region L*  ---- */

        for (nv = NFLX; nv--; ){
          state->flux[i][nv] = fL[i][nv] + SL[i]*(usL[nv] - uL[nv]);
        }
        state->press[i] = ptL[i];

      }else if (S1R <= 0.0) {    /*  ----  Region R*  ---- */
    
        for (nv = NFLX; nv--; ){
          state->flux[i][nv] = fR[i][nv] + SR[i]*(usR[nv] - uR[nv]);
        }
        state->press[i] = ptR[i];
         
      } else {   /* -- This state exists only if B_x != 0 -- */
      
  /* ---------------------------
           Compute U**
     --------------------------- */

        ussl[RHO] = usL[RHO];
        ussr[RHO] = usR[RHO];
 
        EXPAND(                           ,
       
               vss  = sqrL*vsL + sqrR*vsR + (usR[BXt] - usL[BXt])*sBx;       
               vss /= sqrL + sqrR;        ,
            
               wss  = sqrL*wsL + sqrR*wsR + (usR[BXb] - usL[BXb])*sBx;
               wss /= sqrL + sqrR;)
           
        EXPAND(ussl[MXn] = ussl[RHO]*SM;
               ussr[MXn] = ussr[RHO]*SM;    ,
     
               ussl[MXt] = ussl[RHO]*vss;
               ussr[MXt] = ussr[RHO]*vss;  ,
           
               ussl[MXb] = ussl[RHO]*wss;
               ussr[MXb] = ussr[RHO]*wss;)           
    
        EXPAND(ussl[BXn] = ussr[BXn] = Bx;   ,

               ussl[BXt]  = sqrL*usR[BXt] + sqrR*usL[BXt] + sqrL*sqrR*(vsR - vsL)*sBx;
               ussl[BXt] /= sqrL + sqrR;        
               ussr[BXt]  = ussl[BXt];        ,
           
               ussl[BXb]  = sqrL*usR[BXb] + sqrR*usL[BXb] + sqrL*sqrR*(wsR - wsL)*sBx;
               ussl[BXb] /= sqrL + sqrR;        
               ussr[BXb]  = ussl[BXb];)
          
        scrhL  = EXPAND(SM*Bx, +  vsL*usL [BXt], +  wsL*usL [BXb]);
        scrhL -= EXPAND(SM*Bx, +  vss*ussl[BXt], +  wss*ussl[BXb]);

        scrhR  = EXPAND(SM*Bx, +  vsR*usR [BXt], +  wsR*usR [BXb]);
        scrhR -= EXPAND(SM*Bx, +  vss*ussr[BXt], +  wss*ussr[BXb]);

        ussl[ENG] = usL[ENG] - sqrL*scrhL*sBx;
        ussr[ENG] = usR[ENG] + sqrR*scrhR*sBx;

        #ifdef GLM_MHD
         ussl[PSI_GLM] = ussr[PSI_GLM] = vL[PSI_GLM];
        #endif
    
  /* --------------------------------------
      verify consistency condition 
     -------------------------------------- */

/*
      for (nv = 0; nv < NFLX; nv++){
        scrh = (SR[i] - S1R)*usR[nv]  + (S1R - SM)*ussr[nv] +
               (SM - S1L)*ussl[nv] + (S1L - SL[i])*usL[nv] -
               SR[i]*UR[nv] + SL[i]*UL[nv] + fr[i][nv] - fl[i][nv];

        if (fabs(scrh) > 1.e-2){
          printf (" ! Consistency condition violated, pt %d, nv %d, %12.6e \n", 
                   i,nv,scrh);
          printf ("%f %f %f %f %f %f\n",UL[BXn], usL[BXn], ussl[BXn],
                                   ussr[BXn], usR[BXn], UR[BXn]);
          printf ("%f %f %f %f %f\n",SL[i], S1L, SM, S1R, SR[i]);
          printf ("%f %f  %d\n",fr[i][nv],fl[i][nv],BXn);
          exit(1);
        }
      }
*/

        if (SM >= 0.0){           /*  ----  Region L**  ---- */

          for (nv = NFLX; nv--; ){
            state->flux[i][nv] = fL[i][nv] + S1L*(ussl[nv] - usL[nv])
                                           + SL[i]*(usL[nv]  - uL[nv]);
          }
          state->press[i] = ptL[i];
        }else{                   /*  ----  Region R**  ---- */

          for (nv = NFLX; nv--; ){
            state->flux[i][nv] = fR[i][nv] + S1R*(ussr[nv] - usR[nv])
                                           + SR[i]*(usR[nv]  - uR[nv]);
          }
          state->press[i] = ptR[i];
        }
      }
    } 
  }
/*
#ifdef GLM_MHD
state->uL = uL0;
state->uR = uR0;
for (i = beg; i <= end; i++){
 state->vL[i][BXn] = BxL[i];
 state->vR[i][BXn] = BxR[i];
 state->vL[i][PSI_GLM] = psiL[i];
 state->vL[i][PSI_GLM] = psiR[i];
}
#endif
*/
}
#endif

#if EOS == ISOTHERMAL
/* ***************************************************************************** */
void HLLD_Solver (const State_1D *state, int beg, int end, 
                  real *cmax, Grid *grid)
/* 
 *
 *
 * NAME
 *
 *   HLLD_SOLVER 
 *
 *
 * PURPOSE
 *
 *   Solve riemann problem for the isothermal MHD equations using 
 *   the threee-state HLLD 
 * 
 *   Reference:    "A simple and accurate Riemann Solver for 
 *                  isothermal MHD"
 *                 A. Mignone, JCP (2005)
 *
 * ARGUMENTS
 *
 *
 * SWITCHES
 *
 *
 *
 * LAST_MODIFIED
 *
 *   June 8, 2007 by Andrea Mignone  (mignone@to.astro.it)
 *              
 *
 ******************************************************************************** */
{
  int  nv, i;
  int  revert_to_hll;
  double scrh;
  double usL[NFLX], *SL;
  double usR[NFLX], *SR;
  double usc[NFLX];
  
  double scrhL, S1L, duL;
  double scrhR, S1R, duR;
  double Bx, SM, sBx, rho, sqrho;

  double *vL, *vR, *uL, *uR;
  static double *ptL, *ptR, *a2L, *a2R;
  static double **fL, **fR;
  static double **VL, **VR, **UL, **UR;

  if (fL == NULL){
    fL  = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR  = ARRAY_2D(NMAX_POINT, NFLX, double);

    ptR = ARRAY_1D(NMAX_POINT, double);
    ptL = ARRAY_1D(NMAX_POINT, double);

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
   print ("! Background field splitting not allowed with HLLD solver\n");
   QUIT_PLUTO(1);
  #endif
  #if MHD_FORMULATION == EIGHT_WAVES
   print ("! hlld Riemann solver does not work with Powell\n");
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

  Flux (UL, VL, a2L, NULL, fL, ptL, beg, end);
  Flux (UR, VR, a2R, NULL, fR, ptR, beg, end);

/* ----------------------------------------
      get max and min signal velocities
   ---------------------------------------- */
             
  SL = state->SL; SR = state->SR;
  HLL_Speed (VL, VR, a2L, a2R, NULL, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {
    
  /* ----------------------------------------
      get max propagation speed for dt comp.
     ---------------------------------------- */             

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    vL = VL[i]; uL = UL[i];
    vR = VR[i]; uR = UR[i];

/* ---------------------------------------------------------- 
                COMPUTE FLUXES and STATES
   ---------------------------------------------------------- */

    if (SL[i] >= 0.0){                     /*  ----  Region L  ---- */

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
      state->press[i] = ptL[i];

    }else if (SR[i] <= 0.0) {              /*  ----  Region R   ---- */

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
      state->press[i] = ptR[i];
 
    } else {

      scrh = 1.0/(SR[i] - SL[i]);
      duL = SL[i] - vL[VXn];
      duR = SR[i] - vR[VXn];

      Bx   = (SR[i]*vR[BXn] - SL[i]*vL[BXn])*scrh; 

      rho                = (uR[RHO]*duR - uL[RHO]*duL)*scrh;
      state->flux[i][RHO] = (SL[i]*uR[RHO]*duR - SR[i]*uL[RHO]*duL)*scrh;
           
  /* ---------------------------
          compute S*
     --------------------------- */

      sqrho = sqrt(rho);

      SM  = state->flux[i][RHO]/rho;
      S1L = SM - fabs(Bx)/sqrho;
      S1R = SM + fabs(Bx)/sqrho;

    /* ---------------------------------------------
        Prevent degeneracies when S1L -> SL or 
        S1R -> SR. Revert to HLL if necessary.
       --------------------------------------------- */

      revert_to_hll = 0;

      if ( (S1L - SL[i]) <  1.e-4*(SR[i] - SL[i]) ) revert_to_hll = 1;
      if ( (S1R - SR[i]) > -1.e-4*(SR[i] - SL[i]) ) revert_to_hll = 1;

      if (revert_to_hll){
        scrh = 1.0/(SR[i] - SL[i]);
        for (nv = NFLX; nv--; ){
          state->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                               SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
          state->flux[i][nv] *= scrh;
        }
        state->press[i] = (SR[i]*ptL[i] - SL[i]*ptR[i])*scrh;
        continue;
      }

      state->flux[i][MXn] = (SR[i]*fL[i][MXn] - SL[i]*fR[i][MXn] 
                            + SR[i]*SL[i]*(uR[MXn] - uL[MXn]))*scrh;

      state->press[i]    = (SR[i]*ptL[i] - SL[i]*ptR[i])*scrh;
      #ifdef GLM_MHD
       state->flux[i][BXn]      = fL[i][BXn];
       state->flux[i][PSI_GLM] = fL[i][PSI_GLM];
      #else
       state->flux[i][BXn] = SR[i]*SL[i]*(uR[BXn] - uL[BXn])*scrh;
      #endif

  /* ---------------------------
             Compute U*  
     --------------------------- */
       
      scrhL = 1.0/((SL[i] - S1L)*(SL[i] - S1R));
      scrhR = 1.0/((SR[i] - S1L)*(SR[i] - S1R));

      EXPAND(                                                    ;  ,
             usL[MXt] = rho*vL[VXt] - Bx*uL[BXt]*(SM - vL[VXn])*scrhL;
             usR[MXt] = rho*vR[VXt] - Bx*uR[BXt]*(SM - vR[VXn])*scrhR;  ,

             usL[MXb] = rho*vL[VXb] - Bx*uL[BXb]*(SM - vL[VXn])*scrhL;
             usR[MXb] = rho*vR[VXb] - Bx*uR[BXb]*(SM - vR[VXn])*scrhR;)
         
      EXPAND(                                                   ; ,
             usL[BXt] = uL[BXt]/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL; 
             usR[BXt] = uR[BXt]/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR; ,

             usL[BXb] = uL[BXb]/rho*(uL[RHO]*duL*duL - Bx*Bx)*scrhL;           
             usR[BXb] = uR[BXb]/rho*(uR[RHO]*duR*duR - Bx*Bx)*scrhR;)           


      if (S1L >= 0.0){       /*  ----  Region L*  ---- */

        EXPAND(                                                    ;  ,
          state->flux[i][MXt] = fL[i][MXt] + SL[i]*(usL[MXt] - uL[MXt]);  ,
          state->flux[i][MXb] = fL[i][MXb] + SL[i]*(usL[MXb] - uL[MXb]);  
        ) 
        EXPAND(                                                    ;  ,
          state->flux[i][BXt] = fL[i][BXt] + SL[i]*(usL[BXt] - uL[BXt]);  ,
          state->flux[i][BXb] = fL[i][BXb] + SL[i]*(usL[BXb] - uL[BXb]);  
        ) 

      }else if (S1R <= 0.0) {    /*  ----  Region R*  ---- */
    
        EXPAND(                                                    ;  ,
          state->flux[i][MXt] = fR[i][MXt] + SR[i]*(usR[MXt] - uR[MXt]);  ,
          state->flux[i][MXb] = fR[i][MXb] + SR[i]*(usR[MXb] - uR[MXb]);  
        ) 
        EXPAND(                                                    ;  ,
          state->flux[i][BXt] = fR[i][BXt] + SR[i]*(usR[BXt] - uR[BXt]);  ,
          state->flux[i][BXb] = fR[i][BXb] + SR[i]*(usR[BXb] - uR[BXb]);  
        ) 
         
      } else {
                      
       /* ---------------------------
               Compute U** = Uc
          --------------------------- */

        sBx = (Bx > 0.0 ? 1.0 : -1.0);

        EXPAND(                                                                 ;  ,
               usc[MXt] = 0.5*(usR[MXt] + usL[MXt] + (usR[BXt] - usL[BXt])*sBx*sqrho);  ,     
               usc[MXb] = 0.5*(usR[MXb] + usL[MXb] + (usR[BXb] - usL[BXb])*sBx*sqrho);)
           
        EXPAND(                                                                 ;  ,
               usc[BXt] = 0.5*(usR[BXt] + usL[BXt] + (usR[MXt] - usL[MXt])*sBx/sqrho);  ,
               usc[BXb] = 0.5*(usR[BXb] + usL[BXb] + (usR[MXb] - usL[MXb])*sBx/sqrho);)


        EXPAND(                                            ;  ,
               state->flux[i][MXt] = usc[MXt]*SM - Bx*usc[BXt];  ,
               state->flux[i][MXb] = usc[MXb]*SM - Bx*usc[BXb]; )
        EXPAND(                                                ;  ,
               state->flux[i][BXt] = usc[BXt]*SM - Bx*usc[MXt]/rho;  ,
               state->flux[i][BXb] = usc[BXb]*SM - Bx*usc[MXb]/rho;)

    /* --------------------------------------
          verify consistency condition 
       -------------------------------------- */

        #if VERIFY_CONSISTENCY_CONDITION == YES
         for (nv = NFLX; nv--; ){
           if (nv == DN || nv == MXn || nv == BXn) continue;
           scrh = (S1L - SL[i])*usL[nv]  + (S1R - S1L)*usc[nv] +
                  (SR[i] - S1R)*usR[nv] -
                  SR[i]*uR[nv] + SL[i]*uL[nv] + fR[i][nv] - fL[i][nv];

           if (fabs(scrh) > 1.e-6){
             printf (" ! Consistency condition violated, pt %d, nv %d, %12.6e \n", 
                     i,nv,scrh);
             printf (" scrhL = %12.6e   scrhR = %12.6e\n",scrhL, scrhR);
             printf (" SL = %12.6e, S1L = %12.6e, S1R = %12.6e, SR = %12.6e\n",
                     SL[i],S1L,S1R, SR[i]);
             Show(state->vL,i);
             Show(state->vR,i);

             exit(1);
           }
         }
        #endif	  
                 
      }
    }
  }
}
#endif /* end #if on EOS  */

#undef VERIFY_CONSISTENCY_CONDITION 
