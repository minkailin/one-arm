/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Compute viscosity fluxes.      */
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* *********************************************************************** */
void ViscousFlux (Data_Arr V, real **ViF, real **ViS, double **dcoeff, 
                  int beg, int end, Grid *grid)

/*!
 *    \brief Computes viscous fluxes and source terms for the HD/MHD equations. 
 *
 *    \param [in]      V  data array containing cell-centered quantities
 *    \param [in,out]  ViF pointer to viscous fluxes
 *    \param [in,out]  ViS pointer to viscous source terms
 *    \param [in,out]  dcoeff  pointer to diffusion coefficient for dt calculation
 *    \param [in]      beg     integer, index for loop beg
 *    \param [in]      end     integer, index for loop end
 *    \param [in]      grid  pointer to array of Grid structures 
 *
 *    \return This function has no return value.
 *
 *    \notes Calculates the stress tensor components (at the i+1/2 face of each cell) and
 *     adds explicit viscous terms to the energy and momentum equation. It is 
 *     called in the during the sweep integrators. The stress tensor is given by
 *
 *     \f[                                                      
 *     T = \left(
 *     \begin{array}{ccc}
 *                     T_{xx} & T_{xy} & T_{xz} \\
 *                     T_{yx} & T_{yy} & T_{yz} \\
 *                     T_{zx} & T_{zy} & T_{zz} 
 *     \end{array}\right)
 *     \f]
 *    
 *     where \f$ T_{ij} = T_{ji}\f$ and the components are given by 
 *     
 *    \f$ T_{ij} = 2\,m\,e(ij) + (l - 2/3 m) \nabla\cdot {\bf V}\, \delta_{ij} \f$
 *    
 *    where \f$ e(ij)= e_{ij}/(h_ih_j)\f$ and  \f$ e_{ij} = 0.5( V_{i;j} + V_{j;i} ) \f$      
 *    whereas m,l are the first and second parameter of viscosity respectively
 *
 *
 *  \author Petros Tzeferacos (petros.tzeferacos@ph.unito.it) \& Andrea Mignone
 *  \date   Nov 10th, 2010
 ************************************************************************* */
                   
#define D_DX_I(q)  (q[k][j][i + 1] - q[k][j][i])
#define D_DY_J(q)  (q[k][j + 1][i] - q[k][j][i])
#define D_DZ_K(q)  (q[k + 1][j][i] - q[k][j][i])

#define D_DY_I(q)  (  0.25*(q[k][j + 1][i] + q[k][j + 1][i + 1]) \
                    - 0.25*(q[k][j - 1][i] + q[k][j - 1][i + 1]))

#define D_DZ_I(q)  (  0.25*(q[k + 1][j][i] + q[k + 1][j][i + 1])  \
                    - 0.25*(q[k - 1][j][i] + q[k - 1][j][i + 1]))

#define D_DX_J(q)  (  0.25*(q[k][j][i + 1] + q[k][j + 1][i + 1]) \
                    - 0.25*(q[k][j][i - 1] + q[k][j + 1][i - 1]))

#define D_DZ_J(q)  (  0.25*(q[k + 1][j][i] + q[k + 1][j + 1][i]) \
                    - 0.25*(q[k - 1][j][i] + q[k - 1][j + 1][i]))

#define D_DX_K(q)  (  0.25*(q[k][j][i + 1] + q[k + 1][j][i + 1]) \
                    - 0.25*(q[k][j][i - 1] + q[k + 1][j][i - 1]))

#define D_DY_K(q)  (  0.25*(q[k][j + 1][i] + q[k + 1][j + 1][i]) \
                    - 0.25*(q[k][j - 1][i] + q[k + 1][j - 1][i]))                
{
  int i,j,k,n,nv;
  real eta1_visc,eta2_visc,vi[NVAR];
  real div_v;
  real dx_1, dy_1, dz_1;
  real dx1, dx2, dx3;
  real x1, x2, x3;
  real dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz;
  static double *tau1_xx, *tau1_xy, *tau1_xz, *tau1_yx, 
       *tau1_yy, *tau1_yz, *tau1_zx, *tau1_zy,
       *tau1_zz; /* -- Visc Stress Tensor components for sweep Idir-- */
  static double *tau2_xx, *tau2_xy, *tau2_xz, *tau2_yx, 
       *tau2_yy, *tau2_yz, *tau2_zx, *tau2_zy, 
       *tau2_zz; /* -- Visc Stress Tensor components for sweep Jdir-- */
  static double *tau3_xx, *tau3_xy, *tau3_xz, *tau3_yx, 
       *tau3_yy, *tau3_yz, *tau3_zx, *tau3_zy, 
       *tau3_zz; /* -- Visc Stress Tensor components for sweep Kdir-- */
  real ***Vx, ***Vy, ***Vz;
  real scrh, scrh1, scrh2, scrh3;
  real *r, *th, r_1, dr, dr_1, s_1, tan_1, dy;
  real dVxi,dVyi,dVzi;
  real dVxj,dVyj,dVzj;
  real dVxk,dVyk,dVzk;
  static real *one_dVr, *one_dmu; /*auxillary volume components for r_1 singularity @ cylindrical and spherical*/
  
  EXPAND(Vx  = V[VX1]; , Vy  = V[VX2]; , Vz  = V[VX3];)
  
  if (tau1_xx == NULL){
    tau1_xx = ARRAY_1D(NX1_TOT, double);
    tau1_xy = ARRAY_1D(NX1_TOT, double);
    tau1_xz = ARRAY_1D(NX1_TOT, double);
    tau1_yx = ARRAY_1D(NX1_TOT, double);
    tau1_yy = ARRAY_1D(NX1_TOT, double);
    tau1_yz = ARRAY_1D(NX1_TOT, double);
    tau1_zx = ARRAY_1D(NX1_TOT, double);
    tau1_zy = ARRAY_1D(NX1_TOT, double);
    tau1_zz = ARRAY_1D(NX1_TOT, double);
    
    tau2_xx = ARRAY_1D(NX2_TOT, double);
    tau2_xy = ARRAY_1D(NX2_TOT, double);
    tau2_xz = ARRAY_1D(NX2_TOT, double);
    tau2_yx = ARRAY_1D(NX2_TOT, double);
    tau2_yy = ARRAY_1D(NX2_TOT, double);
    tau2_yz = ARRAY_1D(NX2_TOT, double);
    tau2_zx = ARRAY_1D(NX2_TOT, double);
    tau2_zy = ARRAY_1D(NX2_TOT, double);
    tau2_zz = ARRAY_1D(NX2_TOT, double);
          
    tau3_xx = ARRAY_1D(NX3_TOT, double);
    tau3_xy = ARRAY_1D(NX3_TOT, double);
    tau3_xz = ARRAY_1D(NX3_TOT, double);
    tau3_yx = ARRAY_1D(NX3_TOT, double);
    tau3_yy = ARRAY_1D(NX3_TOT, double);
    tau3_yz = ARRAY_1D(NX3_TOT, double);
    tau3_zx = ARRAY_1D(NX3_TOT, double);
    tau3_zy = ARRAY_1D(NX3_TOT, double);
    tau3_zz = ARRAY_1D(NX3_TOT, double);
  }

  r  = grid[IDIR].x; th = grid[JDIR].x;
  if (one_dVr == NULL){
    one_dVr = ARRAY_1D(NX1_TOT, double);  /* -- intercell (i) and (i+1) volume -- */
    one_dmu = ARRAY_1D(NX2_TOT, double);  /* -- intercell (j) and (j+1) volume -- */
    for (i = 0; i < NX1_TOT - 1; i++){
      one_dVr[i] = r[i+1]*fabs(r[i + 1]) - r[i]*fabs(r[i]);
      one_dVr[i] = 2.0/one_dVr[i];
    }
    for (j = 0; j < NX2_TOT - 1; j++){
      one_dmu[j] = 1.0 - cos(th[j + 1]) - (1.0 - cos(th[j]))*(th[j] > 0.0 ? 1.0:-1.0);
      one_dmu[j] = 1.0/one_dmu[j];
    }
  }
                    
  if (g_dir == IDIR){     /* ------- X 1     S W E E P ------------- */
                        /* ------ CALCULATE DERIVATIVES ---------- */

    j = *g_j; k = *g_k;
    dy_1 = 1.0/grid[JDIR].dx[j]; dz_1 = 1.0/grid[KDIR].dx[k];
    
    for (i = beg; i <= end; i++){

      dx_1 = 1.0/grid[IDIR].dx[i];

      tau1_xx[i]= tau1_xy[i]= tau1_xz[i]= tau1_yx[i]= 
      tau1_yy[i]= tau1_yz[i]= tau1_zx[i]= tau1_zy[i]= 
      tau1_zz[i]=0.0;

      dVxi = dVxj = dVxk = dVyi = dVyj = dVyk = dVzi = dVzj = dVzk = 0.0;

      dxVx = dxVy = dxVz = dyVx = dyVy = dyVz = dzVx = dzVy = dzVz = 0.0;

      EXPAND (dVxi = D_DX_I(Vx);, dVyi = D_DX_I(Vy);, dVzi = D_DX_I(Vz);)
      #if DIMENSIONS >= 2
       EXPAND (dVxj = D_DY_I(Vx);, dVyj = D_DY_I(Vy);, dVzj = D_DY_I(Vz);)
       #if DIMENSIONS == 3
        dVxk = D_DZ_I(Vx); dVyk = D_DZ_I(Vy); dVzk = D_DZ_I(Vz);   
       #endif
      #endif

      EXPAND(dxVx = dVxi*dx_1; , dxVy = dVyi*dx_1; , dxVz = dVzi*dx_1; )
      #if DIMENSIONS >= 2 
       EXPAND(dyVx = dVxj*dy_1; , dyVy = dVyj*dy_1; , dyVz = dVzj*dy_1; )
       #if DIMENSIONS == 3
        dzVx = dVxk*dz_1; dzVy = dVyk*dz_1; dzVz = dVzk*dz_1;
       #endif
      #endif

      #if GEOMETRY == CARTESIAN
      
       /* -- calculate the div U -- */
       
       div_v = D_EXPAND(dxVx, + dyVy , + dzVz); 

      /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
         -------------------------------------------------------------------- */
      x1  = grid[IDIR].xr[i]; x2  = grid[JDIR].x[j]; x3  = grid[KDIR].x[k];

      for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j][i + 1]);      
      eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

      /* --------------------------------------------------------------------
         find stress tensor components(not all are needed for flux computation)             
         -------------------------------------------------------------------- */                                

      tau1_xx[i] = 2.0*eta1_visc*dxVx + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v; /*2eta1 dxVx + (eta2 - 2/3 eta1) divV*/
      tau1_xy[i] = eta1_visc*(dyVx + dxVy);   /*eta1 (dyVx + dxVy)*/
      tau1_xz[i] = eta1_visc*(dzVx + dxVz);   /*eta1 (dzVx + dxVz)*/
      tau1_yx[i] = tau1_xy[i];                /*tau1 xy*/
      tau1_zx[i] = tau1_xz[i];                /*tau1 xz*/
       
      /* --------------------------------------------------------------------
                             compute fluxes
         -------------------------------------------------------------------- */

      EXPAND(ViF[i][MX1] = -tau1_xx[i]; ,
             ViF[i][MX2] = -tau1_xy[i]; ,
             ViF[i][MX3] = -tau1_xz[i]; )

      #if EOS != ISOTHERMAL
       ViF[i][ENG] = EXPAND(- 0.5*(Vx[k][j][i] + Vx[k][j][i + 1])*tau1_xx[i],
                           - 0.5*(Vy[k][j][i] + Vy[k][j][i + 1])*tau1_yx[i],
                           - 0.5*(Vz[k][j][i] + Vz[k][j][i + 1])*tau1_zx[i]);
      #endif                                               

      /* --------------------------------------------------------------------
                             compute source terms
         -------------------------------------------------------------------- */

      EXPAND(ViS[i][MX1] = 0.0; ,
             ViS[i][MX2] = 0.0; ,  
             ViS[i][MX3] = 0.0; )                              
                 
     #elif GEOMETRY == CYLINDRICAL

      r = grid[IDIR].x; dr = grid[IDIR].dx[i];
      dr_1 = 1.0/dr; dz_1 = 0.0;
      r_1 = 1.0/grid[IDIR].x[i];
       
       /* -- calculate the div U -- */

      dxVx = (Vx[k][j][i + 1]*r[i + 1] - Vx[k][j][i]*fabs(r[i]))*one_dVr[i];  /* trick @ axis for dxVx*/
      div_v = D_EXPAND(  dxVx, + dVyj*dy_1, + 0.0      ); 
      dxVx = dVxi*dx_1;
      /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
         -------------------------------------------------------------------- */

      x1  = grid[IDIR].xr[i]; x2  = grid[JDIR].x[j]; x3  = grid[KDIR].x[k];
      
      for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j][i + 1]);      
      eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

      /* --------------------------------------------------------------------
         find stress tensor components(not all are needed for flux computation )             
         -------------------------------------------------------------------- */                                

      tau1_xx[i] = 2.0*eta1_visc*dxVx + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v; /* tau_rr = 2 eta1 drVr + (eta2 - 2/3 eta1) divV    */
      tau1_xy[i] = eta1_visc*(dyVx + dxVy);   /*tau_rz = eta1 (dzVr + drVz)*/
      EXPAND(tau1_xz[i] = eta1_visc*(0.0);,
             tau1_xz[i] = eta1_visc*(0.0);,   
             tau1_xz[i] = eta1_visc*0.5*(r[i]+r[i+1])*dr_1*((1./r[i+1])*Vz[k][j][i+1] - (1./r[i])*Vz[k][j][i]);)   
                               /*tau_rphi = eta1 (1/r dphiVr + drVphi - 1/r Vphi)= eta1 (r dr(1/r Vphi) )*/    
      tau1_yx[i] = tau1_xy[i]; /*tau_zr*/
      tau1_zx[i] = tau1_xz[i]; /*tau_phir*/
      
      /* -- we calculate at the center cause we don't need it for flux but for src, avoiding 1/r->inf at r=r_f=0 -- */

      x1  = grid[IDIR].x[i]; x2  = grid[JDIR].x[j]; x3  = grid[KDIR].x[k];

      for (nv = 0; nv < NVAR; nv++) vi[nv] =(V[nv][k][j][i]);      
      eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

      div_v = D_EXPAND( 0.5*(Vx[k][j][i+1]-Vx[k][j][i-1])*dx_1 + Vx[k][j][i]*r_1,
                      + 0.5*(Vy[k][j + 1][i]-Vy[k][j - 1][i])*dy_1, 
                      + 0.0 ); 
 
      tau1_zz[i] = 2.0*eta1_visc*(r_1*(Vx[k][j][i])) + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;
       
       /* --------------------------------------------------------------------
                             compute fluxes
         -------------------------------------------------------------------- */

      EXPAND(ViF[i][MX1] = -tau1_xx[i]; ,
             ViF[i][MX2] = -tau1_xy[i]; ,
             ViF[i][MX3] = -tau1_xz[i]; )
             
      #if EOS != ISOTHERMAL
       ViF[i][ENG] = EXPAND(- 0.5*(Vx[k][j][i] + Vx[k][j][i + 1])*tau1_xx[i],
                           - 0.5*(Vy[k][j][i] + Vy[k][j][i + 1])*tau1_yx[i],
                           - 0.5*(Vz[k][j][i] + Vz[k][j][i + 1])*tau1_zx[i]);
      #endif

      /* --------------------------------------------------------------------
                             compute source terms
         -------------------------------------------------------------------- */
     
      EXPAND(ViS[i][MX1] = -(tau1_zz[i])*r_1; ,
             ViS[i][MX2] = 0.0;                ,  
             ViS[i][MX3] = 0.0; ) 

     #elif GEOMETRY == POLAR 
      
      r = grid[IDIR].xr; th = grid[JDIR].x;
      dr   = grid[IDIR].dx[i]; dr_1 = 1.0/dr;
      r_1  = 1.0/grid[IDIR].xr[i];

      /* -- calculate the div U -- */
      div_v =  D_EXPAND(r_1*0.5*(Vx[k][j][i] + Vx[k][j][i + 1]) + dVxi*dr_1 ,
                       + r_1*dVyj*dy_1                                      , 
                       + dVzk*dz_1                                          ); 


      /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
         -------------------------------------------------------------------- */

      x1  = grid[IDIR].xr[i]; x2  = grid[JDIR].x[j]; x3  = grid[KDIR].x[k];

      for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j][i + 1]);      
      eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

      /* --------------------------------------------------------------------
         find stress tensor components(not all are needed for flux computation ) 
         -------------------------------------------------------------------- */                                

      tau1_xx[i] = 2.0*eta1_visc*dxVx + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v; /*tau_rr as in cylindrical */
      EXPAND(tau1_xy[i] = eta1_visc*(dxVy );,
             tau1_xy[i] = eta1_visc*(r_1*dyVx + dxVy - r_1*0.5*(Vy[k][j][i] + Vy[k][j][i + 1]));,
             tau1_xy[i] = eta1_visc*(r_1*dyVx + dxVy - r_1*0.5*(Vy[k][j][i] + Vy[k][j][i + 1]));)
                      /*tau_rphi = eta1 (1/r dphiVr + drVphi - 1/r Vphi) */
      tau1_xz[i] = eta1_visc*(dzVx + dxVz);  /*tau_rz as in cylindrical*/
      tau1_yx[i] = tau1_xy[i];               /*tau_phir = tau_rphi*/
      EXPAND(tau1_yy[i] = 2.0*eta1_visc*(r_1*dyVy )+ (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;,
             tau1_yy[i] = 2.0*eta1_visc*(r_1*dyVy  
                        + r_1*0.5*(Vx[k][j][i] + Vx[k][j][i + 1]))+ (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;,
             tau1_yy[i] = 2.0*eta1_visc*(r_1*dyVy 
                        + r_1*0.5*(Vx[k][j][i] + Vx[k][j][i + 1]))+ (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;)
                        /*tau_phiphi = 2 eta1 (1/r dphiVphi + 1/r Vr) + (eta2 - 2/3 eta1) divV*/
      tau1_zx[i] = tau1_xz[i]; /*tau_zr = tau_rz*/
      /* --------------------------------------------------------------------
                             compute fluxes
         -------------------------------------------------------------------- */

      EXPAND(ViF[i][MX1] = -tau1_xx[i]; ,
             ViF[i][MX2] = -tau1_xy[i]; ,
             ViF[i][MX3] = -tau1_xz[i]; )

      #if EOS != ISOTHERMAL
       ViF[i][ENG] = EXPAND(- 0.5*(Vx[k][j][i] + Vx[k][j][i + 1])*tau1_xx[i],
                           - 0.5*(Vy[k][j][i] + Vy[k][j][i + 1])*tau1_yx[i],
                           - 0.5*(Vz[k][j][i] + Vz[k][j][i + 1])*tau1_zx[i]);
      #endif

      /* --------------------------------------------------------------------
                             compute source terms
         -------------------------------------------------------------------- */
      r_1  = 1.0/grid[IDIR].x[i];
      
      EXPAND(ViS[i][MX1] = -0.5*(tau1_yy[i - 1] + tau1_yy[i])*r_1; ,
             ViS[i][MX2] = 0.0;                                    ,
             ViS[i][MX3] = 0.0;                                    )
                
                                            
     #elif GEOMETRY == SPHERICAL

      r = grid[IDIR].xr; th = grid[JDIR].x;
      dr_1 = 1.0/grid[IDIR].dx[i]; r_1  = 1.0/grid[IDIR].xr[i];
      tan_1= 1.0/tan(th[j]); s_1 = 1.0/sin(th[j]);

      /* -- calculate the div U -- */
      div_v = D_EXPAND(  2.0*r_1*0.5*(Vx[k][j][i] + Vx[k][j][i + 1]) + dVxi*dr_1 ,
                      + r_1*dVyj*dy_1 + r_1*tan_1*0.5*(Vy[k][j][i] + Vy[k][j][i + 1]), 
                      + r_1*s_1*dVzk*dz_1                                        ); 
       

      /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
         -------------------------------------------------------------------- */
      x1  = grid[IDIR].xr[i]; x2  = grid[JDIR].x[j]; x3  = grid[KDIR].x[k];

      for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j][i + 1]);      
      eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

      /* --------------------------------------------------------------------
         find stress tensor components(not all are needed for flux computation )       
         -------------------------------------------------------------------- */                                

      tau1_xx[i] = 2.0*eta1_visc*dxVx + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v; /*tau_rr = 2eta1 drVr + + (eta2 - 2/3 eta1) divV*/
      EXPAND(tau1_xy[i] = eta1_visc*(r_1*dyVx + dxVy );,
             tau1_xy[i] = eta1_visc*(r_1*dyVx + dxVy - 0.5*r_1*(Vy[k][j][i] + Vy[k][j][i + 1]));,
             tau1_xy[i] = eta1_visc*(r_1*dyVx + dxVy - 0.5*r_1*(Vy[k][j][i] + Vy[k][j][i + 1]));)
                      /*tau_rthita = eta1 (1/r dthitaVr + drVthita -1/r Vthita)*/
      EXPAND(tau1_xz[i] = eta1_visc*(r_1*s_1*dzVx + dxVz) ;,  
             tau1_xz[i] = eta1_visc*(r_1*s_1*dzVx + dxVz );,  
             tau1_xz[i] = eta1_visc*(r_1*s_1*dzVx + dxVz - 0.5*r_1*(Vz[k][j][i] + Vz[k][j][i + 1]));)
                      /*tau_rphi = eta1 (1/r dthitaVr + drVthita -1/r Vthita)*/
             tau1_yx[i] = tau1_xy[i]; /*tau_thitar*/
             tau1_yy[i] = 2.0*eta1_visc*(r_1*dyVy + 0.5*r_1*(Vx[k][j][i] + Vx[k][j][i + 1])) 
                        + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v; 
                        /*tau_thitathita= 2 eta1 (1/r dthitaVthita + 1/r Vr) + (eta2 - 2/3 eta1)divV*/
      tau1_zx[i] = tau1_xz[i]; /*tau_phir*/
      EXPAND(tau1_zz[i] = 2.0*eta1_visc*(r_1*s_1*dzVz + 0.5*r_1*(Vx[k][j][i] + Vx[k][j][i + 1]) ) 
                        + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;,
             tau1_zz[i] = 2.0*eta1_visc*(r_1*s_1*dzVz + 0.5*r_1*(Vx[k][j][i] + Vx[k][j][i + 1]) + tan_1*r_1*0.5*(Vy[k][j][i] + Vy[k][j][i + 1])) 
                        + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;,
             tau1_zz[i] = 2.0*eta1_visc*(r_1*s_1*dzVz + 0.5*r_1*(Vx[k][j][i] + Vx[k][j][i + 1]) + tan_1*r_1*0.5*(Vy[k][j][i] + Vy[k][j][i + 1])) 
                        + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;)
                        /*tau_phiphi= 2 eta1 (1/rs dphiVphi + 1/r Vr + 1/r cot Vthita) + (eta2 - 2/3 eta1)divV*/
      /* --------------------------------------------------------------------
                             compute fluxes
         -------------------------------------------------------------------- */

      EXPAND(ViF[i][MX1] = -tau1_xx[i]; ,
             ViF[i][MX2] = -tau1_xy[i]; ,
             ViF[i][MX3] = -tau1_xz[i]; )

      #if EOS != ISOTHERMAL
       ViF[i][ENG] = EXPAND(- 0.5*(Vx[k][j][i] + Vx[k][j][i + 1])*tau1_xx[i],
                           - 0.5*(Vy[k][j][i] + Vy[k][j][i + 1])*tau1_yx[i],
                           - 0.5*(Vz[k][j][i] + Vz[k][j][i + 1])*tau1_zx[i]);
      #endif

      /* --------------------------------------------------------------------
                              compute source terms
         -------------------------------------------------------------------- */

      r_1  = 1.0/grid[IDIR].x[i];
      EXPAND(ViS[i][MX1] = -0.5*(tau1_yy[i - 1] + tau1_yy[i] + tau1_zz[i - 1] + tau1_zz[i])*r_1;,
             ViS[i][MX2] = 0.0;,
             ViS[i][MX3] = 0.0;)       
               
     #endif  /* -- end #if GEOMETRY -- */
/*
     *nu_max = MAX(*nu_max,eta1_visc/vi[RHO]);
     *nu_max = MAX(*nu_max,eta2_visc/vi[RHO]);
*/
      dcoeff[i][MX1]  = MAX(eta1_visc, eta2_visc);
      dcoeff[i][MX1] /= vi[RHO];
    }

  }else if (g_dir == JDIR){    /* ------- X 2     S W E E P  ------------- */

    i = *g_i; k = *g_k;
    dx_1 = 1.0/grid[IDIR].dx[i]; dz_1 = 1.0/grid[KDIR].dx[k];

    for (j = beg ; j <= end; j++){
      dy_1 = 1.0/grid[JDIR].dx[j];
      
      tau2_xx[j]= tau2_xy[j]= tau2_xz[j]= tau2_yx[j]= tau2_yy[j]= 
      tau2_yz[j]= tau2_zx[j]= tau2_zy[j]= tau2_zz[j]=0.0;
      
      dVxi = dVxj = dVxk = dVyi = dVyj = dVyk = dVzi = dVzj = dVzk = 0.0;
      dxVx = dxVy = dxVz = dyVx = dyVy = dyVz = dzVx = dzVy = dzVz = 0.0;

    /* ------ CALCULATE DERIVATIVES ---------- */

      EXPAND (dVxi = D_DX_J(Vx);, dVyi = D_DX_J(Vy);, dVzi = D_DX_J(Vz);)
      EXPAND (dVxj = D_DY_J(Vx);, dVyj = D_DY_J(Vy);, dVzj = D_DY_J(Vz);)
      #if DIMENSIONS == 3
       dVxk = D_DZ_J(Vx); dVyk = D_DZ_J(Vy); dVzk = D_DZ_J(Vz);   
      #endif
      
      EXPAND(dxVx = dVxi*dx_1; , dxVy = dVyi*dx_1; , dxVz = dVzi*dx_1; )
      EXPAND(dyVx = dVxj*dy_1; , dyVy = dVyj*dy_1; , dyVz = dVzj*dy_1; )
      #if DIMENSIONS == 3
       dzVx = dVxk*dz_1; dzVy = dVyk*dz_1; dzVz = dVzk*dz_1;
      #endif
     
      #if GEOMETRY == CARTESIAN
   
       /* -- calculate the div U -- */
       
       div_v = D_EXPAND(  dxVx  , + dyVy  , + dzVz  ); 
       
      /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
         -------------------------------------------------------------------- */

       x1  = grid[IDIR].x[i]; x2  = grid[JDIR].xr[j]; x3  = grid[KDIR].x[k];

       for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j + 1][i]);      
       eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

      /* --------------------------------------------------------------------
         find stress tensor components(not all are needed for flux computation )                
         -------------------------------------------------------------------- */                                

       tau2_xy[j] = eta1_visc*(dyVx + dxVy);
       tau2_yx[j] = tau2_xy[j];
       tau2_yy[j] = 2.0*eta1_visc*dyVy + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;
       tau2_yz[j] = eta1_visc*(dyVz + dzVy);
       tau2_zy[j] = tau2_yz[j];
       
      /* --------------------------------------------------------------------
                             compute fluxes
         -------------------------------------------------------------------- */

       EXPAND(ViF[j][MX1] = -tau2_yx[j]; ,
              ViF[j][MX2] = -tau2_yy[j]; ,
              ViF[j][MX3] = -tau2_yz[j]; )

       #if EOS != ISOTHERMAL
        ViF[j][ENG] = EXPAND(- 0.5*(Vx[k][j][i] + Vx[k][j + 1][i])*tau2_xy[j],
                            - 0.5*(Vy[k][j][i] + Vy[k][j + 1][i])*tau2_yy[j],
                            - 0.5*(Vz[k][j][i] + Vz[k][j + 1][i])*tau2_zy[j]);
       #endif

       /* --------------------------------------------------------------------
                              compute source terms
          -------------------------------------------------------------------- */

       EXPAND(ViS[j][MX1] = 0.0; ,
              ViS[j][MX2] = 0.0; ,  
              ViS[j][MX3] = 0.0; )                              

      #elif GEOMETRY == CYLINDRICAL
       
       r = grid[IDIR].x; th = grid[KDIR].x; dy_1 = 1.0/grid[JDIR].dx[j];
       dr = grid[IDIR].dx[i]; dr_1 = 1.0/dr;
       r_1 = 1.0/grid[IDIR].x[i]; dz_1 = 0.0; 
       /* -- calculate the div U -- */
       
       div_v = D_EXPAND(  r_1*0.5*(Vx[k][j][i] + Vx[k][j + 1][i]) + dVxi*dr_1,
                        + dVyj*dy_1                                          , 
                        + 0.0                                                ); 

       /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
          -------------------------------------------------------------------- */

       x1  = grid[IDIR].x[i]; x2  = grid[JDIR].xr[j]; x3  = grid[KDIR].x[k];

       for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j + 1][i]);      
       eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

       /* --------------------------------------------------------------------
          find stress tensor components(not all are needed for flux computation )               
          -------------------------------------------------------------------- */                                
      
       tau2_xy[j] = eta1_visc*(dyVx + dxVy); /*tau_rz = eta1 (dzVr + drVz)*/
       tau2_yx[j] = tau2_xy[j];  /*tau_zr */
       tau2_yy[j] = 2.0*eta1_visc*dyVy + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v; /*tau_zz = 2eta1 dzVz + (eta2 - 2/3 eta1)divV*/
       EXPAND(tau2_yz[j] = 0.0;, tau2_yz[j] = 0.0;, tau2_yz[j] = eta1_visc*(dyVz);) /*tau_zphi = eta1(1/r dphiVz + dzVphi)*/
       tau2_zy[j] = tau2_yz[j]; /*tau_phiz*/
       
      /* --------------------------------------------------------------------
                             compute fluxes
         -------------------------------------------------------------------- */

       EXPAND(ViF[j][MX1] = -tau2_yx[j]; ,
              ViF[j][MX2] = -tau2_yy[j]; ,
              ViF[j][MX3] = -tau2_yz[j]; )

       #if EOS != ISOTHERMAL
        ViF[j][ENG] = EXPAND(- 0.5*(Vx[k][j][i] + Vx[k][j + 1][i])*tau2_xy[j],
                            - 0.5*(Vy[k][j][i] + Vy[k][j + 1][i])*tau2_yy[j],
                            - 0.5*(Vz[k][j][i] + Vz[k][j + 1][i])*tau2_zy[j]);
       #endif

       /* --------------------------------------------------------------------
                              compute source terms
          -------------------------------------------------------------------- */

       EXPAND(ViS[j][MX1] = 0.0; ,
              ViS[j][MX2] = 0.0; ,  
              ViS[j][MX3] = 0.0; )                              

      #elif GEOMETRY == POLAR 

       r = grid[IDIR].x; th = grid[JDIR].xr; dr   = grid[IDIR].dx[i];
       dr_1 = 1.0/dr; r_1  = 1.0/grid[IDIR].x[i]; dy_1 = 1.0/grid[JDIR].dx[j];
       
       /* -- calculate the div U -- */
       
       div_v = D_EXPAND(  r_1*0.5*(Vx[k][j][i] + Vx[k][j + 1][i]) + dVxi*dr_1  ,
                        + r_1*dVyj*dy_1                , 
                        + dVzk*dz_1                    ); 

       /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
          -------------------------------------------------------------------- */

       x1  = grid[IDIR].x[i]; x2  = grid[JDIR].xr[j]; x3  = grid[KDIR].x[k];
      
       for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j + 1][i]);      
       eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

       /* --------------------------------------------------------------------
          find stress tensor components(not all are needed for flux computation )      
          -------------------------------------------------------------------- */                                

       tau2_xy[j] = eta1_visc*(r_1*dyVx + dxVy - r_1*0.5*(Vy[k][j][i] + Vy[k][j + 1][i])); 
                    /*tau_rphi = eta1(1/r dphiVr + drVphi -1/r Vphi)*/
       tau2_yx[j] = tau2_xy[j]; /*tau_phir*/
       tau2_yy[j] = 2.0*eta1_visc*(r_1*dyVy 
                  + r_1*0.5*(Vx[k][j][i] + Vx[k][j + 1][i]))+ (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;
                  /*tau_phiphi = 2 eta1 (1/r dphiVphi + 1/rVr) + (eta2 - 2/3 eta1)divV*/
       tau2_yz[j] = eta1_visc*(r_1*dyVz + dzVy); /* tau_phiz = eta1(dzVphi + 1/r dphiVz)*/
       tau2_zy[j] = tau2_yz[j]; /*tau_zphi*/
     
       /* --------------------------------------------------------------------
                              compute fluxes
          -------------------------------------------------------------------- */

       EXPAND(ViF[j][MX1] = -tau2_yx[j]; ,
              ViF[j][MX2] = -tau2_yy[j]; ,
              ViF[j][MX3] = -tau2_yz[j]; )

       #if EOS != ISOTHERMAL
        ViF[j][ENG] = EXPAND(- 0.5*(Vx[k][j][i] + Vx[k][j + 1][i])*tau2_xy[j],
                            - 0.5*(Vy[k][j][i] + Vy[k][j + 1][i])*tau2_yy[j],
                            - 0.5*(Vz[k][j][i] + Vz[k][j + 1][i])*tau2_zy[j]);
       #endif   
       r_1  = 1.0/grid[IDIR].x[i];

       /* --------------------------------------------------------------------
                              compute source terms
          -------------------------------------------------------------------- */
       
       EXPAND(ViS[j][MX1] = 0.0; ,
              ViS[j][MX2] = 0.0; ,
              ViS[j][MX3] = 0.0; )
      
      #elif GEOMETRY == SPHERICAL
       
       r = grid[IDIR].x; th = grid[JDIR].xr; dy_1 = 1.0/grid[JDIR].dx[j];
       dr_1 = 1.0/grid[IDIR].dx[i]; r_1  = 1.0/grid[IDIR].x[i];
       tan_1= 1.0/tan(th[j]); s_1 = 1.0/sin(th[j]);
       /*------------------------------------
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         Hotfix for tan_1 at the axis: 
        since we have terms tan_1*smth that 
        cannot be treated with the volume
        trick, we set an if condition 
        for this to go to zero, as it is
        the correct behaviour of the  
        term in question there.
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       --------------------------------------*/
       if (fabs(tan(th[j]))< 1.e-12) tan_1 = 0.0;      
       
       /* -- calculate the div U -- */

       th = grid[JDIR].x;
       
       div_v = D_EXPAND( 2.0*r_1*0.5*(Vx[k][j][i] + Vx[k][j + 1][i]) + dVxi*dr_1    ,
                       + ( sin(th[j + 1])*Vy[k][j + 1][i] - fabs(sin(th[j]))*Vy[k][j][i])*r_1*one_dmu[j], 
                       + r_1*s_1*dVzk*dz_1 ); 
      
       th = grid[JDIR].xr;
        /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
         -------------------------------------------------------------------- */

       x1  = grid[IDIR].x[i]; x2  = grid[JDIR].xr[j]; x3  = grid[KDIR].x[k];

      for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j + 1][i]);      
      eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

      /* --------------------------------------------------------------------
         find stress tensor components(not all are needed for flux computation )                 
         -------------------------------------------------------------------- */                                
      tau2_xy[j] = eta1_visc*(r_1*dyVx + dxVy - 0.5*r_1*(Vy[k][j][i] + Vy[k][j + 1][i])); 
                  /*tau_rthita = eta1(1/r dthitaVr + drVthita -1/r Vthita)*/
      tau2_yx[j] = tau2_xy[j]; /*tau_thitar*/
      tau2_yy[j] = 2.0*eta1_visc*(r_1*dyVy + 0.5*r_1*(Vx[k][j][i] + Vx[k][j + 1][i])) 
                 + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;
                 /*tau_thitathita= 2 eta1 (1/r dthitaVthita + 1/r Vr) + (eta2 - 2/3 eta1)divV*/
      #if DIMENSIONS <= 2                  
       EXPAND(tau2_yz[j] = eta1_visc*(r_1*dyVz);,
              tau2_yz[j] = eta1_visc*(r_1*dyVz);,
              tau2_yz[j] = eta1_visc*(r_1*dyVz - 0.5*tan_1*r_1*(Vz[k][j][i] + Vz[k][j + 1][i])));
      #endif       
      #if DIMENSIONS == 3                  
       tau2_yz[j] = eta1_visc*(s_1*r_1*dzVy + r_1*dyVz - 0.5*tan_1*r_1*(Vz[k][j][i] + Vz[k][j + 1][i]));
      #endif       /*tau_thitaphi= eta1 (1/rs dphiVthita + 1/r dthitaVphi -1/r cot Vphi)*/
      tau2_zy[j] = tau2_yz[j]; /*tau_phithita*/
      #if DIMENSIONS <= 2                  
       EXPAND(tau2_zz[j] = 2.0*eta1_visc*(0.5*r_1*(Vx[k][j][i] + Vx[k][j + 1][i])) 
                         + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;,
              tau2_zz[j] = 2.0*eta1_visc*(0.5*r_1*(Vx[k][j][i] + Vx[k][j + 1][i]) 
                         + tan_1*r_1*0.5*(Vy[k][j][i] + Vy[k][j + 1][i])) 
                         + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;,
              tau2_zz[j] = 2.0*eta1_visc*(0.5*r_1*(Vx[k][j][i] + Vx[k][j + 1][i]) 
                         + tan_1*r_1*0.5*(Vy[k][j][i] + Vy[k][j + 1][i])) 
                         + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;)
      #endif       
      #if DIMENSIONS == 3                  
       tau2_zz[j] = 2.0*eta1_visc*(r_1*s_1*dzVz + 0.5*r_1*(Vx[k][j][i] + Vx[k][j + 1][i]) 
                  + tan_1*r_1*0.5*(Vy[k][j][i] + Vy[k][j + 1][i])) 
                  + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;
      #endif       /*tau_phiphi = 2eta1(1/rs dphiVphi + 1/r Vr + 1/r cot Vthita) + (eta2 -2/3 eta1)divV */
                  
      /* --------------------------------------------------------------------
                             compute fluxes
         -------------------------------------------------------------------- */

      EXPAND(ViF[j][MX1] = -tau2_yx[j]; ,
             ViF[j][MX2] = -tau2_yy[j]; ,
             ViF[j][MX3] = -tau2_yz[j]; )

      #if EOS != ISOTHERMAL
       ViF[j][ENG] = EXPAND(- 0.5*(Vx[k][j][i] + Vx[k][j + 1][i])*tau2_xy[j],
                           - 0.5*(Vy[k][j][i] + Vy[k][j + 1][i])*tau2_yy[j],
                           - 0.5*(Vz[k][j][i] + Vz[k][j + 1][i])*tau2_zy[j]);
      #endif

      /* --------------------------------------------------------------------
                             compute source terms
         -------------------------------------------------------------------- */

      th = grid[JDIR].x; tan_1= 1.0/tan(th[j]);
      
      EXPAND(ViS[j][MX1] = 0.0;                                                                                 ,
             ViS[j][MX2] = 0.5*(tau2_yx[j - 1] + tau2_yx[j])*r_1 - tan_1*0.5*(tau2_zz[j - 1] + tau2_zz[j])*r_1; ,
             ViS[j][MX3] = 0.0; )

     #endif  /* -- end #if GEOMETRY -- */
      
/*
     *nu_max = MAX(*nu_max,eta1_visc/vi[RHO]);
     *nu_max = MAX(*nu_max,eta2_visc/vi[RHO]);
*/
     dcoeff[j][MX1]  = MAX(eta1_visc,eta2_visc);
     dcoeff[j][MX1] /= vi[RHO];
    }      

  }else if (g_dir == KDIR){       /* ------- X 3     S W E E P  ------------- */
   
    i = *g_i; j = *g_j;
   
    dx_1 = 1.0/grid[IDIR].dx[i]; dy_1 = 1.0/grid[JDIR].dx[j];

    for (k = beg; k <= end; k++){
      dz_1 = 1.0/grid[KDIR].dx[k];

      tau3_xx[k]= tau3_xy[k]= tau3_xz[k]= tau3_yx[k]= tau3_yy[k]= tau3_yz[k]= 
      tau3_zx[k]= tau3_zy[k]= tau3_zz[k]=0.0;
      dVxi = dVxj = dVxk = dVyi = dVyj = dVyk = dVzi = dVzj = dVzk = 0.0;

      dVxi = D_DX_K(Vx); dVyi = D_DX_K(Vy); dVzi = D_DX_K(Vz);
      dVxj = D_DY_K(Vx); dVyj = D_DY_K(Vy); dVzj = D_DY_K(Vz);
      dVxk = D_DZ_K(Vx); dVyk = D_DZ_K(Vy); dVzk = D_DZ_K(Vz);
  
      /* ------ CALCULATE DERIVATIVES ---------- */

      dxVx = dVxi*dx_1; dxVy = dVyi*dx_1; dxVz = dVzi*dx_1; 
      dyVx = dVxj*dy_1; dyVy = dVyj*dy_1; dyVz = dVzj*dy_1;
      dzVx = dVxk*dz_1; dzVy = dVyk*dz_1; dzVz = dVzk*dz_1;
         
      #if GEOMETRY == CARTESIAN
        
       /* -- calculate the div U -- */
       
       div_v = dxVx + dyVy + dzVz;
      
       /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
          -------------------------------------------------------------------- */

       x1  = grid[IDIR].x[i]; x2  = grid[JDIR].x[j]; x3  = grid[KDIR].xr[k];
       for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k + 1][j][i]);      
       eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

       /* --------------------------------------------------------------------
          find stress tensor components(not all are needed for flux computation )
          -------------------------------------------------------------------- */                                

       tau3_xz[k] = eta1_visc*(dzVx + dxVz);  
       tau3_yz[k] = eta1_visc*(dyVz + dzVy);
       tau3_zx[k] = tau3_xz[k];
       tau3_zy[k] = tau3_yz[k];
       tau3_zz[k] = 2.0*eta1_visc*dzVz + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;

       /* --------------------------------------------------------------------
                             compute fluxes
          -------------------------------------------------------------------- */

       ViF[k][MX1] = -tau3_zx[k];
       ViF[k][MX2] = -tau3_zy[k]; 
       ViF[k][MX3] = -tau3_zz[k]; 

       #if EOS != ISOTHERMAL
        ViF[k][ENG] = - 0.5*(Vx[k][j][i] + Vx[k + 1][j][i])*tau3_xz[k]
                     - 0.5*(Vy[k][j][i] + Vy[k + 1][j][i])*tau3_yz[k]
                     - 0.5*(Vz[k][j][i] + Vz[k + 1][j][i])*tau3_zz[k];
       #endif

       /* --------------------------------------------------------------------
                              compute source terms
          -------------------------------------------------------------------- */

       EXPAND(ViS[k][MX1] = 0.0; ,
              ViS[k][MX2] = 0.0; ,  
              ViS[k][MX3] = 0.0; )                              

      #elif GEOMETRY == POLAR 

       r = grid[IDIR].x; th = grid[JDIR].x;
       dr   = grid[IDIR].dx[i]; dr_1 = 1.0/dr;
       r_1  = 1.0/grid[IDIR].x[i]; dy_1 = 1.0/grid[JDIR].dx[j];

       /* -- calculate the div U -- */
       
       div_v = r_1*0.5*(Vx[k][j][i] + Vx[k + 1][j][i]) + dVxi*dr_1
             + r_1*dVyj*dy_1
             + dVzk*dz_1; 

       /* --------------------------------------------------------------------
                           compute viscosity (face centered) 
          -------------------------------------------------------------------- */

       x1  = grid[IDIR].x[i]; x2  = grid[JDIR].x[j]; x3  = grid[KDIR].xr[k];
      
       for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k + 1][j][i]);      
       eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

       /* --------------------------------------------------------------------
          find stress tensor components(not all are needed for flux computation )         
          -------------------------------------------------------------------- */                                

       tau3_xz[k] = eta1_visc*(dzVx + dxVz);  /*tau_rz = eta1(drVz + dzVr)*/
       tau3_yy[k] = 2.0*eta1_visc*(r_1*dyVy + r_1*0.5*(Vx[k][j][i] + Vx[k + 1][j][i]))
                  + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;
                  /*tau_phiphi = 2eta1(1/r dphiVphi + 1/r V) + (eta2- 2/3 eta1)divV*/
       tau3_yz[k] = eta1_visc*(r_1*dyVz + dzVy); /*tau_phiz = eta1(dzVphi + 1/r dphiVz)*/
       tau3_zx[k] = tau3_xz[k]; /*tau_zr*/
       tau3_zy[k] = tau3_yz[k]; /*tau_zphi*/
       tau3_zz[k] = 2.0*eta1_visc*dzVz + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;
                  /*tau_zz = 2eta1(dzVz) + (eta2- 2/3 eta1)divV*/

       /* --------------------------------------------------------------------
                              compute fluxes
          -------------------------------------------------------------------- */

       ViF[k][MX1] = -tau3_zx[k]; 
       ViF[k][MX2] = -tau3_zy[k]; 
       ViF[k][MX3] = -tau3_zz[k]; 

       #if EOS != ISOTHERMAL
        ViF[k][ENG] = - 0.5*(Vx[k][j][i] + Vx[k + 1][j][i])*tau3_xz[k]
                     - 0.5*(Vy[k][j][i] + Vy[k + 1][j][i])*tau3_yz[k]
                     - 0.5*(Vz[k][j][i] + Vz[k + 1][j][i])*tau3_zz[k];
       #endif

       /* --------------------------------------------------------------------
                              compute source terms
          -------------------------------------------------------------------- */

       EXPAND(ViS[k][MX1] = 0.0; ,
              ViS[k][MX2] = 0.0; ,  
              ViS[k][MX3] = 0.0; )                              
                   
      #elif GEOMETRY == SPHERICAL
       r = grid[IDIR].x; th = grid[JDIR].x;
       dy_1 = 1.0/grid[JDIR].dx[j]; dr_1 = 1.0/grid[IDIR].dx[i];
       r_1  = 1.0/grid[IDIR].x[i]; tan_1= 1.0/tan(th[j]); s_1 = 1.0/sin(th[j]);

       /* -- calculate the div U -- */

       div_v = 2.0*r_1*0.5*(Vx[k][j][i] + Vx[k + 1][j][i]) + dVxi*dr_1 
             + r_1*dVyj*dy_1 + r_1*tan_1*0.5*(Vy[k][j][i] + Vy[k][j + 1][i])
             + r_1*s_1*dVzk*dz_1; 
       /* --------------------------------------------------------------------
                            compute viscosity (face centered) 
          -------------------------------------------------------------------- */

       x1  = grid[IDIR].x[i]; x2  = grid[JDIR].x[j]; x3  = grid[KDIR].xr[k];
      
       for (nv = 0; nv < NVAR; nv++) vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k + 1][j][i]);      
       eta_visc_func(vi, x1, x2, x3, &eta1_visc, &eta2_visc);

       /* --------------------------------------------------------------------
          find stress tensor components(not all are needed for flux computation )               
          -------------------------------------------------------------------- */                                

       tau3_xz[k] = eta1_visc*(r_1*s_1*dzVx + dxVz - 0.5*r_1*(Vz[k][j][i] + Vz[k + 1][j][i]));  
                    /*tau_rphi = eta1(drVphi + 1/rs dphiVr - 1/r Vphi)*/ 
       tau3_yz[k] = eta1_visc*(s_1*r_1*dzVy + r_1*dyVz - 0.5*tan_1*r_1*(Vz[k][j][i] + Vz[k + 1][j][i]));
                    /*tau_thitaphi = eta1(1/rs dphiVthita + 1/r dthitaVphi - 1/r cot Vphi)*/ 
       tau3_zx[k] = tau3_xz[k]; /* tau_phir */
       tau3_zy[k] = tau3_yz[k]; /* tau_phithita */
       tau3_zz[k] = 2.0*eta1_visc*(r_1*s_1*dzVz + 0.5*r_1*(Vx[k][j][i] + Vx[k + 1][j][i]) + tan_1*r_1*0.5*(Vy[k][j][i] + Vy[k + 1][j][i])) 
                  + (eta2_visc - (2.0/3.0)*eta1_visc)*div_v;
                    /*tau_phiphi = 2eta1(1/rs dphiVphi +1/r Vr +1/r cot Vthita) + (eta2 - 2/3 eta1)divV*/
       /* --------------------------------------------------------------------
                              compute fluxes
          -------------------------------------------------------------------- */

       ViF[k][MX1] = -tau3_zx[k]; 
       ViF[k][MX2] = -tau3_zy[k]; 
       ViF[k][MX3] = -tau3_zz[k]; 

       #if EOS != ISOTHERMAL
        ViF[k][ENG] = - 0.5*(Vx[k][j][i] + Vx[k + 1][j][i])*tau3_xz[k]
                     - 0.5*(Vy[k][j][i] + Vy[k + 1][j][i])*tau3_yz[k]
                     - 0.5*(Vz[k][j][i] + Vz[k + 1][j][i])*tau3_zz[k];
       #endif

       /* --------------------------------------------------------------------
                              compute source terms
          -------------------------------------------------------------------- */

       EXPAND(ViS[k][MX1] = 0.0; ,
              ViS[k][MX2] = 0.0; ,  
              ViS[k][MX3] = 0.0; )                              
                                     
      #endif  /* -- end #if GEOMETRY -- */
/*
      *nu_max = MAX(*nu_max,eta1_visc/vi[RHO]);
      *nu_max = MAX(*nu_max,eta2_visc/vi[RHO]);
*/
      dcoeff[k][MX1]  = MAX(eta1_visc, eta2_visc);
      dcoeff[k][MX1] /= vi[RHO];
    }/*loop*/
  }/*sweep*/

/* --------------------------------------------------
    now reduce the current (inverse) time step 
    for this sweep only
   -------------------------------------------------- */

/*  return REDUCE_LOCAL_DT (nu, 2, beg, end, grid); */

}/*function*/
