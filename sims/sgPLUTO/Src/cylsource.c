#include"pluto.h"

static void CONS_CYLSOURCE   (double *u, double *);
static double ROS34_SOLVE (real *v0, real *k1, real *v4th, real dt, real tol);

static void CYLSOURCE   (double *u, double *);
static void CYLJACOBIAN (double *v, double **);

static double rad;

/* ***************************************************** */
void CONS_CYLSOLVE (Data_Arr VV, double dt, Grid *grid)
/*
 *
 *
 ******************************************************* */
{
  int i,j,k,nv;
  double u0[NVAR], u1[NVAR];
  double k1[NVAR], k2[NVAR], k3[NVAR], k4[NVAR];
  static double **u, **v;
  static unsigned char *flag;

  #if GEOMETRY != CARTESIAN
   return;
  #endif

  if (u == NULL){
    u = ARRAY_2D(2*NMAX_POINT,NVAR, double);
    v = ARRAY_2D(2*NMAX_POINT,NVAR, double);
    flag = ARRAY_1D(2*NMAX_POINT, unsigned char);
  }

  KTOT_LOOP(k){ 
  JTOT_LOOP(j){

  /* --------------------------------------
      convert to conservative vars
     -------------------------------------- */

    ITOT_LOOP(i){
      for (nv = NVAR; nv--;  ) v[i][nv] = VV[nv][k][j][i];
    }
    PrimToCons(v, u, 0, NX1_TOT-1);

  /* --------------------------------------
       Integrate
     -------------------------------------- */

    ITOT_LOOP(i){      
      rad = grid[IDIR].x[i];
      for (nv = 0; nv < NVAR; nv++) u0[nv] = u[i][nv];
      CONS_CYLSOURCE(u0, k1);
      for (nv = NVAR; nv--; ) u1[nv] = u0[nv] + 0.5*dt*k1[nv]; 

      CONS_CYLSOURCE(u1, k2);
      for (nv = NVAR; nv--; ){
        u[i][nv] = u0[nv] + dt*k2[nv];
      }

if (u[i][RHO] < 0.0){
  printf ("! CYLSOURCE: negative density\n");
  QUIT_PLUTO(1);
}
if (u[i][ENG] < 0.0){
  printf ("! CYLSOURCE: negative energy\n");
  QUIT_PLUTO(1);
}

      continue;


      for (nv = 0; nv < NVAR; nv++) u0[nv] = u[i][nv];
  
      /* -- step 1 -- */

      CONS_CYLSOURCE(u0, k1);
      for (nv = NVAR; nv--; ) u1[nv] = u0[nv] + 0.5*dt*k1[nv]; 

      /* -- step 2 -- */

      CONS_CYLSOURCE(u1, k2);
      for (nv = NVAR; nv--; ) u1[nv] = u0[nv] + 0.5*dt*k2[nv]; 

      /* -- step 3 -- */

      CONS_CYLSOURCE(u1, k3);
      for (nv = NVAR; nv--; ) u1[nv] = u0[nv] + dt*k3[nv]; 

      /* -- step 4 -- */

      CONS_CYLSOURCE(u1, k4);

      for (nv = 0; nv < NVAR; nv++) {
        u[i][nv] = u0[nv] + dt*(k1[nv] + 2.0*(k2[nv] + k3[nv]) + k4[nv])/6.0;
      }
    }

  /* --------------------------------------
      convert to primitive vars
     -------------------------------------- */

    ConsToPrim(u,v, 0, NX1_TOT-1, flag);
    ITOT_LOOP(i){
      for (nv = NVAR; nv--;  )  VV[nv][k][j][i] = v[i][nv];
    }
  }}

}
/* *********************************************************** */
void CONS_CYLSOURCE(double *u, double *S)
/*
 *
 * Provide the righ hand side of the system of ODE
 * in conservative variables, i.e.
 *
 *  S = -F/r
 *
 *
 ************************************************************* */
{
  int nv;
  double rho, tau, p;
  double vr, vz, vphi;
  double Br=0.0, Bz=0.0, Bphi=0.0, Kin, pm=0.0, vB=0.0;

  rho = u[RHO];
  tau = 1.0/rho;
  vr   = u[MX1]*tau;
  vz   = u[MX2]*tau;
  vphi = u[MX3]*tau;

  #if PHYSICS == MHD
   Br   = u[BX1];
   Bz   = u[BX2];
   Bphi = u[BX3];
  #endif

  pm  = 0.5*(Br*Br + Bz*Bz + Bphi*Bphi);
  Kin = 0.5*(vr*vr + vz*vz + vphi*vphi)*rho;
  p   = (g_gamma - 1.0)*(u[ENG] - Kin - pm);
  vB  = vr*Br + vz*Bz + vphi*Bphi;

  S[RHO] = rho*vr;  
  S[MX1] = rho*(vr*vr - vphi*vphi) - Br*Br + Bphi*Bphi;
  S[MX2] = rho*vr*vz - Br*Bz;
  S[MX3] = 2.0*(rho*vr*vphi - Bphi*Br);
  S[ENG] = vr*(u[ENG] + p + pm) - Br*vB;
  #if PHYSICS == MHD
   S[BX1] = 0.0;
   S[BX2] = vr*Bz - vz*Br;
   S[BX3] = 0.0;
  #endif

  #if NVAR != NFLX 
   for (nv = NFLX; nv < NVAR; nv++) S[nv] = u[nv]*vr;
  #endif

  for (nv = 0; nv < NVAR; nv++) S[nv] *= -1.0/rad;
}










/* ***************************************************** */
void CYLSOLVE (const Data *d, double dt, Grid *grid)
/*
 *
 * Integrate the system of ODE
 *
 *    dU/dt = S = -F/r
 *
 * arising when writing the axisymmetric MHD equations 
 * using a Cartesian+Source term formalism.
 * The F's are the fluxes (except for pressure, which 
 * must be discretized as gradient).
 * This routine must be called prior to an integrator
 * (in CARTESIAN coordinates):
 *
 * -->  CYLSOLVE (&d, g_dt, grid);
 *
 *
 * To make the integration process simpler, we switch
 * to primitive variables, using
 *
 *  dU/dV * dV/dt = S  -->  dV/dt = dV/dU*S 
 * 
 * The following MAPLE script does the algebra:
 *
--------------------------------------------------------
restart;
with(linalg);
K  := v[r]^2 + v[y]^2 + v[phi]^2;
pm := B[r]^2 + B[y]^2 + B[phi]^2;
vB := v[r]*B[r] + v[y]*B[y] + v[phi]*B[phi];

E := p/G1 + rho*K/2 + pm/2;

U := vector( [rho, v[r]*rho, v[y]*rho, v[phi]*rho, 
                   B[r], B[y], B[phi], E, psi, rho*T]);
V := vector( [rho, v[r], v[y], v[phi], 
                   B[r], B[y], B[phi], p, psi,T]);

dUdV := jacobian(U, V);
dVdU := inverse(dUdV);

S := vector( [rho*v[r], 
      rho*v[r]*v[r] - B[r]*B[r] + B[phi]^2 - rho*v[phi]^2, 
      rho*v[y]*v[r] - B[r]*B[y], 
      2*(rho*v[phi]*v[r] - B[r]*B[phi]),
      0, v[r]*B[y] - v[y]*B[r], 0, 
      (E+p+pm)*v[r] - B[r]*vB, 
       rho*T*v[r]]);

SV := multiply(dVdU,S);
dSdV := simplify(jacobian(SV, V));

simplify(dSdV[8,5]);
--------------------------------------------------------

 *
 *
 *
 ******************************************************* */
{
  int i,j,k,nv;
  double v0[NVAR], v1[NVAR], S[NVAR];

  #if GEOMETRY != CARTESIAN
   return;
  #endif

  KTOT_LOOP(k){ 
  JTOT_LOOP(j){
    ITOT_LOOP(i){
      rad = grid[IDIR].x[i];
      for (nv = 0; nv < NVAR; nv++) {
        v0[nv] = d->Vc[nv][k][j][i];
      }
      CYLSOURCE(v0, S);
      ROS34_SOLVE(v0, S, v1, dt, 1.e-5);
      for (nv = 0; nv < NVAR; nv++) {
        d->Vc[nv][k][j][i] = v1[nv];
      }
    }
  }}
}

/* *********************************************************** */
void CYLSOURCE(double *v, double *S)
/*
 *
 *
 * Provide the righ hand side of the system of ODE
 * in primitive variables, i.e.
 *
 *  S = -dV/dU * F/r
 *
 *
 ************************************************************* */
{
  int nv;
  double g1, vB, pm;

  g1 = g_gamma - 1.0;

  EXPAND(S[RHO] =  v[RHO]*v[VX1];  ,
                                ,
         S[VX1] = -v[VX3]*v[VX3];)
  EXPAND(S[VX2] =  0.0;           ,
                                 ,
         S[VX3] =  v[VX3]*v[VX1];)
  S[PRS] = g_gamma*v[PRS]*v[VX1];

  #if PHYSICS == MHD
   S[VX1] -= (EXPAND(v[BX1]*v[BX1],    , - v[BX3]*v[BX3]))/v[RHO];
   S[VX2] -=  v[BX1]*v[BX2];
   #if COMPONENTS == 3
    S[VX3] -= 2.0*v[BX1]*v[BX3]/v[RHO];
   #endif

   EXPAND(S[BX1] = 0.0;                         ,
          S[BX2] = 0.0;                         ,
          S[BX3] = v[VX1]*v[BX2] - v[BX1]*v[VX2];)

   vB     = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
   pm     = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
   S[PRS] += g1*(vB*v[BX1] + 0.5*v[VX1]*pm);
  #endif

  #if NVAR != NFLX 
   for (nv = NFLX; nv < NVAR; nv++) S[nv] = 0.0;
  #endif

  for (nv = 0; nv < NVAR; nv++) S[nv] *= -1.0/rad;

}

/* ************************************************************ */
void CYLJACOBIAN (double *v, double **dSdV)
/*
 *
 *
 * Provide the Jacobian dS/dV
 *
 *
 ************************************************************* */
{
  int ni, nj;
  double g1, tau, tau2;
  double vB, pm;

  g1 = g_gamma - 1.0;

  for (ni = 0; ni < NVAR; ni++){
  for (nj = 0; nj < NVAR; nj++){
    dSdV[ni][nj] = 0.0;
  }}

  tau  = 1.0/v[RHO];
  tau2 = tau*tau;

  EXPAND(dSdV[RHO][RHO] = v[VX1];
         dSdV[RHO][VX1] = v[RHO];    ,
                                  ,
         dSdV[VX1][VX3] =  -2.0*v[VX3];)
  #if PHYSICS == MHD
   EXPAND(dSdV[VX1][RHO] = (v[BX1]*v[BX1] - v[BX3]*v[BX3])*tau2;
          dSdV[VX1][BX1] =  -2.0*v[BX1]*tau;                   ,
                                                            ,
          dSdV[VX1][BX3] =   2.0*v[BX3]*tau;)
 
  #endif

  #if PHYSICS == MHD
   dSdV[VX2][RHO] =  v[BX1]*v[BX2]*tau2;
   dSdV[VX2][BX1] = -v[BX2]*tau;
   dSdV[VX2][BX2] = -v[BX1]*tau;
  #endif

  #if COMPONENTS == 3
   dSdV[VX3][VX1] =  v[VX3];
   dSdV[VX3][VX3] =  v[VX1];
   #if PHYSICS == MHD
    dSdV[VX3][RHO] = 2.0*v[BX1]*v[BX3]*tau2;
    dSdV[VX3][BX1] = -2.0*v[BX3]*tau;
    dSdV[VX3][BX3] = -2.0*v[BX1]*tau;
   #endif
  #endif
  
  #if PHYSICS == MHD
   dSdV[BX2][VX1] =  v[BX2];
   dSdV[BX2][VX2] = -v[BX1];
   dSdV[BX2][VX3] =  0.0;
   dSdV[BX2][BX1] = -v[VX2];
   dSdV[BX2][BX2] =  v[VX1];
  #endif

  dSdV[PRS][VX1] = g_gamma*v[PRS];
  dSdV[PRS][PRS] = g_gamma*v[VX1];

  #if PHYSICS == MHD
   vB = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
   pm = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
   EXPAND(dSdV[PRS][VX1] += 0.5*g1*(pm + 2.0*v[BX1]*v[BX1]);
          dSdV[PRS][VX2] += g1*v[BX1]*v[BX2];                 ,
                                                          ,
          dSdV[PRS][VX3] += g1*v[BX1]*v[BX3];)
   EXPAND(dSdV[PRS][BX1] += g1*(vB + 2.0*v[VX1]*v[BX1]);      ,
          dSdV[PRS][BX2] += g1*(v[VX1]*v[BX2] + v[BX1]*v[VX2]); , 
          dSdV[PRS][BX3] += g1*(v[VX1]*v[BX3] + v[BX1]*v[VX3]);)
  #endif

/* ---------------------------------
         divide by radius 
   --------------------------------- */

  for (ni = 0; ni < NVAR; ni++){
  for (nj = 0; nj < NVAR; nj++){
    dSdV[ni][nj] *= -1.0/rad;
  }}
}



#define GAM (1.0/2.0)
#define A21 2.0
#define A31 (48.0/25.0)
#define A32 (6.0/25.0)
#define C21 -8.0
#define C31 (372.0/25.0)
#define C32 (12.0/5.0)
#define C41 (-112.0/125.0)
#define C42 (-54.0/125.0)
#define C43 (-2.0/5.0)
#define B1 (19.0/9.0)
#define B2 (1.0/2.0)
#define B3 (25.0/108.0)
#define B4 (125.0/108.0)
#define E1 (17.0/54.0)
#define E2 (7.0/36.0)
#define E3 0.0
#define E4 (125.0/108.0)

#define C1X (1.0/2.0)
#define C2X (-3.0/2.0)
#define C3X (121.0/50.0)
#define C4X (29.0/250.0)
#define A2X 1.0
#define A3X (3.0/5.0)

/* ****************************************************************************************** */
double ROS34_SOLVE (real *v0, real *k1, real *v4th, real dt, real tol)
/*
 * 
 *   Solve using Semi-Implicit Rosenbrock Method
 *
 ********************************************************************************************* */
{
  int    j, nv, ksub, kfail;
  static int *indx;
  real   err, scrh, t, tend;
  real   v1[NVAR], vscal[NVAR], k2[NVAR], k3[NVAR];
  real   tsub[4096], dt_grow, dt_shrink;
  static real  **a, **J, **J2;
  static real  *g1, *g2, *g3, *g4;

real vbeg[NVAR];

/* -------------------------------------------
    copy ALL variables so that density is
    defined when calling Radiat.
   ------------------------------------------- */
                                                                                                                             
  for (nv = 0; nv < NVAR; nv++) v1[nv] = vbeg[nv] = v0[nv];
  
  if (indx == NULL){
    indx  = ARRAY_1D (NVAR, int);
    a     = ARRAY_2D (NVAR, NVAR, double);
    J     = ARRAY_2D (NVAR, NVAR, double);
    J2    = ARRAY_2D (NVAR, NVAR, double);
    g1    = ARRAY_1D (NVAR, double);
    g2    = ARRAY_1D (NVAR, double);
    g3    = ARRAY_1D (NVAR, double);
    g4    = ARRAY_1D (NVAR, double);
  }

  for (nv = 0; nv < NVAR; nv++) vscal[nv] = fabs(v0[nv]);
  vscal[VX1] = vscal[VX2] = vscal[VX3] = 10.0;
  #if PHYSICS == MHD
   vscal[BX1] = vscal[BX2] = vscal[BX3] = 1.0;
  #endif   
  #if NVAR != NFLX
   for (nv = NFLX; nv < NVAR; nv++) vscal[nv] = 1.0;
  #endif

  CYLJACOBIAN (v0, J);   

  t    = 0.0;
  tend = dt;

  ksub = kfail = 0;
  for (;;){

  /* --  Compute (I - hJ)  -- */

    for (nv = 0; nv < NVAR; nv++) {   
      for (j = 0; j < NVAR; j++) a[nv][j] = -J[nv][j];
      a[nv][nv] += 1.0/(GAM*dt);
    }
    LUDecomp (a, NVAR, indx, &scrh);    /*    LU decomposition of the matrix. */

  /* -- set right hand side for g1 -- */

    for (nv = NVAR; nv--; ) g1[nv] = k1[nv];

  /* -- solve for g1 -- */

    LUBackSubst (a, NVAR, indx, g1); 
    for (nv = NVAR; nv--; ) v1[nv] = v0[nv] + A21*g1[nv];

    CYLSOURCE (v1, k2);    

  /* -- set right hand side for g2 -- */

    for (nv = NVAR; nv--; ) g2[nv] = k2[nv] + C21*g1[nv]/dt;

  /* -- solve for g2 -- */

    LUBackSubst (a, NVAR, indx, g2);  
    for (nv = NVAR; nv--;  ) v1[nv] = v0[nv] + A31*g1[nv] + A32*g2[nv];

    CYLSOURCE(v1, k3);
    
  /* -- set right hand side for g3 -- */

    for (nv = NVAR; nv--; ) 
      g3[nv] = k3[nv] + (C31*g1[nv] + C32*g2[nv])/dt;

  /* -- solve for g3 -- */

    LUBackSubst (a, NVAR, indx, g3);  

  /* -- set right hand side for g4 -- */

    for (nv = NVAR; nv--;   ) { 
      g4[nv] = k3[nv] + (C41*g1[nv] + C42*g2[nv] + C43*g3[nv])/dt;
    }

  /* -- solve for g4  -- */ 

    LUBackSubst (a, NVAR, indx, g4);    

  /* --  4th order solution & error estimation -- */

    err = 0.0;
    for (nv = NVAR; nv--; ) {   
      v4th[nv] = v0[nv] + B1*g1[nv] + B2*g2[nv] 
                        + B3*g3[nv] + B4*g4[nv];
      scrh     = E1*g1[nv] + E2*g2[nv] + E3*g3[nv] + E4*g4[nv];
      err      = MAX(err, fabs(scrh)/vscal[nv]);

/*
if (fabs(scrh)/vscal[nv] > 1.e-4){
printf ("nv = %d, v = %12.6e,  rhs = %12.6e,  err = %12.6e\n",
         nv, v4th[nv], k1[nv], fabs(scrh)/vscal[nv]);
}
*/

    }

   
    err /= tol;
    if (err < 1.0) {

      ksub++;      
      err = MAX(0.0, 1.e-18);

      t          += dt;
      tsub[ksub]  = t;

    /* -- provide an estimate for next dt -- */

      dt_grow = 0.9*dt*pow(err, -0.25);
      dt_grow = MIN(dt_grow, 5.0*dt);
      dt      = dt_grow;

    /* -- exit loop if final time has been reached -- */

      if (fabs(t/tend - 1.0) < 1.e-9) break;

    /* -- adjust dt not to exceed tend -- */

      if (dt > (tend - t)) dt = tend - t;

    /* -- initialize solution vector continuing -- */
   
      for (nv = NVAR; nv--;   ) v0[nv] = v4th[nv];
      CYLSOURCE   (v0, k1);
      CYLJACOBIAN (v0, J);   

      if (ksub > 1000){
        print (" ! COOL_SOLVER_ROS34: Number of substeps too large (%d)\n",ksub);
        QUIT_PLUTO(1);
      }
                                                                                                                             
    }else{   /* -- shrink dt and redo time step -- */

      kfail++;
      dt_shrink = 0.9*dt*pow(err, -1.0/3.0);
      dt        = MAX(dt_shrink, 0.05*dt);                                                                                                                        
    }
  }

  if (ksub > 2) {
    int i;
    print (" ! COOL_SOLVER_CK45: number of substeps is %d\n", ksub);
/*
    for (i = 1; i <= ksub; i++){
      printf ("# %d, dt = %12.6e, t/dt = %f,  tend = %12.6e\n",
           i, tsub[i]-tsub[i-1],tsub[i]/tend, tend);
    }
    printf ("kfail = %d\n",kfail);
    QUIT_PLUTO(1);
*/
  }

  return (dt);
}

#undef GAM 
#undef A21 
#undef A31 
#undef A32 
#undef C21 
#undef C31 
#undef C32 
#undef C41 
#undef C42 
#undef C43 
#undef B1 
#undef B2 
#undef B3 
#undef B4 
#undef E1 
#undef E2 
#undef E3 
#undef E4 
#undef C1X 
#undef C2X 
#undef C3X 
#undef C4X 
#undef A2X 
#undef A3X 


/*
restart;
with(linalg);
unprotect(Gamma);
K  := v[r]^2 + v[y]^2 + v[phi]^2;
pm := B[r]^2 + B[y]^2 + B[phi]^2;
vB := v[r]*B[r] + v[y]*B[y] + v[phi]*B[phi];

E := p/(Gamma - 1) + rho*K/2 + pm/2;

U := vector( [rho, v[r]*rho, v[y]*rho, v[phi]*rho, 
                   B[r], B[y], B[phi], E, rho*T]);
V := vector( [rho, v[r], v[y], v[phi], 
                   B[r], B[y], B[phi], p, T]);

dUdV := jacobian(U, V);
dVdU := inverse(dUdV);

S := vector( [rho*v[r], 
      rho*v[r]*v[r] - B[r]*B[r] + B[phi]^2 - rho*v[phi]^2, 
      rho*v[y]*v[r] - B[r]*B[y], 
      2*(rho*v[phi]*v[r] - B[r]*B[phi]),
      0, v[r]*B[y] - v[y]*B[r], 0, 
      (E+p+pm)*v[r] - B[r]*vB, 
       rho*T*v[r]]);

SV := multiply(dVdU,S);
dSdV := simplify(jacobian(SV, V));


#
#  HD
#
restart;
with(linalg);
unprotect(Gamma);
U := vector( [rho, m[r], m[y], m[phi], E]);

v[r]   := m[r]/rho:
v[phi] := m[phi]/rho:
v[y]   := m[y]/rho:

K      := m[r]^2 + m[phi]^2 + m[y]^2;

p := (Gamma-1)*(E - K/(2*rho));

S := vector( [rho*v[r], 
      rho*v[r]*v[r] - rho*v[phi]^2, 
      rho*v[y]*v[r], 
      2*rho*v[phi]*v[r],
      (E+p)*v[r]]);
dSdU := simplify(jacobian(S, U));


#
#  MHD
#

restart;
with(linalg);
unprotect(Gamma);
U := vector( [rho, 
              m[r], m[y], m[phi], 
              B[r], B[y], B[phi],  
              E, 
              T]);

v[r]   := m[r]/rho:
v[phi] := m[phi]/rho:
v[y]   := m[y]/rho:

vB     := v[r]*B[r] + v[y]*B[y] + v[phi]*B[phi];
pm     := Fm(B[r], B[phi], B[y]); # B[r]*B[r] + B[phi]*B[phi] + B[y]*B[y];
K      := FK(m[r], m[phi], m[y]); # m[r]^2 + m[phi]^2 + m[y]^2;

p := (Gamma-1)*(E - K/(2*rho) - pm/2);

S := vector( [rho*v[r], 
      rho*v[r]*v[r] - B[r]*B[r] + B[phi]^2 - rho*v[phi]^2, 
      rho*v[y]*v[r] - B[r]*B[y], 
      2*(rho*v[phi]*v[r] - B[r]*B[phi]),
      0, v[r]*B[y] - v[y]*B[r], 0, 
      FE(rho, m[r], m[phi], m[y], B[r], B[phi], B[y]), #(E+p+pm)*v[r] - B[r]*vB, 
       rho*T*v[r]]);
dSdU := simplify(jacobian(S, U));
#inverse(dUdV);

*/
