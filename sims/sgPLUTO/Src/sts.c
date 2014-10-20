/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Super Time Stepping driver for integration of diffusion terms.

  Take one step in the solution of the diffusion equation 
  \f$ dV/dt = R  \f$ where R is a nonlinear right hand side involving 
  second derivatives. 
  The super step is taken to be equal to the current time step
  ::g_dt and the number of steps Nsts is given by solving Eq. (??)
  of Alexiades et al. with the explicit parabolic time step being computed
  from
  \f[
    \frac{2}{N_d} \max\left[  \frac{\eta_x}{\Delta x^2} 
                            + \frac{\eta_y}{\Delta y^2} 
                            + \frac{\eta_z}{\Delta z^2} \right] 
    \Delta t_{exp} = C_p < \frac{1}{N_d}
  \f]
  where \f$C_p\f$ is the parabolic Courant number, \f$ N_d \f$ is the number
  of spatial dimensions and the maximum of the square
  bracket is computed during the call to ::ParabolicRHS.
  
  This function is called in an operator-split way before/after advection has 
  been carried out.
 
  \b References
     - Alexiades, V., Amiez, A., \& Gremaud E.-A. 1996, 
       Com. Num. Meth. Eng., 12, 31

  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date    Oct 29, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

#ifndef STS_nu 
 #define STS_nu  0.01
#endif

#define STS_MAX_STEPS 1024

static void   STS_ComputeSubSteps(double, double tau[], int);
static double STS_FindRoot(double, double, double);
static double STS_CorrectTimeStep(int, double);

/* ********************************************************************* */
void STS (const Data *d, Time_Step *Dts, Grid *grid)
/*!
 * Solve diffusion equation using Super-Time-Stepping.
 *
 * \param [in,out]  d    pointer to Data structure
 * \param [in,out]  Dts  pointer to Time_Step structure  
 * \param [in]     grid  pointer to an array of Grid structures
 *
 *
 *********************************************************************** */
{
  int i, j, k, nv, n, m;
  static unsigned char *flag;
  double N, ts[STS_MAX_STEPS];
  double dt_par, tau, tsave;
  static double   **v;
  static Data_Arr UU, rhs;

  if (UU == NULL) { 
    UU   = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    rhs  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double); 
    flag = ARRAY_1D(NMAX_POINT, unsigned char);
    v    = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

/* --------------------------------------------------
    time step should not exceed the stability limit.
    Adjust the parabolic CFL number to control it.
   -------------------------------------------------- */

  g_intStage = 1;
  Boundary(d, ALL_DIR, grid);
  #if SHOCK_FLATTENING == MULTID
   FindShock (d, grid);
  #endif
  tsave = g_time;

/* -------------------------------------------------------------------
    obtain the conservative vector UU from the primitive one in order 
    to start cycle. This will be useless if the data structure 
    contains Vc as well as Uc (for future improvements).
   ------------------------------------------------------------------- */

  KDOM_LOOP(k){
  JDOM_LOOP(j){
    IDOM_LOOP(i){
    for (nv = NVAR; nv--;   ){
      v[i][nv] = d->Vc[nv][k][j][i];
    }}
    PrimToCons (v, UU[k][j], IBEG, IEND);
  }}

/* -------------------------------------------------
    Compute the parabolic time step by calling 
    the RHS function once outside the main loop.
   ------------------------------------------------- */

/*  Dts->inv_dtp = STS_TC_RHS(d->Vc, rhs, grid); */
  Dts->inv_dtp = ParabolicRHS(d, rhs, 1.0, grid); 

/* ----------------------------------------------------
    compute (explicit) parabolic time step.
    Restriction on explicit time step should be

      +-                            -+
      | 2*eta_x    2*eta_y   2*eta_z |
      | -------- + ------- + ------- | dt < 1
      |  dx*dx      dy*dy     dz*dz  |
      +-                            -+

    which in PLUTO is expressed through

      +-                      -+
      |   2        2       2   | dt    1
      | ------ + ----- + ----- | -- = ---  = cfl_par
      | dtp_x    dtp_y   dtp_z | Nd    Nd
      +-                      -+

    where Nd is the number of spatial dimensions.
   ---------------------------------------------------- */

  Dts->inv_dtp /= (double) DIMENSIONS;  
  #ifdef PARALLEL
   MPI_Allreduce (&Dts->inv_dtp, &tau, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   Dts->inv_dtp = tau;
  #endif
  dt_par = Dts->cfl_par/(2.0*Dts->inv_dtp); /* -- explicit parabolic time step -- */

/* -------------------------------------------------
    Compute the number of steps needed to fit 
    the supertimestep with the advection one.
   ------------------------------------------------- */

  N = STS_FindRoot(1.0, dt_par, g_dt);
  N = floor(N+1.0);
  n = (int)N;

  Dts->Nsts = n;
  if (n > STS_MAX_STEPS){
    print1 ("! STS: the number of substeps (%d) is > %d\n",n, STS_MAX_STEPS); 
    QUIT_PLUTO(1);
  }

  if (Dts->Nsts > 1){

  /* ---------------------------------------------------------
      Adjust explicit timestep to make Nsts an integer number
     --------------------------------------------------------- */

    dt_par = STS_CorrectTimeStep(Dts->Nsts, g_dt);

  /* -----------------------------------------
             Compute the substeps
     ----------------------------------------- */

    STS_ComputeSubSteps(dt_par, ts, Dts->Nsts);
  }

/* -------------------------------------------
    first step is done outside the main loop
    to save computational time
   ------------------------------------------- */

  if (Dts->Nsts == 1) tau = g_dt;
  else                tau = ts[n-1];

/* VERIFY 
if (n > 1){
  double sum=ts[n-1];
  for (m = 1; m < n; m++){
    sum += ts[n-m-1];
  }
  if (fabs(sum-g_dt) > 1.e-14) {
    print ("PlutoError IN STS\n");
    QUIT_PLUTO(1);
  }
}
*/
  DOM_LOOP(k,j,i){
    #if VISCOSITY == SUPER_TIME_STEPPING
     EXPAND(UU[k][j][i][MX1] += tau*rhs[k][j][i][MX1];  ,
            UU[k][j][i][MX2] += tau*rhs[k][j][i][MX2];  ,
            UU[k][j][i][MX3] += tau*rhs[k][j][i][MX3];)
    #endif
    #if RESISTIVE_MHD == SUPER_TIME_STEPPING
     EXPAND(UU[k][j][i][BX1] += tau*rhs[k][j][i][BX1];  ,
            UU[k][j][i][BX2] += tau*rhs[k][j][i][BX2];  ,
            UU[k][j][i][BX3] += tau*rhs[k][j][i][BX3];)
    #endif
    #if EOS != ISOTHERMAL 
     #if (THERMAL_CONDUCTION == SUPER_TIME_STEPPING) || \
         (RESISTIVE_MHD      == SUPER_TIME_STEPPING) || \
         (VISCOSITY          == SUPER_TIME_STEPPING) 
      UU[k][j][i][ENG] += tau*rhs[k][j][i][ENG];
     #endif
    #endif
  }

  g_dir = IDIR;
  g_i = &i; g_j = &j; g_k = &k;
  KDOM_LOOP(k){
  JDOM_LOOP(j){
    ConsToPrim (UU[k][j], v, IBEG, IEND, flag);
    IDOM_LOOP(i){
    for (nv = NVAR; nv--;   ){
      d->Vc[nv][k][j][i] = v[i][nv];
    }}
  }}   

/* ------------------------------------------------------------
               Main STS Loop starts here
   ------------------------------------------------------------ */

  for (m = 1; m < n; m++){

    g_intStage = m;
    Boundary(d, ALL_DIR, grid);

    tau = ts[n-m-1];

    ParabolicRHS(d, rhs, 1.0, grid); 
    DOM_LOOP (k,j,i){
      #if VISCOSITY == SUPER_TIME_STEPPING
       EXPAND(UU[k][j][i][MX1] += tau*rhs[k][j][i][MX1];  ,
              UU[k][j][i][MX2] += tau*rhs[k][j][i][MX2];  ,
              UU[k][j][i][MX3] += tau*rhs[k][j][i][MX3];)
      #endif
      #if RESISTIVE_MHD == SUPER_TIME_STEPPING
       EXPAND(UU[k][j][i][BX1] += tau*rhs[k][j][i][BX1];  ,
              UU[k][j][i][BX2] += tau*rhs[k][j][i][BX2];  ,
              UU[k][j][i][BX3] += tau*rhs[k][j][i][BX3];)
      #endif
      #if EOS != ISOTHERMAL 
       #if (THERMAL_CONDUCTION == SUPER_TIME_STEPPING) || \
           (RESISTIVE_MHD      == SUPER_TIME_STEPPING) || \
           (VISCOSITY          == SUPER_TIME_STEPPING) 
        UU[k][j][i][ENG] += tau*rhs[k][j][i][ENG]; 
/*  UU[k][j][i][ENG] += tau*rhs[ENG][k][j][i];  */
       #endif
      #endif
    }

  /* ------------------------------------
      convert conservative variables to
      primitive for  next iteration 
     ------------------------------------ */

    g_dir = IDIR;
    g_i = &i; g_j = &j; g_k = &k;
    KDOM_LOOP(k){
    JDOM_LOOP(j){
      ConsToPrim (UU[k][j], v, IBEG, IEND, flag);
      IDOM_LOOP(i){
      for (nv = NVAR; nv--;   ){
        d->Vc[nv][k][j][i] = v[i][nv];
      }}
    }}   
    g_time += ts[n-m-1];
  }
  g_time = tsave;

}

/* ********************************************************************* */
void STS_ComputeSubSteps(double dtex, double tau[], int ssorder)
/*
 *
 *********************************************************************** */
{
  int i;
  double S;

  S = 0.0;

  for (i = 0; i < ssorder; i++) {
    tau[i] = dtex / ((-1.0+STS_nu)*cos(((2.0*i+1.0)*CONST_PI)/(2.0*ssorder)) + 1.0 + STS_nu);
    S += tau[i];
  }
}

/* ********************************************************************* */
double STS_FindRoot(double x0, double dtr, double dta)
/*
 *
 *********************************************************************** */
{
  double a,b,c;
  double Ns, Ns1;
  int n;

  n = 0;

  Ns = x0+1.0;
  Ns1 = x0;
   
  while(fabs(Ns-Ns1) >= 1.0e-5){
    Ns = Ns1;
    a = (1.0-sqrt(STS_nu))/(1.0+sqrt(STS_nu));
    b = pow(a,2.0*Ns);
    c = (1.0-b)/(1.0+b);
    Ns1 = Ns + (dta - dtr*Ns/(2.0*sqrt(STS_nu))*c)/(dtr/(2.0*sqrt(STS_nu))*(c-2.0*Ns*b*log(a)*(1.0+c)/(1.0+b)));
    /*printf("%d\n", n);*/
    n += 1;
  }
  return(Ns);
}

/* ********************************************************************* */
double STS_CorrectTimeStep(int n0, double dta)
/*
 *
 *********************************************************************** */
{
  double a,b,c;
  double dtr;

  a = (1.0-sqrt(STS_nu))/(1.0+sqrt(STS_nu));
  b = pow(a,2.0*n0);
  c = (1.0-b)/(1.0+b);

  dtr = dta*2.0*sqrt(STS_nu)/(n0*c);
  return(dtr);
}

#undef STS_MAX_STEPS 
