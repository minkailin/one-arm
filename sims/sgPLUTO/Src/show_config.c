/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Print useful information about the current computations.

  Detailed description goes here.         

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
static void CheckConfig();

/* ********************************************************************* */
void ShowConfig ()
/*!
 *  Write a summary of the selected options
 *
      ___  __   __  ____________ 
     / _ \/ /  / / / /_  __/ __ \
    / ___/ /__/ /_/ / / / / /_/ /
   / /  /____/\____/ /_/  \____/ 
   ============================== 
                                  
 *
 * 
 *                        
 *
 *********************************************************************** */
{
  FILE *fconf;
  time_t time_now;
  char  str1[128], str2[128], str3[128];

  CheckConfig();

  print1 ("\n");
  print1("   ___  __   __  ____________   \n");
  print1("  / _ \\/ /  / / / /_  __/ __ \\ \n");
  print1(" / ___/ /__/ /_/ / / / / /_/ /  \n");
  print1("/ /  /____/\\____/ /_/  \\____/   \n");
  print1("=============================    v. %s  \n", PLUTO_VERSION);
  print1 ("\n> System:\n\n");

  if ( (fconf = fopen("sysconf.out","r")) != NULL){

    while (fscanf (fconf, "%s %s %s\n", str1, str2, str3) != EOF) {

      if (!strcmp(str1,"USER")) 
        print1 ("  %s:             %s\n",str1, str3);
      else if (!strcmp(str1,"WORKING_g_dir"))
        print1 ("  %s:      %s\n",str1, str3);
      else if (!strcmp(str1,"SYSTEM_NAME"))
        print1 ("  %s:      %s\n",str1, str3);
      else if (!strcmp(str1,"NODE_NAME"))
        print1 ("  %s:        %s\n",str1, str3);
      else if (!strcmp(str1,"ARCH"))
        print1 ("  %s:             %s\n",str1, str3);
      else if (!strcmp(str1,"BYTE_ORDER"))
        print1 ("  %s:       %s\n\n",str1, str3);
    }

  }else{
    print1 ("! sysconf.out file not found \n\n");
  }
 
  time(&time_now);
  print1("> Local time:    %s\n",asctime(localtime(&time_now)));
      
  if (COMPONENTS < DIMENSIONS) {
    print1 ("Sorry, but the number of vector components can\n");
    print1 ("not be less than the dimension number.\n");
    print1 ("Please edit definitions.h and fix this.\n");
    QUIT_PLUTO(1);
  }

  print1 ("> Configuration:\n\n");
  if (PHYSICS == HD)   print1 ("  PHYSICS:          HD\n");
  if (PHYSICS == RHD)  print1 ("  PHYSICS:          RHD\n");
  if (PHYSICS == MHD)  print1 ("  PHYSICS:          MHD [div.B: ");
  if (PHYSICS == RMHD) print1 ("  PHYSICS:          RMHD [div.B: ");

  #if PHYSICS == MHD || PHYSICS == RMHD
   #if MHD_FORMULATION == NONE
    print1 ("None]\n");
   #elif MHD_FORMULATION == EIGHT_WAVES
    print1 ("Powell's 8wave]\n");
   #elif MHD_FORMULATION == DIV_CLEANING
    print1 ("Divergence Cleaning]\n");
   #elif MHD_FORMULATION == CONSTRAINED_TRANSPORT
     print1 ("CT/");
    #if CT_EMF_AVERAGE == ARITHMETIC
     print1 ("Ar. average]\n");
    #elif CT_EMF_AVERAGE == UCT_CONTACT
     print1 ("UCT_CONTACT]\n");
    #elif CT_EMF_AVERAGE == UCT0
     print1 ("UCT0]\n");
    #elif CT_EMF_AVERAGE == UCT_HLL
     print1 ("UCT_HLL]\n");
    #endif
   #endif
  #endif
  print1 ("  DIMENSIONS:       %d\n", DIMENSIONS);
  print1 ("  COMPONENTS:       %d\n", COMPONENTS);
  print1 ("  TRACERS:          %d\n", NTRACER);
  print1 ("  VARIABLES:        %d\n", NVAR);

  print1 ("  LOADED MODULES:\n");
  #if PHYSICS == MHD
   #ifdef SHEARINGBOX
    print1 ("\n  o [SHEARINGBOX]\n");
    print1 ("     - Order:             %d\n", SB_ORDER);
    print1 ("     - Sym Hydro Flux:    %s\n", 
             (SB_SYMMETRIZE_HYDRO == YES ? "YES":"NO"));
    print1 ("     - Sym Ey:            %s\n", 
             (SB_SYMMETRIZE_EY == YES ? "YES":"NO"));
    print1 ("     - Sym Ez:            %s\n", 
             (SB_SYMMETRIZE_EZ == YES ? "YES":"NO"));
    print1 ("     - Force EMF periods: %s\n", 
             (SB_FORCE_EMF_PERIODS == YES ? "YES":"NO"));
   #endif
  #endif
  #ifdef FARGO
   print1 ("\n  o [FARGO]\n");
   print1 ("     - Order:         %d\n", FARGO_ORDER);
   print1 ("     - Average Speed: %s\n", 
            (FARGO_AVERAGE_VELOCITY == YES ? "YES":"NO"));
   print1 ("     - Av. Frequency: %d\n", FARGO_NSTEP_AVERAGE);

  #endif
  print1 ("\n");

  print1 ("  GEOMETRY:         ");
  if (GEOMETRY == CARTESIAN)    print1 ("Cartesian\n");
  if (GEOMETRY == CYLINDRICAL)  print1 ("Cylindrical\n");
  if (GEOMETRY == POLAR)        print1 ("Polar\n");
  if (GEOMETRY == SPHERICAL)    print1 ("Spherical\n");

  print1 ("  BODY_FORCE:       ");
  print1 (BODY_FORCE == NO ? "NO\n":"EXPLICIT\n");

  print1 ("  ROTATION:         ");
  print1(ROTATING_FRAME == YES ? "YES\n":"NO\n");

  print1 ("  EOS:              ");
  if (EOS == IDEAL)        print1 ("Ideal\n");
  if (EOS == BAROTROPIC)   print1 ("Barotropic\n");
  if (EOS == ISOTHERMAL)   print1 ("Isothermal\n");
  if (EOS == TAUB)         print1 ("Taub - TM\n");

  print1 ("  TIME INTEGRATOR:  ");
  if (TIME_STEPPING == EULER)            print1 ("Euler\n");
  if (TIME_STEPPING == RK2)              print1 ("Runga-Kutta II\n");
  if (TIME_STEPPING == RK3)              print1 ("Runga_Kutta III\n");
  if (TIME_STEPPING == CHARACTERISTIC_TRACING)
                                         print1 ("Characteristic Tracing\n");
  if (TIME_STEPPING == HANCOCK)          print1 ("Hancock\n");

  print1 ("  DIM. SPLITTING:   ");
  if (DIMENSIONAL_SPLITTING == YES)  print1 ("Yes\n");
  else                               print1 ("No\n");
  
  print1 ("  INTERPOLATION:    ");
  #ifndef FINITE_DIFFERENCE
   if (INTERPOLATION == FLAT)          print1 ("Flat");
   if (INTERPOLATION == LINEAR)        print1 ("Linear TVD");
   if (INTERPOLATION == LINEAR_MULTID) print1 ("Linear_Multid");
   if (INTERPOLATION == LimO3)         print1 ("LimO3");
   if (INTERPOLATION == WENO3)         print1 ("WENO 3rd order");
   if (INTERPOLATION == PARABOLIC)     print1 ("Parabolic");
   #ifdef CHAR_LIMITING
    if (CHAR_LIMITING == YES) print1 (" (Characteristic lim)\n");
    else                      print1 (" (Primitive lim)\n");
   #endif
  #endif

  #ifdef FINITE_DIFFERENCE
   if (INTERPOLATION == LIMO3_FD)     print1 ("LimO3 (FD), 3rd order\n");
   if (INTERPOLATION == WENO3_FD)     print1 ("WENO3 (FD), 3rd order\n");
   if (INTERPOLATION == WENOZ_FD)     print1 ("WENOZ (FD) 5th order\n");
   if (INTERPOLATION == MP5_FD)       print1 ("MP5 (FD), 5th order\n");
  #endif

  #if PARABOLIC_FLUX != NO
   print1 ("  DIFFUSION TERMS:");
   #if (RESISTIVE_MHD == EXPLICIT) 
    print1 ("  Resistivity  [EXPLICIT]\n");  
   #elif (RESISTIVE_MHD == SUPER_TIME_STEPPING)
    print1 ("  Resistivity  [STS]\n");  
   #endif

   #if (THERMAL_CONDUCTION == EXPLICIT) 
    print1 ("  Thermal Conduction [EXPLICIT]\n");  
   #elif (THERMAL_CONDUCTION == SUPER_TIME_STEPPING)
    print1 ("  Thermal Conduction [STS]\n");  
   #endif

   #if (VISCOSITY == EXPLICIT) 
    print1 ("  Viscosity [EXPLICIT]\n");  
   #elif (VISCOSITY == SUPER_TIME_STEPPING)
    print1 ("  Viscosity [STS]\n");  
   #endif
  #endif

  print1 ("\n");
}

/* ********************************************************************* */
void ShowUnits ()
/*!
 *  Show units when cooling is enabled.
 *
 *
 *********************************************************************** */
{

  #if COOLING != NO
   print1 ("> Cooling Module:    ");
   if (COOLING == SNEq)  print1 (" SNEq\n");
   if (COOLING == MINEq) print1 (" MINEq\n");
   if (COOLING == TABULATED) print1 (" TABULATED\n");
  #endif

  print1 ("> Normalization Units:\n\n");
  print1 (" [Density]:      %8.3e (gr/cm^3), %8.3e (1/cm^3)\n",
          g_unitDensity,g_unitDensity/CONST_mp);
  print1 (" [Pressure]:     %8.3e (dyne/cm^2)\n",
          g_unitDensity*g_unitVelocity*g_unitVelocity);
  print1 (" [Velocity]:     %8.3e (cm/s)\n",g_unitVelocity);
  print1 (" [Length]:       %8.3e (cm)\n",g_unitLength);
  print1 (" [Temperature]:  %8.3e X (p/rho*mu) (K)\n",KELVIN);
  print1 (" [Time]:         %8.3e (sec), %8.3e (yrs) \n",
       g_unitLength/g_unitVelocity, g_unitLength/g_unitVelocity/86400./365.);
  #if PHYSICS == MHD || PHYSICS == RMHD
   print1 (" [Mag Field]:    %8.3e (Gauss)\n",
            g_unitVelocity*sqrt(4.0*CONST_PI*g_unitDensity));
  #endif

  print1 (" \n");
    
}

/* ********************************************************************* */
void ShowDomainDecomposition(int nprocs, Grid *GXYZ)
/*!
 * Show the parallel domain decomposition by having each processor print
 * its own computational domain.
 * This is activated with the -show-dec command line argument.
 * It may be long for several thousand processors. 
 *
 * \param [in] nprocs   the total number of processors
 * \param [in] GXYZ     a pointer to an array of grid structures
 *********************************************************************** */
{

  int i, j, k, ngh;
  int i0, i1, j0, j1, k0, k1;
  int nxp, nyp, nzp;
  int nghx, nghy, nghz;

  static int *i0_proc, *i1_proc;
  static int *j0_proc, *j1_proc;
  static int *k0_proc, *k1_proc;

  double x0, x1, y0, y1, z0, z1;

  double *x0_proc, *x1_proc;
  double *y0_proc, *y1_proc;
  double *z0_proc, *z1_proc;

  Grid *Gx, *Gy, *Gz;
  
  Gx = GXYZ;
  Gy = GXYZ + 1;
  Gz = GXYZ + 2;

/* ---- Allocate memory ---- */

  i0_proc = ARRAY_1D(nprocs, int); i1_proc = ARRAY_1D(nprocs, int);
  j0_proc = ARRAY_1D(nprocs, int); j1_proc = ARRAY_1D(nprocs, int);
  k0_proc = ARRAY_1D(nprocs, int); k1_proc = ARRAY_1D(nprocs, int);
  x0_proc = ARRAY_1D(nprocs, double); x1_proc = ARRAY_1D(nprocs, double);
  y0_proc = ARRAY_1D(nprocs, double); y1_proc = ARRAY_1D(nprocs, double);
  z0_proc = ARRAY_1D(nprocs, double); z1_proc = ARRAY_1D(nprocs, double);

#ifdef PARALLEL  
  nxp = Gx->np_tot;
  nyp = Gy->np_tot;
  nzp = Gz->np_tot;

  i0 = nghx = Gx->nghost;
  j0 = nghy = Gy->nghost;
  k0 = nghz = Gz->nghost;

  i1 = i0 + Gx->np_int - 1;    
  j1 = j0 + Gy->np_int - 1;
  k1 = k0 + Gz->np_int - 1;

  x0 = Gx->xl[i0]; x1 = Gx->xr[i1];
  y0 = Gy->xl[j0]; y1 = Gy->xr[j1];
  z0 = Gz->xl[k0]; z1 = Gz->xr[k1];

  i0  = Gx->beg; i1 += Gx->beg - nghx;
  j0  = Gy->beg; j1 += Gy->beg - nghy;
  k0  = Gz->beg; k1 += Gz->beg - nghz;

  D_EXPAND(
    MPI_Gather (&x0, 1, MPI_DOUBLE, x0_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&x1, 1, MPI_DOUBLE, x1_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  ,
    
    MPI_Gather (&y0, 1, MPI_DOUBLE, y0_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&y1, 1, MPI_DOUBLE, y1_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  ,
    
    MPI_Gather (&z0, 1, MPI_DOUBLE, z0_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&z1, 1, MPI_DOUBLE, z1_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  )

  D_EXPAND(
    MPI_Gather (&i0, 1, MPI_INT, i0_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather (&i1, 1, MPI_INT, i1_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);  ,
    
    MPI_Gather (&j0, 1, MPI_INT, j0_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather (&j1, 1, MPI_INT, j1_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);  ,
    
    MPI_Gather (&k0, 1, MPI_INT, k0_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather (&k1, 1, MPI_INT, k1_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  )

  print1 ("> Domain Decomposition (%d procs):\n\n", nprocs);

  for (k = 0; k < nprocs; k++){ 
    print1 ("  - Proc # %d, X1: [%f, %f], %d pt\n",k, x0_proc[k], x1_proc[k], Gx->np_int);
    print1 ("              X2: [%f, %f], %d pt\n",    y0_proc[k], y1_proc[k], Gy->np_int);
    #if DIMENSIONS == 3
     print1 ("              X3: [%f, %f], %d pt\n\n",  z0_proc[k], z1_proc[k], Gz->np_int);
    #endif
  }

  MPI_Barrier (MPI_COMM_WORLD);
#endif

/* ---- Free memory ---- */

  FreeArray1D((void *) i0_proc); FreeArray1D((void *) i1_proc);
  FreeArray1D((void *) j0_proc); FreeArray1D((void *) j1_proc);
  FreeArray1D((void *) k0_proc); FreeArray1D((void *) k1_proc);
  FreeArray1D((void *) x0_proc); FreeArray1D((void *) x1_proc);
  FreeArray1D((void *) y0_proc); FreeArray1D((void *) y1_proc);
  FreeArray1D((void *) z0_proc); FreeArray1D((void *) z1_proc);

  print1 ("\n");
}

/* ********************************************************************* */
void CheckConfig()
/*
 *
 *
 * Check if the selected configuration is 
 * allowed.
 *
 *
 *********************************************************************** */
{
  #if DIMENSIONS == 3 

   #if GEOMETRY  == CYLINDRICAL 
    print1 ("\n! Cylindrical coordinates are only 2D.\n");
    print1 ("! Use polar instead.\n");
    QUIT_PLUTO(1);
   #endif

   #if GEOMETRY == SPHERICAL 
    #ifdef SINGLE_STEP
     print1 ("\n ! Spherical 3D only works with RK integrators\n");
     QUIT_PLUTO(1);
    #endif
   #endif

  #endif

  #if DIMENSIONAL_SPLITTING == NO && DIMENSIONS == 1
   #ifndef CH_SPACEDIM
    print1 ("! Cannot integrate a 1-D problem with an unsplit method \n");
    QUIT_PLUTO(1);
   #endif
  #endif

/*
  #if GEOMETRY == SPHERICAL || GEOMETRY == POLAR
   #if TIME_STEPPING != RK2 || TIME_STEPPING != RK3
    print1 (" ! Spherical and Polar geometries work with RK integrators\");
    QUIT_PLUTO(1);
   #endif
  #endif
*/

  #if PARABOLIC_FLUX != 0
   #if ENTROPY_SWITCH == YES
    print1("! Entropy switch not compatible with diffusion operators.\n");
    QUIT_PLUTO(1);
   #endif
  #endif
   
}
