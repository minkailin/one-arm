/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Enforce conservation at the X1 boundaries in the shearing-box
         module.

  This file provides a set of functions for storing and then modifying 
  the upwind fluxes computed during the Riemann solver at the leftmost 
  and righmost physical boundaries.
  These tasks are performed only when either SB_SYMMETRIZE_HYDRO, 
  SB_SYMMETRIZE_EY, SB_SYMMETRIZE_EZ flags are set to YES in 
  Src/MHD/shearingbox.h and are useful to avoid loss of conservation 
  in the hydrodynamical variables (density, momentum, energy) 
  and/or magnetic fields. 
  
  This is first achieved by calling the SB_SAVE_FLUXES() function during
  the time stepping scheme, with the purpose of storing the leftmost and 
  rightmost conservative fluxes in the x-direction into the static
  arrays FluxL[] and FluxR[]. 

  These fluxes are then subsequently used by SB_CORRECT_FLUXES() which 
  interpolates the fluxes and properly correct leftmost and rightmost 
  cell-centered flow quantities to ensure conservation.
  
  The treatment of staggered magnetic field is done similarly by 
  SB_CORRECT_EMF().
  
  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.musicanisi@cineca.it)\n
           G. Bodo (bodo@oato.inaf.it)
           
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static double ****FluxL; /**< Array of fluxes at the left x-boundary. */ 
static double ****FluxR; /**< Array of fluxes at the right x-boundary. */ 

#if MHD_FORMULATION == DIV_CLEANING
/* #define NVLAST  (BX-1) */
 #define NVLAST  (NVAR-1)
#define SKIPNV()  if (nv == PSI_GLM) continue
#else  
 #define NVLAST  (NVAR-1) 
 #define SKIPNV()  
#endif

/* ********************************************************************* */
void SB_SaveFluxes (State_1D *state, Grid *grid)
/*!
 * Store leftmost and rightmost conservative fluxes in the x direction 
 * into FluxL and FluxR. 
 *
 * \param [in] state pointer to a State_1D structure
 * \param [in] grid pointer to an array of Grid structures
 *
 * \return This function has no return value.
 *********************************************************************** */
{
  int i, j, k, nv;

  #if SB_SYMMETRIZE_HYDRO == NO
   return;
  #endif

  if (FluxL == NULL){
    FluxL = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, 1, double);
    FluxR = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, 1, double);
  }

/* ---------------------------------------------------- 
    Save Fluxes on the LEFT physical boundary
   ---------------------------------------------------- */
  
  if (grid[IDIR].lbound != 0){
    if (g_dir == IDIR){
      state->flux[IBEG - 1][MX1] += state->press[IBEG - 1];
      for (nv = 0; nv <= NVLAST; nv++){
        FluxL[nv][*g_k][*g_j][0] = state->flux[IBEG - 1][nv];
SKIPNV();
        state->flux[IBEG - 1][nv] = 0.0;
      }
      state->press[IBEG - 1] = 0.0;  
    } 
  }

/* ---------------------------------------------------- 
    Save Fluxes on the RIGHT pysical boundary
   ---------------------------------------------------- */
  
  if (grid[IDIR].rbound != 0){
    if (g_dir == IDIR){
      state->flux[IEND][MX1] += state->press[IEND];    
      for (nv = 0; nv <= NVLAST; nv++){
        FluxR[nv][*g_k][*g_j][0] = state->flux[IEND][nv];
SKIPNV();
        state->flux[IEND][nv] = 0.0;
      }
      state->press[IEND] = 0.0;  
    }
  }
}

/* ********************************************************************* */
void SB_CorrectFluxes (Data_Arr U, double t, double dt, Grid *grid)
/*!
 * Interpolate x-fluxes FLuxL and FluxR and properly correct leftmost
 * and rightmost cells to ensure conservation.
 *
 * \param [in,out] U data array containing cell-centered quantities
 * \param [in] t  the time step at which fluxes have  to be computed
 * \param [in] dt the time step being used in the integrator stage 
 * \param [in] grid pointer to array of Grid structures
 *********************************************************************** */
{
  int    i, j, k, nv;
  double dtdx;
  static double **fL, **fR;
  static State_1D state;

  #if SB_SYMMETRIZE_HYDRO == NO
   return;
  #endif

  if (fL == NULL){
    fL = ARRAY_2D(NVAR, NMAX_POINT, double);
    fR = ARRAY_2D(NVAR, NMAX_POINT, double);
  }

  dtdx = dt/grid[IDIR].dx[IBEG];

  #ifdef PARALLEL
   ExchangeX (FluxL[0][0][0], FluxR[0][0][0], NVAR*NX2_TOT*NX3_TOT, grid);
  #endif

/* --------------------------------------------------------------------- */
/*! \note 
    In parallel and if there's more than one processor in the
    x direction, we modify the values of FluxL and FluxR
    independently since the two boundary sides are handled by
    different processors.
    Conversely, when a processor owns both sides, we need to
    copy at least one of the two arrays before doing interpolation
    since their original values are lost when setting boundary
    conditions.                                                          */
/* --------------------------------------------------------------------- */

  if (grid[IDIR].lbound != 0){  /* ---- Left x-boundary ---- */

    RBox box;
    static double ****Flux1;

  /* -- set grid ranges of FluxR and exchange b.c. -- */

    box.ib = 0; box.ie = 0; 
    box.jb = 0; box.je = NX2_TOT-1; 
    box.kb = 0; box.ke = NX3_TOT-1;

  /* -- make a copy of FluxR to avoid overwriting one nproc_x = 1 -- */

    if (grid[IDIR].nproc == 1){
      if (Flux1 == NULL) Flux1 = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, 1, double);
      for (nv = 0; nv <= NVLAST; nv++)
        BOX_LOOP((&box), k, j, i) Flux1[nv][k][j][i] = FluxR[nv][k][j][i];
    }else{
      Flux1 = FluxR;
    }
    
  /* -- set boundary conditions on Flux1 -- */

    for (nv = 0; nv <= NVLAST; nv++)
      SB_SetBoundaryVar(Flux1[nv], &box, X1_BEG, t, grid);

    for (k = KBEG; k <= KEND; k++){
      for (nv = 0; nv <= NVLAST; nv++){
        for (j = JBEG; j <= JEND; j++)
          fL[nv][j] = 0.5*(Flux1[nv][k][j][0] + FluxL[nv][k][j][0]);
      }
      for (j = JBEG; j <= JEND; j++){
        #if EOS == IDEAL
         fL[ENG][j] += 0.5*sb_vy*(fL[MX2][j] + 0.5*sb_vy*fL[RHO][j]);
        #endif
        #ifdef GLM_MHD
         fL[BX2][j] -= 0.5*sb_vy*fL[PSI_GLM][j]/glm_ch/glm_ch;
/*         fL[BX2][j] -= 0.5*sb_vy*FluxL[PSI_GLM][k][j]/glm_ch/glm_ch; */
        #endif
        fL[MX2][j] += 0.5*sb_vy*fL[RHO][j];
        for (nv = 0; nv <= NVLAST; nv++){
SKIPNV();
          U[k][j][IBEG][nv] += dtdx*fL[nv][j];
        }
      }
    }
  }   

  if (grid[IDIR].rbound != 0){  /* ---- Right x-boundary ---- */

  /* -- set grid ranges of FluxL and exchange b.c. -- */

    RBox box;
    box.ib = 0; box.ie = 0;
    box.jb = 0; box.je = NX2_TOT-1;
    box.kb = 0; box.ke = NX3_TOT-1;

    for (nv = 0; nv <= NVLAST; nv++)
      SB_SetBoundaryVar(FluxL[nv], &box, X1_END, t, grid);
    
    for (k = KBEG; k <= KEND; k++){
      for (nv = 0; nv <= NVLAST; nv++){
        for (j = JBEG; j <= JEND; j++)
          fR[nv][j] = 0.5*(FluxL[nv][k][j][0] + FluxR[nv][k][j][0]);
      }
      for (j = JBEG; j <= JEND; j++){
        #if EOS == IDEAL
         fR[ENG][j] += 0.5*sb_vy*(-fR[MX2][j] + 0.5*sb_vy*fR[RHO][j]);
        #endif
        #ifdef GLM_MHD
         fR[BX2][j] += 0.5*sb_vy*fR[PSI_GLM][j]/glm_ch/glm_ch;
/*         fR[BX2][j] += 0.5*sb_vy*FluxR[PSI_GLM][k][j]/glm_ch/glm_ch; */
        #endif
        fR[MX2][j] -= 0.5*sb_vy*fR[RHO][j];
        for (nv = 0; nv <= NVLAST; nv++){
SKIPNV();
          U[k][j][IEND][nv] -= dtdx*fR[nv][j];
        }
      }
    }  
  }
}

#ifdef STAGGERED_MHD
/* ********************************************************************* */
void SB_CorrectEMF (EMF *emf, Data_Arr Vs, double dt, Grid *grid)
/*! 
 *
 * \param [in,out] emf pointer to EMF structure
 * \param [in]     Vs  4D array of staggered fields
 * \param [in]     dt  time step
 * \param [in]   grid  pointer to array of Grid structures
 *********************************************************************** */
{
  int    i, j, k, nghost;
  int    dimx[3] = {1, 0, 0};
  int    dimy[3] = {0, 1, 0};
  int    dimz[3] = {0, 0, 1};
  double   bxS, bxN, esym;
  double   t, ddt;
  static double **bxL, ***dbxL_lim, *dbxL, ***eyL, ***ezL;
  static double **bxR, ***dbxR_lim, *dbxR, ***eyR, ***ezR;

  #if    (SB_SYMMETRIZE_EY == NO) && (SB_SYMMETRIZE_EZ == NO) \
      && (SB_FORCE_EMF_PERIODS == NO)
   return;
  #endif

  nghost = grid[IDIR].nghost;
  t      = g_time + g_dt;
  #ifdef SINGLE_STEP
   ddt = 1.0 - sb_vy*dt/grid[JDIR].dx[JBEG];
  #else
   ddt = 1.0;
  #endif

  if (ezL == NULL){
    eyL    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    eyR    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);

    ezL    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    ezR    = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);

    bxL    = ARRAY_2D(NX3_TOT, NX2_TOT, double);
    bxR    = ARRAY_2D(NX3_TOT, NX2_TOT, double);

    dbxL     = ARRAY_1D(NMAX_POINT, double);
    dbxR     = ARRAY_1D(NMAX_POINT, double);
    dbxL_lim = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
    dbxR_lim = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
  }

/* -- compute Bx slopes on left side -- */

  if (grid[IDIR].lbound != 0){

  /* -- Store Ey, Ez, Bx on the left and compute slopes -- */

    KTOT_LOOP(k) JTOT_LOOP(j){
      D_EXPAND(                                         ,
               ezL[k][j][0] = emf->ez[k][j][IBEG - 1];
               bxL[k][j]    = Vs[BX1s][k][j][IBEG - 1];  ,
               eyL[k][j][0] = emf->ey[k][j][IBEG - 1];)
    }

    for (k = KBEG; k <= KEND; k++){
      for (j = 1; j <= JEND + nghost; j++) dbxL[j] = bxL[k][j] - bxL[k][j-1]; 
      for (j = JBEG-1; j <= JEND+1; j++){
        dbxL_lim[k][j][0] = VAN_LEER(dbxL[j+1], dbxL[j]);
      }
    } 
  }

/* -- compute Bx slopes on right side -- */

  if (grid[IDIR].rbound != 0){

  /* -- Store Ey, Ez, Bx on the right and compute slopes -- */

    KTOT_LOOP(k) JTOT_LOOP(j){
      D_EXPAND(                                     ,
               ezR[k][j][0] = emf->ez[k][j][IEND];
               bxR[k][j]    = Vs[BX1s][k][j][IEND];  ,
               eyR[k][j][0] = emf->ey[k][j][IEND]; )
    }

    for (k = KBEG; k <= KEND; k++){
      for (j = 1; j <= JEND + nghost; j++) dbxR[j] = bxR[k][j] - bxR[k][j-1];
      for (j = JBEG-1; j <= JEND+1; j++){
        dbxR_lim[k][j][0] = VAN_LEER(dbxR[j+1], dbxR[j]);
      }
    }
  }

/* ---------------------------------------------------------------------
            exchange data between processors
   --------------------------------------------------------------------- */

  #ifdef PARALLEL
   D_EXPAND(                                                           ,
            ExchangeX (ezL[0][0], ezR[0][0], NX3_TOT*NX2_TOT, grid); 
            ExchangeX (dbxL_lim[0][0], dbxR_lim[0][0], NX3_TOT*NX2_TOT, grid);   ,
            ExchangeX (eyL[0][0], eyR[0][0], NX3_TOT*NX2_TOT, grid); )
  #endif

/* ----------------------------------------------------------------
    Symmetrize Ey and Ez on the left
    Similarly to the hydro fluxes, if a processor owns both sides,
    we copy one of the two arrays before doing interpolation since
    its original value is lost after b.c.
   ---------------------------------------------------------------- */

  if (grid[IDIR].lbound != 0){
    RBox box;
    static double ***e1, ***dbx1_lim;
   
    box.ib = 0; box.ie = 0;
    box.jb = 0; box.je = NX2_TOT-1;
    box.kb = 0; box.ke = NX3_TOT-1;
    #if SB_SYMMETRIZE_EY == YES

   /* -- make a copy of FluxR to avoid overwriting one nproc_x = 1 -- */

     if (grid[IDIR].nproc == 1){
       if (e1 == NULL) e1 = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
       BOX_LOOP((&box), k, j, i) e1[k][j][i] = eyR[k][j][i];
     }else{
       e1 = eyR;
     }

     SB_SetBoundaryVar(e1, &box, X1_BEG, t, grid);
     for (k = KBEG - 1; k <= KEND; k++){
       for (j = JBEG; j <= JEND; j++) 
         emf->ey[k][j][IBEG - 1] = 0.5*(eyL[k][j][0] + e1[k][j][0]);
     }
    #endif

    #if SB_SYMMETRIZE_EZ == YES
     if (grid[IDIR].nproc == 1){
       if (e1       == NULL) e1       = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
       if (dbx1_lim == NULL) dbx1_lim = ARRAY_3D(NX3_TOT, NX2_TOT, 1, double);
       BOX_LOOP((&box), k, j, i) {
         dbx1_lim[k][j][i] = dbxR_lim[k][j][i];
               e1[k][j][i] =      ezR[k][j][i];
       }
     }else{
       dbx1_lim = dbxR_lim;
             e1 = ezR;
     }

     SB_SetBoundaryVar(dbx1_lim, &box, X1_BEG, g_time, grid);
     SB_SetBoundaryVar(      e1, &box, X1_BEG,         t, grid);
     for (k = KBEG; k <= KEND; k++){ 
       for (j = JBEG - 1; j <= JEND + 1; j++)
         dbxL[j] = 0.5*(dbx1_lim[k][j][0] + dbxL_lim[k][j][0]);

       for (j = JBEG - 1; j <= JEND; j++){
         bxS = bxL[k][j] + 0.5*dbxL[j]*ddt;
         emf->ez[k][j][IBEG-1] = 0.5*(ezL[k][j][0] + e1[k][j][0] + sb_vy*bxS);
       }
     }
    #endif
  }

/* ------------------------------------------------------
         Symmetrize Ey, and Ez on the right
   ------------------------------------------------------ */

  if (grid[IDIR].rbound != 0){
    RBox box;
    box.ib = 0; box.ie = 0;
    box.jb = 0; box.je = NX2_TOT-1;
    box.kb = 0; box.ke = NX3_TOT-1;

    #if SB_SYMMETRIZE_EY == YES
     SB_SetBoundaryVar(eyL, &box, X1_END, t, grid);
     for (k = KBEG - 1; k <= KEND; k++){
       for (j = JBEG; j <= JEND; j++) 
         emf->ey[k][j][IEND] = 0.5*(eyR[k][j][0] + eyL[k][j][0]);
     } 
    #endif

    #if SB_SYMMETRIZE_EZ == YES
     SB_SetBoundaryVar(dbxL_lim, &box, X1_END, g_time, grid);
     SB_SetBoundaryVar(     ezL, &box, X1_END,         t, grid);

     for (k = KBEG; k <= KEND; k++){ 
       for (j = JBEG - 1; j <= JEND + 1; j++){
         dbxR[j] = 0.5*(dbxL_lim[k][j][0] + dbxR_lim[k][j][0]);
       }

       for (j = JBEG - 1; j <= JEND; j++){
         bxN = bxR[k][j + 1] - 0.5*dbxR[j + 1]*ddt;
         emf->ez[k][j][IEND] = 0.5*(ezR[k][j][0] + ezL[k][j][0] - sb_vy*bxN);
       }
     }
    #endif
  }

/* --------------------------------------------------
                Force Periodicity 
   -------------------------------------------------- */

  #if DIMENSIONS == 3 && SB_FORCE_EMF_PERIODS == YES
               
  /* -- Ex at Z faces: force periodicty  -- */
                                                                     
   #ifdef PARALLEL
    MPI_Barrier (MPI_COMM_WORLD);
    AL_Exchange_dim (emf->ex[0][0], dimz, SZ);
    MPI_Barrier (MPI_COMM_WORLD);
   #else
    for (j = JBEG - 1; j <= JEND; j++){
    for (i = IBEG    ; i <= IEND; i++){
      esym = 0.5*(emf->ex[KBEG - 1][j][i] + emf->ex[KEND][j][i]);
      emf->ex[KBEG - 1][j][i] = esym;
      emf->ex[KEND    ][j][i] = esym;
    }}
   #endif
                                                                                                                                                      
   /*  Ex at Y faces: force periodicity  */
                                                                                                                                                      
   for (k = KBEG - 1; k <= KEND; k++){
   for (i = IBEG    ; i <= IEND; i++){
     esym = 0.5*(emf->ex[k][JBEG - 1][i] + emf->ex[k][JEND][i]);
     emf->ex[k][JBEG - 1][i] = esym;
     emf->ex[k][JEND    ][i] = esym;
   }}
     
   /*  Ey at Z faces: force periodicity   */
                                                                                                                                                      
   #ifdef PARALLEL
    MPI_Barrier (MPI_COMM_WORLD);
    AL_Exchange_dim (emf->ey[0][0], dimz, SZ);
    MPI_Barrier (MPI_COMM_WORLD);
   #else
    for (j = JBEG    ; j <= JEND; j++){
    for (i = IBEG - 1; i <= IEND; i++){
      esym = 0.5*(emf->ey[KBEG - 1][j][i] + emf->ey[KEND][j][i]);
      emf->ey[KBEG - 1][j][i] = esym;
      emf->ey[KEND    ][j][i] = esym;
    }}
   #endif

   /*  Ez at Y faces: force periodicity   */
 
   for (k = KBEG    ; k <= KEND; k++){
   for (i = IBEG - 1; i <= IEND; i++){
     esym = 0.5*(ez[k][JBEG - 1][i] + ez[k][JEND][i]);
     ez[k][JBEG - 1][i] = esym;
     ez[k][JEND    ][i] = esym;
   }}
 
  #endif 
}
#endif

#ifdef PARALLEL
/* ********************************************************************* */
void ExchangeX (double *bufL, double *bufR, int nel, Grid *grid)
/*!
 * Send bufL owned by the processor at X1_BEG
 * to the processor at X1_END;
 * Send bufR owned by the processor at X1_END
 * to the processor at X1_BEG;
 *
 *********************************************************************** */
{
  static int dest = -1;
  int stag = 1, rtag = 1;
  MPI_Comm cartcomm;
  MPI_Status istat;
  MPI_Request req;
  static  int nprocs[3], periods[3], coords[3];

  AL_Get_cart_comm(SZ, &cartcomm);
  MPI_Barrier (MPI_COMM_WORLD);

/* --------------------------------------
     get rank of the processor lying 
     on the opposite side of x-domain
   -------------------------------------- */

  if (dest == -1){
    if (grid[IDIR].lbound != 0){
      MPI_Cart_get(cartcomm, 3, nprocs, periods, coords);
      coords[0] += nprocs[0] - 1;
      MPI_Cart_rank (cartcomm, coords, &dest);
    }

    if (grid[IDIR].rbound != 0){
      MPI_Cart_get(cartcomm, 3, nprocs, periods, coords);
      coords[0] += - nprocs[0] + 1;
      MPI_Cart_rank (cartcomm, coords, &dest); 
    }
  }

  if (grid[IDIR].lbound != 0){
    if (prank != dest){
      MPI_Sendrecv (bufL, nel, MPI_DOUBLE, dest, stag,
                    bufR, nel, MPI_DOUBLE, dest, rtag,
                    MPI_COMM_WORLD, &istat);
    }
  }

  if (grid[IDIR].rbound != 0){
    if (prank != dest){
      MPI_Sendrecv (bufR, nel, MPI_DOUBLE, dest, stag,
                    bufL, nel, MPI_DOUBLE, dest, rtag,
                    MPI_COMM_WORLD, &istat);
    }
  }

  MPI_Barrier (MPI_COMM_WORLD);
}

#endif
