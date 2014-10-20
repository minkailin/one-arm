/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Initialize PLUTO.

  Initialize() performs a number of initialization tasks before 
  starting the main computation loop.\n
  More precisely, it completes the following sequence of steps:

  - parse command line options 
  - runtime initialization (i.e. pluto.ini)
  - parallel domain decomposition
  - grid generation
  - memory allocation
  - initial conditions
  - set output attributes

  The function GetDecompMode() sets the parallel domain decomposition 
  mode which can be equal to
  
  - AL_AUTO_DECOMP   default;
  - AL_USER_DECOMP   if the -dec n1 [n2] [n3] command line 
                     argument has been given;
                     In this case only, procs[] contains the number
                     of processors in the three directions;
  - AL_MPI_DECOMP   [todo]

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

static int GetDecompMode (Cmd_Line *cmd_line, int procs[]);

/* ********************************************************************* */
void Initialize(int argc, char *argv[], Data *data, 
                Input *input, Grid *grid, Cmd_Line *cmd_line)
/*!
 * Initialize computational grid, domain decomposition and memory
 * allocation. Also, set initial conditions and output attributes.
 *
 * \param [in]    argc  
 * \param [in]    argv
 * \param [in/out] data
 * \param [in]      grid  pointer to an array of Grid structures
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int  nprocs, decomp_mode;
  int  i, j, k, idim, nv;
  int  nx, ny, nz, nghost;
  int  gsize[DIMENSIONS], lsize[DIMENSIONS];
  int  beg[DIMENSIONS], end[DIMENSIONS];
  int  gbeg[DIMENSIONS], gend[DIMENSIONS];
  int  lbeg[DIMENSIONS], lend[DIMENSIONS];
  int  is_gbeg[DIMENSIONS], is_gend[DIMENSIONS];
  int  ghosts[DIMENSIONS];
  int  periods[DIMENSIONS];
  int  pardim[DIMENSIONS], stagdim[DIMENSIONS];
  int  procs[DIMENSIONS];
  char ini_file[128];
  double scrh, dxmin[3], dxming[3];
  Output *output;
  #ifdef PARALLEL
   MPI_Datatype rgb_type;
   MPI_Datatype Float_Vect_type;
  #endif

  sprintf (ini_file,"pluto.ini");

  #ifdef PARALLEL

/* ----------------------------------------
    By default, parallelize all dimensions 
   ---------------------------------------- */

   cmd_line->parallel_dim[IDIR] = YES;
   cmd_line->parallel_dim[JDIR] = YES;
   cmd_line->parallel_dim[KDIR] = YES;

   ParseCmdLineArgs(argc, argv, ini_file, cmd_line);
   ShowConfig();
   if (prank == 0) Setup (input, cmd_line, ini_file);

   MPI_Bcast (input,  sizeof (struct INPUT) , MPI_BYTE, 0, MPI_COMM_WORLD);

   nghost = GetNghost(input);
   MPI_Allreduce (&nghost, &idim, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   nghost = idim;
   
   for (idim = 0; idim < DIMENSIONS; idim++) {
     gsize[idim]   = input->npoint[idim];
     ghosts[idim]  = nghost;
     periods[idim] =    (input->lft_bound_side[idim] == PERIODIC ? 1:0)
                     || (input->lft_bound_side[idim] == SHEARING ? 1:0);
     pardim[idim]  = cmd_line->parallel_dim[idim];
   }

/* -- find parallel decomposition mode -- */

   print1 ("\n> Parallel domain decomposition...\n");
   decomp_mode = GetDecompMode(cmd_line, procs);

/* ---- double distributed array descriptor ---- */

/* SetDistributedArray (SZ, type, gsize, periods, stagdim); 
  return args: beg, end, lsize, lbeg, lend, gbeg, gend, is_gbeg, is_gend
*/

   AL_Sz_init (MPI_COMM_WORLD, &SZ);
   AL_Set_type (MPI_DOUBLE, 1, SZ);
   AL_Set_dimensions (DIMENSIONS, SZ);
   AL_Set_global_dim (gsize, SZ);
   AL_Set_ghosts (ghosts, SZ);
   AL_Set_periodic_dim (periods, SZ);
   AL_Set_parallel_dim (pardim, SZ);

   AL_Decompose (SZ, procs, decomp_mode);
   AL_Get_local_dim (SZ, lsize);
   AL_Get_bounds (SZ, beg, end, ghosts, AL_C_INDEXES);
   AL_Get_lbounds (SZ, lbeg, lend, ghosts, AL_C_INDEXES);
   AL_Get_gbounds (SZ, gbeg, gend, ghosts, AL_C_INDEXES);
   AL_Is_boundary (SZ, is_gbeg, is_gend);

/* ---- float distributed array descriptor ---- */

   AL_Sz_init (MPI_COMM_WORLD, &SZ_float);
   AL_Set_type (MPI_FLOAT, 1, SZ_float);
   AL_Set_dimensions (DIMENSIONS, SZ_float);
   AL_Set_global_dim (gsize, SZ_float);
   AL_Set_ghosts (ghosts, SZ_float);
   AL_Set_periodic_dim (periods, SZ_float);
   AL_Set_parallel_dim (pardim, SZ_float);

   AL_Decompose (SZ_float, procs, decomp_mode);
   AL_Get_local_dim (SZ_float, lsize);
   AL_Get_bounds (SZ_float, beg, end, ghosts, AL_C_INDEXES);
   AL_Get_lbounds (SZ_float, lbeg, lend, ghosts, AL_C_INDEXES);
   AL_Get_gbounds (SZ_float, gbeg, gend, ghosts, AL_C_INDEXES);
   AL_Is_boundary (SZ_float, is_gbeg, is_gend);

/* ---- char distributed array descriptor ---- */

   AL_Sz_init (MPI_COMM_WORLD, &SZ_char);
   AL_Set_type (MPI_CHAR, 1, SZ_char);
   AL_Set_dimensions (DIMENSIONS, SZ_char);
   AL_Set_global_dim (gsize, SZ_char);
   AL_Set_ghosts (ghosts, SZ_char);
   AL_Set_periodic_dim (periods, SZ_char);
   AL_Set_parallel_dim (pardim, SZ_char);

   AL_Decompose (SZ_char, procs, decomp_mode);
   AL_Get_local_dim (SZ_char, lsize);
   AL_Get_bounds (SZ_char, beg, end, ghosts, AL_C_INDEXES);
   AL_Get_lbounds (SZ_char, lbeg, lend, ghosts, AL_C_INDEXES);
   AL_Get_gbounds (SZ_char, gbeg, gend, ghosts, AL_C_INDEXES);
   AL_Is_boundary (SZ_char, is_gbeg, is_gend);

/* ---- Float_Vec distributed array descriptor ---- */

   MPI_Type_contiguous (3, MPI_FLOAT, &Float_Vect_type);
   MPI_Type_commit (&Float_Vect_type);
 
   AL_Sz_init (MPI_COMM_WORLD, &SZ_Float_Vect);
   AL_Set_type (MPI_FLOAT, 3, SZ_Float_Vect);
   AL_Set_dimensions (DIMENSIONS, SZ_Float_Vect);
   AL_Set_global_dim (gsize, SZ_Float_Vect);
   AL_Set_ghosts (ghosts, SZ_Float_Vect);
   AL_Set_periodic_dim (periods, SZ_Float_Vect);
   AL_Set_parallel_dim (pardim, SZ_Float_Vect);

   AL_Decompose (SZ_Float_Vect, procs, decomp_mode);
   AL_Get_local_dim (SZ_Float_Vect, lsize);
   AL_Get_bounds  (SZ_Float_Vect, beg, end, ghosts, AL_C_INDEXES);
   AL_Get_lbounds (SZ_Float_Vect, lbeg, lend, ghosts, AL_C_INDEXES);
   AL_Get_gbounds (SZ_Float_Vect, gbeg, gend, ghosts, AL_C_INDEXES);
   AL_Is_boundary (SZ_Float_Vect, is_gbeg, is_gend);

   for (idim = 0; idim < DIMENSIONS; idim++) {
     grid[idim].nghost      = nghost;
     grid[idim].np_tot      = lsize[idim] + 2*ghosts[idim];
     grid[idim].np_int      = lsize[idim];
     grid[idim].np_tot_glob = input->npoint[idim] + 2*ghosts[idim];
     grid[idim].np_int_glob = input->npoint[idim];
     grid[idim].beg         = beg[idim];
     grid[idim].end         = end[idim];
     grid[idim].gbeg        = gbeg[idim];
     grid[idim].gend        = gend[idim];
     grid[idim].lbeg        = lbeg[idim];
     grid[idim].lend        = lend[idim];
     grid[idim].lbound = input->lft_bound_side[idim]*is_gbeg[idim];
     grid[idim].rbound = input->rgt_bound_side[idim]*is_gend[idim];
     grid[idim].nproc  = procs[idim];
   }

/*  -- Find total number of processors & decomposition mode -- */

   AL_Get_size(SZ, &nprocs);

   #ifdef STAGGERED_MHD

/* ---- x-staggered array descriptor (double) ---- */

    #if DIMENSIONS >= 1
     D_EXPAND(stagdim[IDIR] = AL_TRUE;  ,
              stagdim[JDIR] = AL_FALSE; , 
              stagdim[KDIR] = AL_FALSE;)

     DIM_LOOP(idim) gsize[idim] = input->npoint[idim];
     gsize[IDIR] += 1;

     #ifdef SHEARINGBOX
      periods[IDIR] = 0;
     #endif

     AL_Sz_init (MPI_COMM_WORLD, &SZ_stagx);
     AL_Set_type (MPI_DOUBLE, 1, SZ_stagx);
     AL_Set_dimensions (DIMENSIONS, SZ_stagx);
     AL_Set_global_dim (gsize, SZ_stagx);
     AL_Set_ghosts (ghosts, SZ_stagx);
     AL_Set_staggered_dim(stagdim, SZ_stagx);
     AL_Set_periodic_dim (periods, SZ_stagx);
     AL_Set_parallel_dim (pardim, SZ_stagx);

     AL_Decompose (SZ_stagx, procs, decomp_mode);
     AL_Get_local_dim (SZ_stagx, lsize);
     AL_Get_bounds (SZ_stagx, beg, end, ghosts, AL_C_INDEXES);
     AL_Get_lbounds (SZ_stagx, lbeg, lend, ghosts, AL_C_INDEXES);
     AL_Get_gbounds (SZ_stagx, gbeg, gend, ghosts, AL_C_INDEXES);
     AL_Is_boundary (SZ_stagx, is_gbeg, is_gend);

     #ifdef SHEARINGBOX
      periods[IDIR] = 1;
     #endif
    #endif

    #if DIMENSIONS >= 2
     D_EXPAND(stagdim[IDIR] = AL_FALSE;  ,
              stagdim[JDIR] = AL_TRUE;   , 
              stagdim[KDIR] = AL_FALSE;)

     DIM_LOOP(idim) gsize[idim] = input->npoint[idim];
     gsize[JDIR] += 1;

     AL_Sz_init (MPI_COMM_WORLD, &SZ_stagy);
     AL_Set_type (MPI_DOUBLE, 1, SZ_stagy);
     AL_Set_dimensions (DIMENSIONS, SZ_stagy);
     AL_Set_global_dim (gsize, SZ_stagy);
     AL_Set_ghosts (ghosts, SZ_stagy);
     AL_Set_staggered_dim(stagdim, SZ_stagy);
     AL_Set_periodic_dim (periods, SZ_stagy);
     AL_Set_parallel_dim (pardim, SZ_stagy);

     AL_Decompose (SZ_stagy, procs, decomp_mode);
     AL_Get_local_dim (SZ_stagy, lsize);
     AL_Get_bounds  (SZ_stagy, beg, end, ghosts, AL_C_INDEXES);
     AL_Get_lbounds (SZ_stagy, lbeg, lend, ghosts, AL_C_INDEXES);
     AL_Get_gbounds (SZ_stagy, gbeg, gend, ghosts, AL_C_INDEXES);
     AL_Is_boundary (SZ_stagy, is_gbeg, is_gend);
    #endif

    #if DIMENSIONS == 3
     D_EXPAND(stagdim[IDIR] = AL_FALSE;  ,
              stagdim[JDIR] = AL_FALSE;  , 
              stagdim[KDIR] = AL_TRUE;)

     DIM_LOOP(idim) gsize[idim] = input->npoint[idim];
     gsize[KDIR] += 1;

     AL_Sz_init (MPI_COMM_WORLD, &SZ_stagz);
     AL_Set_type (MPI_DOUBLE, 1, SZ_stagz);
     AL_Set_dimensions (DIMENSIONS, SZ_stagz);
     AL_Set_global_dim (gsize, SZ_stagz);
     AL_Set_ghosts (ghosts, SZ_stagz);
     AL_Set_staggered_dim(stagdim, SZ_stagz);
     AL_Set_periodic_dim (periods, SZ_stagz);
     AL_Set_parallel_dim (pardim, SZ_stagz);

     AL_Decompose (SZ_stagz, procs, decomp_mode);
     AL_Get_local_dim (SZ_stagz, lsize);
     AL_Get_bounds  (SZ_stagz, beg, end, ghosts, AL_C_INDEXES);
     AL_Get_lbounds (SZ_stagz, lbeg, lend, ghosts, AL_C_INDEXES);
     AL_Get_gbounds (SZ_stagz, gbeg, gend, ghosts, AL_C_INDEXES);
     AL_Is_boundary (SZ_stagz, is_gbeg, is_gend);
    #endif

   #endif /* STAGGERED_MHD */

  /* -- find processors coordinates in a Cartesian topology -- */

  {
    int coords[3] = {0,0,0};
    int rank;
    MPI_Comm cartcomm;
    AL_Get_cart_comm(SZ, &cartcomm);
    MPI_Cart_get(cartcomm, DIMENSIONS, procs, periods, coords);
    MPI_Cart_rank(cartcomm, coords, &rank);
    if (rank != prank) {
      print1 ("! Initialize: rank and prank are different\n");
      QUIT_PLUTO(1);
    }
    for (idim = 0; idim < DIMENSIONS; idim++) {
      grid[idim].rank_coord = coords[idim];
    }
  }

  #else  /* if NOT PARALLEL */

/* -----------------------------------------------------
               Serial Initialization
   ----------------------------------------------------- */

   ParseCmdLineArgs (argc, argv, ini_file, cmd_line);
   ShowConfig ();
   Setup (input, cmd_line, ini_file);
   nghost = GetNghost(input);

   for (idim = 0; idim < DIMENSIONS; idim++) {
     grid[idim].nghost  = nghost;
     grid[idim].np_int  = grid[idim].np_int_glob = input->npoint[idim];
     grid[idim].np_tot  = grid[idim].np_tot_glob = input->npoint[idim] + 2*nghost;
     grid[idim].beg     = grid[idim].gbeg = grid[idim].lbeg = nghost;
     grid[idim].end     = grid[idim].gend = grid[idim].lend = (grid[idim].lbeg - 1) + grid[idim].np_int;
     grid[idim].lbound  = input->lft_bound_side[idim];
     grid[idim].rbound  = input->rgt_bound_side[idim];
     grid[idim].nproc   = 1;
   }
   
   nprocs = 1;

  #endif

/* ---------------------------------------------------
                Grid Generation
   --------------------------------------------------- */

  SetGrid (input, grid);
  Where (-1, grid);     /* -- store grid inside the "Where" 
                              function for subsequent calls -- */
  if (cmd_line->makegrid == YES) {
    print1 ("\n> Done < \n");
    QUIT_PLUTO(0);
  }

/* ------------------------------------------------
         initialize global variables
   ------------------------------------------------ */

  g_dt   = input->first_dt;
  g_time = g_maxMach = 0.0;
  g_maxRiemannIter  = 0;
  g_usedMem = 0;

  IBEG = grid[IDIR].lbeg; IEND = grid[IDIR].lend;
  JBEG = grid[JDIR].lbeg; JEND = grid[JDIR].lend;
  KBEG = grid[KDIR].lbeg; KEND = grid[KDIR].lend;

  NX1 = grid[IDIR].np_int;
  NX2 = grid[JDIR].np_int;
  NX3 = grid[KDIR].np_int;

  NX1_TOT = grid[IDIR].np_tot; 
  NX2_TOT = grid[JDIR].np_tot;
  NX3_TOT = grid[KDIR].np_tot;

/* ---------------------------------------
     get the maximum number of points 
     among all directions
   --------------------------------------- */

  NMAX_POINT = MAX(NX1_TOT, NX2_TOT);
  NMAX_POINT = MAX(NMAX_POINT, NX3_TOT);

/* --------------------------------------------------------------------
        FIND THE MINUM PHYSICAL CELL LENGTH IN EACH DIMENSIONS
   -------------------------------------------------------------------- */

  for (idim = 0; idim < DIMENSIONS; idim++)  dxmin[idim] = 1.e30;

  for (i = IBEG; i <= IEND; i++) {
  for (j = JBEG; j <= JEND; j++) {
  for (k = KBEG; k <= KEND; k++) {

    scrh = Length_1(i, j, k, grid);
    dxmin[IDIR] = MIN (dxmin[IDIR], scrh);

    scrh = Length_2(i, j, k, grid);
    dxmin[JDIR] = MIN (dxmin[JDIR], scrh);
 
    scrh = Length_3(i, j, k, grid); 
    dxmin[KDIR] = MIN (dxmin[KDIR], scrh);
     
  }}}

  #ifdef PARALLEL
   MPI_Allreduce (dxmin, dxming, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   dxmin[IDIR] = dxming[IDIR];
   dxmin[JDIR] = dxming[JDIR];
   dxmin[KDIR] = dxming[KDIR];
  #endif

  grid[IDIR].dl_min = dxmin[IDIR];
  grid[JDIR].dl_min = dxmin[JDIR];
  grid[KDIR].dl_min = dxmin[KDIR];

/* ------------------------------------------------------
    Copy user defined parameters into global array g_inputParam 
   ------------------------------------------------------ */

  for (nv = 0; nv < USER_DEF_PARAMETERS; nv++) g_inputParam[nv] = input->aux[nv];

/* ------------------------------------------------------------
          Allocate memory for 3D data arrays
   ------------------------------------------------------------ */

  print1 ("\n> memory allocation\n");
  data->Vc = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
/*  data->Uc = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double); */

  #ifdef STAGGERED_MHD
   data->Vs = ARRAY_1D(DIMENSIONS, double ***);
   D_EXPAND(
     data->Vs[BX1s] = ArrayBox( 0, NX3_TOT-1, 0, NX2_TOT-1,-1, NX1_TOT-1); ,
     data->Vs[BX2s] = ArrayBox( 0, NX3_TOT-1,-1, NX2_TOT-1, 0, NX1_TOT-1); ,
     data->Vs[BX3s] = ArrayBox(-1, NX3_TOT-1, 0, NX2_TOT-1, 0, NX1_TOT-1);)
  #endif  

  #if UPDATE_VECTOR_POTENTIAL == YES 
   D_EXPAND(                                                  ,
     data->Ax3 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  , 
     data->Ax1 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  
     data->Ax2 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  
   )
  #endif

  data->flag = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, unsigned char);

/* ------------------------------------------------------------
              Assign initial conditions
   ------------------------------------------------------------ */

   StartupSG(data, grid);
   Startup (data, grid);
   if(g_inputParam[sg_on]>0.0) SGDensity(data, grid);
   TotalMass(data, grid);
/*   TestPoisson(data, grid); 
   exit(1);
*/

/* ------------------------------------------------------------ 
    Set output attributes (names, images, number of outputs...)
   ------------------------------------------------------------ */

  SetOutput (data, input);

/* -----------------------------------
      print normalization units
   ----------------------------------- */

  ShowUnits();  

  print1 ("> Number of processors: %d\n",nprocs);
  D_EXPAND(print1 ("> Typical proc size:    %d",grid[IDIR].np_int);  ,
           print1 (" X %d", grid[JDIR].np_int);                      ,
           print1 (" X %d", grid[KDIR].np_int);)
  print1 ("\n");
  #ifdef PARALLEL
   print1 ("> Parallel Directions: ");
   D_EXPAND(if (pardim[IDIR]) print1 (" X1");  ,
            if (pardim[JDIR]) print1 ("/X2");  ,
            if (pardim[KDIR]) print1 ("/X3");)
   print1 ("\n");
  #endif   
  if (cmd_line->show_dec) ShowDomainDecomposition (nprocs, grid);
}

#ifdef PARALLEL
/* ********************************************************************* */
int GetDecompMode (Cmd_Line *cmd_line, int procs[])
/*!
 * Returns the parallel domain decomposition mode.
 *
 * \param [in]  cmd_line  pointer to the Cmd_Line structure
 * \param [out] procs     an array of integers giving the number
 *                        of processors in each direction only if
 *                        the -dec command line option has been given
 *
 *  \return  The decomposition mode:
 *  
 *   - AL_AUTO_DECOMP   defaults
 *   - AL_USER_DECOMP   if the -dec n1 [n2] [n3] command line 
 *                      argument has been given;
 *                      In this case only, procs[] contains the number
 *                      of processors in the three directions;
 *   -  AL_MPI_DECOMP   if [todo]
 * \todo AL_MPI_DECOMP mode
 *********************************************************************** */ 
{
  long int npx = 1, npy = 1, npz = 1;
  long int npp;

  D_EXPAND(npx = cmd_line->nproc[IDIR];  ,
           npy = cmd_line->nproc[JDIR];  ,
           npz = cmd_line->nproc[KDIR];)

/* ------------------------------------------------
    if -dec has not been given
    decomposition will be set to be AL_AUTO_DECOMP
   ------------------------------------------------ */

  if (npx == -1 || npy == -1 || npz == -1){
    return AL_AUTO_DECOMP;
  }

/* -----------------------------------------------
    enter in user decomp mode if the number of
    processors along all directions has been given
   ------------------------------------------------ */

  if (npx > 0 && npy > 0 && npz > 0){
    procs[IDIR] = npx;
    procs[JDIR] = npy;
    procs[KDIR] = npz;
    return AL_USER_DECOMP;
  }

  print1 ("! GetDecompMode: invalid decomposition mode");
  QUIT_PLUTO(1);
}
#endif

