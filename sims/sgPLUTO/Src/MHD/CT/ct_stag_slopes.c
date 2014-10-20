#include "pluto.h"

static double MC_LIM2 (double dp, double dm);

/* ********************************************************************* */
void CT_StoreVelSlopes (EMF *emf, const State_1D *state, int beg, int end)
/*!
 *
 *
 *
 *********************************************************************** */
{
  int i, j, k;

/* ----------------------------------------------------
             Allocate static memory areas 
   ---------------------------------------------------- */

  if (emf->dvx_dx == NULL){

    emf->dvx_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    emf->dvx_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    emf->dvy_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    emf->dvy_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    #if DIMENSIONS == 3
     emf->dvx_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dvy_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

     emf->dvz_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dvz_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dvz_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  }

  if (g_dir == IDIR){

    for (i = beg; i <= end; i++) { 
      emf->dvx_dx[*g_k][*g_j][i] = state->vp[i][VX1] - state->vm[i][VX1];
      emf->dvy_dx[*g_k][*g_j][i] = state->vp[i][VX2] - state->vm[i][VX2];
      #if DIMENSIONS == 3
       emf->dvz_dx[*g_k][*g_j][i] = state->vp[i][VX3] - state->vm[i][VX3];
      #endif
    }

  }else if (g_dir == JDIR){

    for (j = beg; j <= end; j++) {
      emf->dvx_dy[*g_k][j][*g_i] = state->vp[j][VX1] - state->vm[j][VX1];
      emf->dvy_dy[*g_k][j][*g_i] = state->vp[j][VX2] - state->vm[j][VX2];
      #if DIMENSIONS == 3
       emf->dvz_dy[*g_k][j][*g_i] = state->vp[j][VX3] - state->vm[j][VX3];
      #endif
    }

  }else if (g_dir == KDIR){

    for (k = beg; k <= end; k++) {
      emf->dvx_dz[k][*g_j][*g_i] = state->vp[k][VX1] - state->vm[k][VX1];
      emf->dvy_dz[k][*g_j][*g_i] = state->vp[k][VX2] - state->vm[k][VX2];
      #if DIMENSIONS == 3
       emf->dvz_dz[k][*g_j][*g_i] = state->vp[k][VX3] - state->vm[k][VX3];
      #endif
    }
  }
}
/* ********************************************************************* */
void CT_GetStagSlopes (const Data_Arr b, EMF *emf, Grid *grid)
/*!
 * Compute slopes of staggered magnetic fields components
 * Exclude normal derivatives (i.e. dbx_dx).
 *
 *
 *
 *********************************************************************** */
{
  int    i,j,k;
  double ***bx, ***by, ***bz;
  #if BACKGROUND_FIELD == YES
   static double ***B0x, ***B0y, ***B0z;
   double *x,  *y,  *z;
   double *xr, *yr, *zr, Bckgrnd[3];
   double dBp, dBm;

   x  = grid[IDIR].x;  y  = grid[JDIR].x;  z  = grid[KDIR].x;
   xr = grid[IDIR].xr; yr = grid[JDIR].xr; zr = grid[KDIR].xr;
   if (B0x == NULL){
     B0x = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     B0y = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     B0z = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     for (k = KOFFSET; k < NX3_TOT - KOFFSET; k++){
     for (j = JOFFSET; j < NX2_TOT - JOFFSET; j++){
     for (i = IOFFSET; i < NX1_TOT - IOFFSET; i++){
       BackgroundField (xr[i], y[j], z[k], Bckgrnd);
       B0x[k][j][i] = Bckgrnd[0];

       BackgroundField (x[i], yr[j], z[k], Bckgrnd);
       B0y[k][j][i] = Bckgrnd[1];

       BackgroundField (x[i], y[j], zr[k], Bckgrnd);
       B0z[k][j][i] = Bckgrnd[2];
     }}}
   }
  #endif

  D_EXPAND(
    bx = b[BX1s];  , 
    by = b[BX2s];  ,
    bz = b[BX3s];
  )

  if (emf->dbx_dy == NULL){
    emf->dbx_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    emf->dby_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #if DIMENSIONS == 3
     emf->dbx_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dby_dz = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dbz_dx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     emf->dbz_dy = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  }
 
  for (k = KOFFSET; k < NX3_TOT - KOFFSET; k++){
  for (j = JOFFSET; j < NX2_TOT - JOFFSET; j++){
  for (i = IOFFSET; i < NX1_TOT - IOFFSET; i++){

    #if BACKGROUND_FIELD == YES
     dBp =  (B0x[k][j+1][i] + bx[k][j+1][i]) 
          - (B0x[k][j][i]   + bx[k][j][i]);
     dBm =  (B0x[k][j][i]   + bx[k][j][i]) 
          - (B0x[k][j-1][i] + bx[k][j-1][i]);
     emf->dbx_dy[k][j][i] = MC_LIM2(dBp, dBm);

     dBp =  (B0y[k][j][i+1] + by[k][j][i+1]) 
          - (B0y[k][j][i]   + by[k][j][i]);
     dBm =  (B0y[k][j][i]   + by[k][j][i]) 
          - (B0y[k][j][i-1] + by[k][j][i-1]);
     emf->dby_dx[k][j][i] = MC_LIM2(dBp, dBm);
    #else
     emf->dbx_dy[k][j][i] = MC_LIM2(bx[k][j+1][i] - bx[k][j][i], 
                                   bx[k][j][i]   - bx[k][j-1][i]);
     emf->dby_dx[k][j][i] = MC_LIM2(by[k][j][i+1] - by[k][j][i], 
                                   by[k][j][i]   - by[k][j][i-1]);
    #endif
    #if DIMENSIONS == 3
     emf->dbx_dz[k][j][i] = MC_LIM2(bx[k+1][j][i] - bx[k][j][i], 
                                   bx[k][j][i]   - bx[k-1][j][i]);
     emf->dby_dz[k][j][i] = MC_LIM2(by[k+1][j][i] - by[k][j][i], 
                                   by[k][j][i]   - by[k-1][j][i]);
     emf->dbz_dx[k][j][i] = MC_LIM2(bz[k][j][i+1] - bz[k][j][i], 
                                   bz[k][j][i]   - bz[k][j][i-1]);
     emf->dbz_dy[k][j][i] = MC_LIM2(bz[k][j+1][i] - bz[k][j][i], 
                                   bz[k][j][i]   - bz[k][j-1][i]);
    #endif
  }}}

}

/* ********************************************************* */
double MC_LIM2 (double dp, double dm)
/*
 *
 *
 *
 *********************************************************** */
{
  double dc, scrh;

  if (dp*dm < 0.0) return(0.0);

  dc   = 0.5*(dp + dm);
  scrh = 2.0*(fabs(dp) < fabs(dm) ? dp:dm);
  return (fabs(dc) < fabs(scrh) ? dc:scrh);
}


