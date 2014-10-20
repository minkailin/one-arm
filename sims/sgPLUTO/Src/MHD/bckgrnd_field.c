#include "pluto.h"

#if BACKGROUND_FIELD == YES
/* ******************************************************************** */
double **GetBackgroundField (int beg, int end, int where, Grid *grid)
/*
 *
 *
 *
 *
 *
 ********************************************************************* */
{
  int i;
  double *x1, *x2, *x3;
  static real **bck_fld;

/* ----------------------------------------------------
         Check for incompatibilities 
   ---------------------------------------------------- */
  
  #if (TIME_STEPPING != RK2) && (TIME_STEPPING != RK3)
   print1 ("! Background field splitting works with RK integrators ONLY \n");
   QUIT_PLUTO(1);
  #endif
  #if MHD_FORMULATION != CONSTRAINED_TRANSPORT
   print1 ("! Background field splitting works with CT ONLY \n");
   QUIT_PLUTO(1);
  #else
   #if (CT_EMF_AVERAGE != ARITHMETIC) && (CT_EMF_AVERAGE != UCT_HLL)
    print1 ("! Background field splitting works with ARITHMETIC or UCT_HLL averages only\n");
    QUIT_PLUTO(1);
   #endif
  #endif
  #if RESISTIVE_MHD != NO 
    print1 ("! Background field splitting works with Ideal MHD only\n");
    QUIT_PLUTO(1);
  #endif  

  if (bck_fld == NULL) {
    bck_fld = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

/* ----------------------------------
     get pointers to coordinates 
   ---------------------------------- */
   
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  if (g_dir == IDIR){
    if (where == FACE_CENTER) x1 = grid[IDIR].xr;
    for (i = beg; i <= end; i++){
      BackgroundField (x1[i],x2[*g_j],x3[*g_k], bck_fld[i] + BX);
    }
  }else if (g_dir == JDIR){
    if (where == FACE_CENTER) x2 = grid[JDIR].xr;
    for (i = beg; i <= end; i++){
      BackgroundField (x1[*g_i],x2[i],x3[*g_k], bck_fld[i] + BX);
    }
  }else{
    if (where == FACE_CENTER) x3 = grid[KDIR].xr;
    for (i = beg; i <= end; i++){
      BackgroundField (x1[*g_i],x2[*g_j],x3[i], bck_fld[i] + BX);
    }
  }    

  return(bck_fld);
}
#endif
