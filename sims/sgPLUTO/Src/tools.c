/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Collection of general-purpose functions.

  This file contains some general Detailed description of the file goes here.

  - LU decomposition functions
  - Debugging / printing / error functions such as Trace(), Show(), 
    Where(), PlutoError
  - function for swapping/detecting endianity
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 20, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

/* ********************************************************************* */
int LUDecomp (double **a, int n, int *indx, double *d)
/*!   
 * Perform LU decomposition. Adapted from Numerical Recipe
 *
 * 
 *********************************************************************** */
#define TINY 1.0e-20;
{
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv;

  vv = ARRAY_1D (n, double);
  *d = 1.0;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if ((temp = fabs (a[i][j])) > big)
        big = temp;
    if (big == 0.0) {
/*      print1 ("! Singular matrix in routine LUDecomp - (i=%d, j=%d)",i,j); */
      return (0);
    }
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = a[i][j];
      for (k = 0; k < i; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++) {
      sum = a[i][j];
      for (k = 0; k < j; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i] * fabs (sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0)
      a[j][j] = TINY;
    if (j != n - 1) {
      dum = 1.0 / (a[j][j]);
      for (i = j + 1; i < n; i++)
        a[i][j] *= dum;
    }
  }
  FreeArray1D(vv);
  return (1); /* -- success -- */
}
#undef TINY

/* ********************************************************************* */
void LUBackSubst (double **a, int n, int *indx, double b[])
/*!
 * Solve a linear system after LUDecomp has been called.
 * Adapted from Numerical Recipe.
 * 
 *********************************************************************** */
{
  int i, ii = 0, ip, j;
  double sum;

/*
  for (i = 0; i < n; i++) {
    if (b[i]!=b[i]) printf(" at input to lubksb: b[%d] = %12.6e\n",i,b[i]);
  }
*/  
  for (i = 0; i < n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii)
      for (j = ii - 1; j <= i - 1; j++)
        sum -= a[i][j] * b[j];
    else if (sum)
      ii = i + 1;
    b[i] = sum;
  }
  for (i = n - 1; i >= 0; i--) {
    sum = b[i];
    for (j = i + 1; j < n; j++) {
      sum -= a[i][j] * b[j];
/* if (sum!=sum) printf(" a[%d][%d] = %12.6e    sum = %12.6e    n = %d\n",i,j,a[i][j],sum,n);*/
    }
    b[i] = sum / a[i][i];
  }
}

/* ********************************************************************* */
void Trace (double xx)
/*!
 * Print a number xx and the number of times it has been called.
 *
 *********************************************************************** */
{
  static int ik;

  printf ("Trace ------> %f ,  %d\n", xx, ++ik);
}
/* ********************************************************************* */
void Show (double **a, int ip)
/*! 
 * Print the component of the array \c a at grid index \c ip  
 *
 *********************************************************************** */
{
  int nv, ix, iy, iz;

  if (g_dir == IDIR) {
    print ("X-sweep");
    ix = ip;
    iy = (g_j == NULL ? -1:*g_j);
    iz = (g_k == NULL ? -1:*g_k);
  } else if (g_dir == JDIR) {
    print ("Y-sweep");
    ix = (g_i == NULL ? -1:*g_i);
    iy = ip;
    iz = (g_k == NULL ? -1:*g_k);
  } else if (g_dir == KDIR) {
    print ("Z-sweep");
    ix = (g_i == NULL ? -1:*g_i);
    iy = (g_j == NULL ? -1:*g_j);
    iz = ip;
  }

  D_SELECT( print (" (%d)> ", ix);     ,
            print (" (%d,%d)> ", ix, iy);  ,
            print (" (%d,%d,%d)> ", ix, iy, iz);  )


  for (nv = 0; nv < NVAR; nv++) {
    print ("%10.4e  ", a[ip][nv]);
  }
  print ("\n");
}

/* ********************************************************************* */
int CheckNaN (double **u, int is, int ie, int id)
/*!
 * Cheeck whether the array \c u contains Not-a-Number
 *  (NaN)
 *
 *********************************************************************** */
{
  int ii, nv, i, j;

  for (ii = is; ii <= ie; ii++) {
  for (nv = 0; nv < NVAR; nv++) {
    if (u[ii][nv] != u[ii][nv]) {
      print (" > NaN found (%d), |", id);
      Show (u, ii);
      QUIT_PLUTO(1);
    }
  }}
  return (0);
}

/* ********************************************************************* */
void Where (int i, Grid *grid)
/*!
 *  Print the location of a particular zone (i,j,k)
 *  in the computational domain.
 *  \note This function must be initialized before using it 
 *        to store grid information. This is done  by calling 
 *        Where(i, grid) the very first time.
 *        Subsequent calls can be then done by simply using 
 *        Where(i,NULL). 
 *
 *********************************************************************** */
{
  int    ii=0, jj=0, kk=0;
  double x1, x2, x3;
  static Grid *grid1, *grid2, *grid3;

/* --------------------------------------------------
    Keep a local copy of grid for subsequent calls
   -------------------------------------------------- */
 
  if (grid != NULL){
    grid1 = grid + IDIR;
    grid2 = grid + JDIR;
    grid3 = grid + KDIR;
    return;
  }

  #ifdef CH_SPACEDIM
   if (g_intStage < 0) return; /* HOT FIX used by CHOMBO
                             (g_intStage = -1) when writing HDF5 file */
  #endif

  if (g_dir == IDIR){
    D_EXPAND(ii = i;, jj = *g_j;, kk = *g_k;)
  }else if (g_dir == JDIR){
    D_EXPAND(ii = *g_i;, jj = i;, kk = *g_k;)
  }else if (g_dir == KDIR){
    D_EXPAND(ii = *g_i;, jj = *g_j;, kk = i;)
  }

  D_EXPAND(
    x1 = grid1->x[ii];  ,
    x2 = grid2->x[jj];  ,
    x3 = grid3->x[kk];
  )

  D_SELECT(
    print ("zone [x1(%d) = %f]",
            ii, grid1->x[ii]);  ,

    print ("zone [x1(%d) = %f, x2(%d) = %f]",
            ii, grid1->x[ii], 
            jj, grid2->x[jj]);  ,

    print ("zone [x1(%d) = %f, x2(%d) = %f, x3(%d) = %f]",
            ii, grid1->x[ii], 
            jj, grid2->x[jj],
            kk, grid3->x[kk]);
  )

  #ifdef CH_SPACEDIM
   print (", Level = %d, %d\n", LEVEL, grid1->level);
   return;
  #endif
  #ifdef PARALLEL
   print (", proc %d\n", prank);
   return;
  #else
   print ("\n");
   return;
  #endif  

}
/* ********************************************************************* */
void PlutoError (int condition, char *str)
/*!
 * If condition is true, issue an error and quit the code.
 *
 *********************************************************************** */
{
  char *str_err="! Error: ";

  if (condition) {
    print (str_err);
    print (str);
    print ("\n");
    QUIT_PLUTO(1);
  }
}

#ifndef CH_SPACEDIM
/* ********************************************************************* */
void print (const char *fmt, ...)
/*!
 * Define print function for the static grid version
 * of PLUTO. The Chombo version is defined in Chombo/amrPLUTO.cpp
 *
 *********************************************************************** */
{
  va_list args;
  va_start(args, fmt);

#if PRINT_TO_FILE == YES
  pluto_log_file = fopen("pluto.log","a"); 
  vfprintf(pluto_log_file, fmt, args);
  fclose(pluto_log_file);
#else
  vprintf(fmt, args);
#endif

  va_end(args);
}
/* ********************************************************************* */
void print1 (const char *fmt, ...)
/*!
 *
 *   Define print1 function
 *
 *********************************************************************** */
{
  va_list args;
  va_start(args, fmt);

  #if PRINT_TO_FILE == YES
   if (prank == 0){
     pluto_log_file = fopen("pluto.log","a"); 
     vfprintf(pluto_log_file,fmt, args);
     fclose(pluto_log_file);
   }
  #else
   if (prank == 0) vprintf(fmt, args);
  #endif
  va_end(args);

}
#endif

/* ********************************************************************* */
void MakeState (State_1D *state)
/*!
 *
 * Allocate memory areas for arrays inside the state
 * structure.
 *
 *
 *********************************************************************** */
{
  state->v       = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->vp      = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->vm      = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->up      = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->um      = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->flux    = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->par_flx = ARRAY_2D(NMAX_POINT, NVAR, double);

  state->src     = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->par_src = ARRAY_2D(NMAX_POINT, NVAR, double);

  state->rhs     = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->press   = ARRAY_1D(NMAX_POINT, double);
  state->bn      = ARRAY_1D(NMAX_POINT, double);
  state->SL      = ARRAY_1D(NMAX_POINT, double);
  state->SR      = ARRAY_1D(NMAX_POINT, double);

/* -- eigenvectors -- */

  state->Lp      = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double);
  state->Rp      = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double);
  state->lambda  = ARRAY_2D(NMAX_POINT, NFLX, double);
  state->lmax    = ARRAY_1D(NVAR, double);

  state->a2   = ARRAY_1D(NMAX_POINT, double);
  state->h    = ARRAY_1D(NMAX_POINT, double);

/*  state->dwlim   = ARRAY_2D(NMAX_POINT, NVAR, double);*/

  state->flag    = ARRAY_1D(NMAX_POINT, unsigned char);

/* --------------------------------------
     define shortcut pointers for
     left and right values with respect
     to the cell center
   -------------------------------------- */
   
  state->vL = state->vp;
  state->vR = state->vm + 1;

  state->uL = state->up;
  state->uR = state->um + 1;

  #ifdef SINGLE_STEP
   state->vh = ARRAY_2D(NMAX_POINT, NVAR, double);
  #else
   state->vh = state->v;
  #endif

/* -- useless ?? -- */

  state->pnt_flx = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->dff_flx = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->vt      = ARRAY_2D(NMAX_POINT, NVAR, double);
}

/* ********************************************************************* */
int IsLittleEndian (void) 
/*!
 * Return 1 if the current architecture has little endian order
 *
 *********************************************************************** */
{
  int TestEndian = 1;
  return *(char*)&TestEndian;
}

/* ********************************************************************* */
void SwapEndian (void *x, const int nbytes) 
/*!
 * Swap the byte order of x.
 * 
 * \param [in] x      pointer to the variable being swapped
 * \param [in] nbytes data type size
 * \return This function has no return value.
 *********************************************************************** */
{
  int k;
  static char Swapped[16];
  char *c;

  c = (char *) x;

  for (k = nbytes; k--; ) Swapped[k] = *(c + nbytes - 1 - k);
  for (k = nbytes; k--; ) c[k] = Swapped[k];

}

/* ********************************************************************* */
void ShowMatrix(double **A, double eps)
/*!
 * Make a nice printing of a 2D matrix \c A
 * Entries with values below eps will display "0.0"
 *
 *
 *********************************************************************** */
{
  int k1,k2;

  print ("----------------------------------------------------------------\n");
  for (k1 = 0; k1 < NFLX; k1++){
    for (k2 = 0; k2 < NFLX; k2++){
      print ("%12.3e   ", fabs(A[k1][k2]) > eps ? A[k1][k2]:0.0);
    }
    printf ("\n");
  }
  print ("----------------------------------------------------------------\n");

}
