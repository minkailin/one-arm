/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Shearing-Box module header file

  The Shearing-Box module header file contains basic macro definitions, 
  function prototypes and declaration of global variables used by the 
  sheraring-box module.
  The variable ::sb_q and ::sb_Omega are the most important ones and 
  \b must be defined and initialized in your init.c in order to 
  configure your shearing-box problem.

  Optionally, the order of interpolation (default is 2) at physical 
  boundaries may be changed using the ::SB_ORDER macro.

  The additional macros ::SB_SYMMETRIZE_HYDRO, ::SB_SYMMETRIZE_EY and
  ::SB_SYMMETRIZE_EZ may be set to YES/NO to enable/disable enforcement
  of conservation at the radial (x) boundaries.

  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)
  \date   Aug 16, 2012
  \todo Check if sb_vy and sb_Ly are really needed as global variables.
*/
/* ///////////////////////////////////////////////////////////////////// */

/*! Sets the order of interpolation at physical boundaries (1, 2 or 3).*/
#if INTERPOLATION == LINEAR || INTERPOLATION == WENO3 || INTERPOLATION == LimO3
 #define SB_ORDER              2    
#else
 #define SB_ORDER              3    
#endif                                 

#ifndef FARGO
 #define SB_SYMMETRIZE_HYDRO  YES  /**< Symmetrize the hydrodynamical fluxes 
         at the left and right x-boundaries in order to enforce conservation of
         hydrodynamic variables like density, momentum and energy
        (no magnetic field). Default is YES.*/

 #define SB_SYMMETRIZE_EY    (YES  && (DIMENSIONS == 3)) /**< Symmetrize the
         y-component of the electric field at the left and right x-boundaries
         to enforce conservation of magnetic field (only in 3D). */

 #define SB_SYMMETRIZE_EZ     YES  /**< Symmetrize the z-component of electric
         field at the left and right x-boundaries to enforce conservation of
         magnetic field. */

 #define SB_FORCE_EMF_PERIODS NO  /**< Force periodicity at y- and z-
                                       boundaries. */
#else
 #define SB_SYMMETRIZE_HYDRO  NO  
 #define SB_SYMMETRIZE_EY    (NO  && (DIMENSIONS == 3)) 
 #define SB_SYMMETRIZE_EZ     NO 
 #define SB_FORCE_EMF_PERIODS NO 
#endif

/* ----------------------------
    Global variables
   ---------------------------- */

extern double sb_q;  /**< The shear parameter, \f$\DS q = -\HALF\frac{d\log
              \Omega^2}{d\log R} \f$. The explicit numerical value and the
              variable definition should be set inside your Init() function. */

extern double sb_Omega; /**< Disk local orbital frequency \f$ \Omega_0 = 
              \Omega(R_0)\f$. The explicit numerical value and the variable
              definition should be set inside your Init() function.  */

extern double sb_vy;  /**< velocity offset (>0), in USERDEF_BOUNDARY */ 

#define sb_A (-0.5*sb_Omega*sb_q)  /**< Short-hand definition for the Oort
                                        constant \f$ A = -q\Omega_0/2 \f$. */

void SB_Boundary (const Data *, int, Grid *) ;
void SB_ShearingInterp (double *, double *, double, int, Grid *);
void SB_CorrectFluxes  (Data_Arr, double, double, Grid *);
#ifdef STAGGERED_MHD
void SB_CorrectEMF (EMF *, Data_Arr, double, Grid *);
#endif
int  SB_JSHIFT (int);
void SB_SaveFluxes (State_1D *, Grid *);
void SB_SetBoundaryVar(double ***, RBox *, int, double, Grid *);
void SB_ShiftBoundaryVar(double ***, RBox *, int, double, Grid *);
void SB_FillBoundaryGhost(double ***, RBox *, int, int, Grid *);

#ifdef PARALLEL
 void ExchangeX (real *, real *, int, Grid *);
#endif

