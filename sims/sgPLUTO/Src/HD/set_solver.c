#include"pluto.h"

/* ****************************************************************** */
Riemann_Solver *SetSolver (const char *solver)
/*
 *
 * PURPOSE
 *
 *   return a pointer to a riemann solver function
 *
 *
 *
 *
 ********************************************************************* */
{
  
  #ifdef FINITE_DIFFERENCE
   return (&FD_Flux);
  #endif

/* ------------------------------------------------------
       Set Pointers for SOLVERS 
   ------------------------------------------------------ */

  #if EOS == IDEAL
   if ( !strcmp(solver, "two_shock"))  return (&TwoShock_Solver);
   else if (!strcmp(solver, "tvdlf"))  return (&LF_Solver);
   else if (!strcmp(solver, "roe"))    return (&Roe_Solver);
   else if (!strcmp(solver, "ausm+"))  return (&AUSMp_Solver);
   else if (!strcmp(solver, "hlle") || 
            !strcmp(solver, "hll"))    return (&HLL_Solver);
   else if (!strcmp(solver, "hllc"))   return (&HLLC_Solver);
  #elif EOS == ISOTHERMAL
   if (!strcmp(solver, "tvdlf"))       return (&LF_Solver);
   else if (!strcmp(solver, "roe"))    return (&Roe_Solver);
   else if (!strcmp(solver, "hlle") || 
            !strcmp(solver, "hll"))    return (&HLL_Solver);
   else if (!strcmp(solver, "hllc"))   return (&HLLC_Solver);
  #endif

  print1 ("\n! SetSolver: '%s' is not available.\n", solver);
  QUIT_PLUTO(1);
}
