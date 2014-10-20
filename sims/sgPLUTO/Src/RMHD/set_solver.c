#include"pluto.h"

/* ****************************************************************** */
Riemann_Solver *SetSolver (const char *solver)
/*
 *
 *
 *
 *
 *
 *
 ********************************************************************* */
{
  
/* ------------------------------------------------------
       Set Pointers for SOLVERS 
   ------------------------------------------------------ */

  if (!strcmp(solver, "tvdlf")) {

    return (&LF_Solver);

  }else if (!strcmp(solver, "hlle") || 
            !strcmp(solver, "hll")) {

    return (&HLL_Solver);

  }else if (!strcmp(solver, "hllc")) {

    return (&HLLC_Solver);

  }else if (!strcmp(solver, "hlld")) {

    return (&HLLD_Solver);

  }

  print1 ("\n ! SetSolver: '%s' is not available.\n", solver);
  QUIT_PLUTO(1);

}
