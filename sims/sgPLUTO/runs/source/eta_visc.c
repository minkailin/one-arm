#include "pluto.h"
#include <stdio.h>

/*---------------------------------------------------------------------------*/
/*---- Specification of explicit first and second viscosity coefficients ----*/
/*---------------------------------------------------------------------------*/

void eta_visc_func(real *v, real x1, real x2, real x3, 
                   real *eta1_visc, real *eta2_visc )
{
  *eta1_visc = v[RHO]*kinematic_visc(x1, x2);
  *eta2_visc = 0.0;
}
