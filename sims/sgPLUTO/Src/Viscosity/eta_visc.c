/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <stdio.h>
/* ************************************************************************** */
void eta_visc_func(real *v, real x1, real x2, real x3, 
                   real *eta1_visc, real *eta2_visc )
/*! 
 * \brief Calculate first and second viscosity coefficients as functions of data and coordinates      
 *
 *    \param [in]      v  pointer to data array containing cell-centered quantities
 *    \param [in]      x1 real, coordinate value 
 *    \param [in]      x2 real, coordinate value 
 *    \param [in]      x3 real, coordinate value 
 *    \param [in, out] eta1_visc pointer to first viscous coefficient
 *    \param [in, out] eta2_visc pointer to second viscous coefficient
 *    \return This function has no return value.
 * ************************************************************************** */

{
 *eta1_visc = 0.0;
 *eta2_visc = 0.0;
}
