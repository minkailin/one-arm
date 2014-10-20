/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file
  \brief Collect various thermodynamic functions for the relativistic
         MHD module.
                    
  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SoundSpeed2 (double **v, double *cs2, double *h, int beg, int end,
                  int pos, Grid *grid)
/*!
 * Define the square of the sound speed.
 * 
 * \param [in]    u   1D array of primitive quantities
 * \param [out] cs2   1D array containing the square of the sound speed
 * \param [in]    h   1D array of enthalpy values
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 * \param [in]  pos   an integer specifying the spatial position 
 *                    inside the cell (only for spatially-dependent EOS)
 * \param [in]  grid  pointer to Grid structure
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int   i;
  double  theta;

  Enthalpy(v, h, beg, end);

  for (i = beg; i <= end; i++) {
    theta = v[i][PRS]/v[i][RHO];
    #if EOS == IDEAL
     cs2[i] = g_gamma*theta/h[i];
    #elif EOS == TAUB
     cs2[i] = theta/(3.0*h[i])*(5.0*h[i] - 8.0*theta)/(h[i] - theta);
    #endif
  }
}

/* ********************************************************************* */
void Enthalpy (double **v, double *h, int beg, int end)
/*!
 * Compute the enthalpy.
 *
 * \param [in]    v   1D array of primitive quantities
 * \param [in]    h   1D array of enthalpy values
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int   i;
  double  theta, gmmr;

  gmmr = g_gamma/(g_gamma - 1.0);
  for (i = beg; i <= end; i++) {
    theta = v[i][PRS]/v[i][RHO];
    #if EOS == IDEAL
     h[i] = 1.0 + gmmr*theta;
    #elif EOS == TAUB
     h[i] = 2.5*theta + sqrt(2.25*theta*theta + 1.0);
    #endif
  }
}

/* ********************************************************************* */
void Entropy (double **v, double *s, int beg, int end)
/*!
 * Compute the entropy.
 * 
 * \param [in]    v   1D array of primitive quantities
 * \param [in]    s   1D array of entropy values
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int   i;
  double rho, th;

  for (i = beg; i <= end; i++) {
    rho = v[i][RHO];
    #if EOS == IDEAL
     s[i] = v[i][PRS]/pow(rho,g_gamma);
    #elif EOS == TAUB
     th   = v[i][PRS]/rho;
     s[i] = v[i][PRS]/pow(rho,5./3.)*(1.5*th + sqrt(2.25*th*th + 1.0));
    #endif
  }
}
