/* ///////////////////////////////////////////////////////////////////// */
/*! 
 \file
 \brief Compute the MHD flux.                                             

  Compute the flux of the conservative MHD equations in the direction 
  given by g_dir.
  This function defines the component of the hyperbolic flux tensor 
  of the standard MHD equations.\n
  In what follows:
  - VXn, MXn, BXn are the velocity, momentum and magnetic field 
    components in the direction given by g_dir (normal, "n")
  - VXt, MXt, BXt and VXb, MXb, BXb are the transverse components
    (tangent "t" and bi-tangent "b").

 \author A. Mignone (mignone@ph.unito.it)
 \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Flux (double **ucons, double **wprim, double *a2, double **bck, 
           double **fx, double *p, int beg, int end)

/*!
 * \param [in]     ucons   1D array of conserved quantities
 * \param [in]     wprim   1D array of primitive quantities
 * \param [in]      a2     1D array of sound speeds
 * \param [in]      bck    1D array of background field values
 * \param [out]     fx     1D array of fluxes (total pressure excluded)
 * \param [out]      p     1D array of pressure values
 * \param [in]      beg    initial index of computation 
 * \param [in]      end    final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int nv, i;
  double vB, ptot;
  double bt1, bt2, bt3;
  double *w, *u;

  for (i = beg; i <= end; i++) {

    w = wprim[i];
    u = ucons[i];
    
    ptot  = 0.5*(EXPAND(w[BX1]*w[BX1] , + w[BX2]*w[BX2], + w[BX3]*w[BX3]));

    #if EOS == IDEAL
     ptot += w[PRS];
    #elif EOS == BAROTROPIC
     ptot += BAROTROPIC_PR(w[RHO]);
    #elif EOS == ISOTHERMAL
     ptot += a2[i]*w[RHO];
    #else
     print ("! FLUX: not defined for this EoS\n");
     QUIT_PLUTO(1);
    #endif

    vB    = EXPAND(w[VX1]*w[BX1] , + w[VX2]*w[BX2], + w[VX3]*w[BX3]);

    #if BACKGROUND_FIELD == YES
     ptot += EXPAND(bck[i][BX1]*w[BX1], + bck[i][BX2]*w[BX2], + bck[i][BX3]*w[BX3]);

     EXPAND(bt1 = w[BXn] + bck[i][BXn];  ,
            bt2 = w[BXt] + bck[i][BXt];  ,
            bt3 = w[BXb] + bck[i][BXb];)

     fx[i][RHO] = u[MXn];
     EXPAND(fx[i][MX1] = w[VXn]*u[MX1] - bt1*w[BX1] - w[BXn]*bck[i][BX1];  ,
            fx[i][MX2] = w[VXn]*u[MX2] - bt1*w[BX2] - w[BXn]*bck[i][BX2];  ,
            fx[i][MX3] = w[VXn]*u[MX3] - bt1*w[BX3] - w[BXn]*bck[i][BX3]; )

     EXPAND(fx[i][BXn] = 0.0;                    ,
            fx[i][BXt] = w[VXn]*bt2 - bt1*w[VXt];  ,
            fx[i][BXb] = w[VXn]*bt3 - bt1*w[VXb]; )
     #if EOS == IDEAL
      fx[i][ENG] = (u[ENG] + ptot)*w[VXn] - bt1*vB;
     #endif
    #else
     fx[i][RHO] = u[MXn];
     EXPAND(fx[i][MX1] = w[VXn]*u[MX1] - w[BXn]*w[BX1];  ,
            fx[i][MX2] = w[VXn]*u[MX2] - w[BXn]*w[BX2];  ,
            fx[i][MX3] = w[VXn]*u[MX3] - w[BXn]*w[BX3]; ) 

     EXPAND(fx[i][BXn] = 0.0;                         ,
            fx[i][BXt] = w[VXn]*w[BXt] - w[BXn]*w[VXt];   ,
            fx[i][BXb] = w[VXn]*w[BXb] - w[BXn]*w[VXb]; )
     #if EOS == IDEAL
      fx[i][ENG] = (u[ENG] + ptot)*w[VXn] - w[BXn]*vB;
     #endif
    #endif

    p[i] = ptot;

    #ifdef GLM_MHD
     fx[i][BXn]      = w[PSI_GLM];
     fx[i][PSI_GLM] = glm_ch*glm_ch*w[BXn];
    #endif
  }
}
