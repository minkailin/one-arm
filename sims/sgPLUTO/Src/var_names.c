#include "pluto.h"

/* ****************************************************************** */
void SetDefaultVarNames(Output *output)
/*
 *
 *  PURPOSE
 *
 *    Set file names for I/O
 *
 *
 ******************************************************************** */
{
  int nv;

/* ----------------------------------------------
    Physics module file names; 
    these pertain to the physics module ONLY
   ---------------------------------------------- */

  output->var_name[RHO] = "rho";
  EXPAND(output->var_name[VX1] = "vx1";  ,
         output->var_name[VX2] = "vx2";  ,
         output->var_name[VX3] = "vx3";)
  #if EOS != ISOTHERMAL && EOS != BAROTROPIC
   output->var_name[PRS] = "prs";
  #endif

  #if PHYSICS == MHD || PHYSICS == RMHD
   EXPAND(output->var_name[BX1] = "bx1";  ,
          output->var_name[BX2] = "bx2";  ,
          output->var_name[BX3] = "bx3";)
  #endif
  
  /* (staggered field names are set in SetOutput) */

  #ifdef GLM_MHD
   output->var_name[PSI_GLM] = "psi_glm";
  #endif
  
/* ------------------------------------------------
                   Tracers 
   ------------------------------------------------ */

  for (nv = TRC; nv < TRC + NTRACER; nv++){
    sprintf (output->var_name[nv],"tr%d",nv - TRC + 1);
  } 

  #if ENTROPY_SWITCH == YES
   sprintf (output->var_name[ENTR],"entropy");
  #endif

/* ------------------------------------------------
               Cooling vars
   ------------------------------------------------ */

  #if COOLING == MINEq
  {

   static char *ion_name[] = {"fHI", "fHeI", "fHeII" 
                       C_EXPAND("fCI","fCII", "fCIII", "fCIV", "fCV")
                       N_EXPAND("fNI","fNII", "fNIII", "fNIV", "fNV")
                       O_EXPAND("fOI","fOII", "fOIII", "fOIV", "fOV")
                      Ne_EXPAND("fNeI","fNeII", "fNeIII", "fNeIV", "fNeV")
                       S_EXPAND("fSI","fSII", "fSIII", "fSIV", "fSV")
                      Fe_EXPAND("fFeI", "fFeII", "fFeIII")};

   for (nv = 0; nv < NIONS; nv++) output->var_name[NFLX+nv] = ion_name[nv];  
/*
   output->var_name[HI]    = "fHI";
   output->var_name[HeI]   = "fHeI";
   output->var_name[HeII]  = "fHeII";

   output->var_name[CI]    = "fCI";  
   output->var_name[CII]   = "fCII"; 
   output->var_name[CIII]  = "fCIII";
   output->var_name[CIV]   = "fCIV";
   output->var_name[CV]    = "fCV";
   output->var_name[NI]    = "fNI";
   output->var_name[NII]   = "fNII";
   output->var_name[NIII]  = "fNIII";
   output->var_name[NIV]   = "fNIV";
   output->var_name[NV]    = "fNV";
   output->var_name[OI]    = "fOI";
   output->var_name[OII]   = "fOII";
   output->var_name[OIII]  = "fOIII";
   output->var_name[OIV]   = "fOIV";
   output->var_name[OV]    = "fOV";
   output->var_name[NeI]   = "fNeI";
   output->var_name[NeII]  = "fNeII";
   output->var_name[NeIII] = "fNeIII";
   output->var_name[NeIV]  = "fNeIV";
   output->var_name[NeV]   = "fNeV";
   output->var_name[SI]    = "fSI";
   output->var_name[SII]   = "fSII";
   output->var_name[SIII]  = "fSIII";
   output->var_name[SIV]   = "fSIV";
   output->var_name[SV]    = "fSV";
*/
/*
   #if INCLUDE_Fe == YES
    output->var_name[FeI]   = "fFeI";
    output->var_name[FeII]  = "fFeII";
    output->var_name[FeIII] = "fFeIII";
   #endif 
*/

  }
  #elif COOLING == SNEq

   output->var_name[FNEUT] = "fneut";

  #endif

}
