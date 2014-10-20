/* ############################################################
      
     FILE:     cooling.h

     PURPOSE:  contains common definitions for the 
               whole CODE

   ############################################################ */

#define NIONS  1
#define DNEUT   NFLX    /* -- this is the 'conservative' index -- */
#define FNEUT   DNEUT   /* -- this is the 'primitive' index -- */

real GetMaxRate (real *, real *, real);
real MeanMolecularWeight  (real *);
double H_MassFrac (void);
real COMP_EUIL (real, real, real *);
void Radiat (real *, real *);










