/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Resistive MHD module header file.

  Contains prototypes for the resistive MHD module.

  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date   Sep 13, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

void ResistiveFlux (Data_Arr, double **, double **, int, int, Grid *);

void GetCurrent (Data_Arr, real **, Grid *);

void ETA_Func (real *, real, real, real, real *);

void ADD_OHM_HEAT (const Data *d, double, Grid *);
