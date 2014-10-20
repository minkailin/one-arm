#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{

  int i, j, k;  
  double ***pot; 

  pot = GetUserVar("pot");

  DOM_LOOP(k,j,i){
   pot[k][j][i] = d->Phi[k][j][i];
  }
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}
