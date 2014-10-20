/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Write tabulated 1D or 2D ascii data files.

  WriteTabArray() provides ascii output for 2-D and 1-D data 
  in tabulated multi-column format:
  \verbatim
      .     .         .         .             .
      .     .         .         .             .
      .     .         .         .             .
    x1(i) x2(j)   var_1(i,j)  var_2(i,j)   var_3(i,j) ...
      .     .         .         .             .
      .     .         .         .             .
      .     .         .         .             .
  \endverbatim
 
  Blocks with different x2 coordinates are separated by a blank row.
  One file (with all variables) per output is written to disk. 
  This format does not work in parallel mode and can be used for simple
  serial runs.
 
  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 21, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void WriteTabArray (Output *output, char *filename, Grid *grid)
/*!
 * Write tabulated array.
 *
 * \param [in] output   the output structure corresponding to a given 
 *                      format
 * \param [in] filename the output file name
 * \param [in] grid     pointer to an array of Grid structures
 *********************************************************************** */
{
  int  nv, i, j, k;
  FILE *fout;

  #ifdef PARALLEL
   print1 ("! Tabulated output disabled in parallel mode\n");
   return;
  #endif

  #if DIMENSIONS == 3
   print1 ("! Tab output not supported in 3D\n");
   return;
  #endif

/* ------------------------------------------
      dump arrays in ascii format to disk 
   ------------------------------------------ */

  fout = fopen (filename, "w");
  k = 0;
  IDOM_LOOP (i){
    JDOM_LOOP(j){
      fprintf (fout, "%f %f ", grid[IDIR].x[i], grid[JDIR].x[j]);
      for (nv = 0; nv < output->nvar; nv++) {
        if (output->dump_var[nv]) 
          fprintf (fout, "%12.6e ", output->V[nv][k][j][i]);
      }
      fprintf (fout, "\n");
    }
    fprintf (fout, "\n");
  }

  fclose (fout);
}

