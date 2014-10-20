/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Main output driver.

  WriteData() is the main driver for writing data arrays in any of
  the available formats (binary, VTK, HDF5, etc...).  

  - For .dbl, .flt or .vtk file formats, access to binary files is 
    provided by the functions in bin_io.c.
  - HDF5 files are handled by hdf5_io.c.
  - image files are handled by write_img.c
  - tabulated ascii files are handled by write_tab.c

  This function also updates the corresponding .out file associated 
  with the output data format.

  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)

  \date   Aug 24, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void WriteData (const Data *d, Output *output, Grid *grid)
/*!
 * Write data to disk using any of the available formats.
 *
 * \param [in] d      pointer to PLUTO Data structre 
 * \param [in] output the output structure corresponding to a given
 *                    format
 * \param [in] grid   pointer to an array of Grid structures
 *********************************************************************** */
{
  int    i, j, k, nv;
  int    single_file;
  size_t dsize;
  char   filename[128], sline[512];
  static int last_computed_var = -1;
  void *Vpt;
  FILE *fout, *fbin;
  time_t tbeg, tend;
  long long offset;

/* -----------------------------------------------
    Convert from conservative to primitive
   ----------------------------------------------- */
/*
  ConsToPrim3D(d);
*/
/* -----------------------------------------------------------
                Increment the file number 
   ----------------------------------------------------------- */

  output->nfile++;

  print1 ("> Writing file #%d (%s) to disk...", output->nfile, output->ext);

  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   if (prank == 0)time(&tbeg);
  #endif

/* --------------------------------------------------------
            Get user var if necessary 
   -------------------------------------------------------- */

  if (last_computed_var != g_stepNumber && d->Vuser != NULL) {
    ComputeUserVar (d, grid);
    last_computed_var = g_stepNumber;
  }

/* --------------------------------------------------------
            Select the output type 
   -------------------------------------------------------- */

  if (output->type == DBL_OUTPUT) {

  /* ------------------------------------------------------------------- */
  /*! - \b DBL output:
        Double-precision data files can be written using single or
        multiple file mode. 
        - for single file, serial: we open the file just once before
          the main variable loop, dump variables and then close.
        - for single file, parallel the distributed array descriptor sz is
          different for cell-centered or staggered data type and we
          thus have to open and close the file after each variable
          has been dumped.
        - when writing multiple files we open, write to and close the
          file one each loop cycle.
        \note In all cases, the pointer to the data array that has to be 
              written must be cast into (void *) and the starting index of 
              the array must be zero.
  */
  /* ------------------------------------------------------------------- */

    int sz;
    single_file = strcmp(output->mode,"single_file") == 0;
    dsize = sizeof(double);

    if (single_file){  /* -- single output file -- */

      sprintf (filename, "data.%04d.%s", output->nfile, output->ext);
      offset = 0;
      #ifndef PARALLEL
       fbin = OpenBinaryFile (filename, 0, "w");
      #endif
      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;

        if      (output->stag_var[nv] == -1) {  /* -- cell-centered data -- */
          sz = SZ;
          Vpt = (void *)output->V[nv][0][0];
        } else if (output->stag_var[nv] == 0) { /* -- x-staggered data -- */
          sz  = SZ_stagx;
          Vpt = (void *)(output->V[nv][0][0]-1);
        } else if (output->stag_var[nv] == 1) { /* -- y-staggered data -- */
          sz = SZ_stagy;
          Vpt = (void *)output->V[nv][0][-1];
        } else if (output->stag_var[nv] == 2) { /* -- z-staggered data -- */
           sz = SZ_stagz;
           Vpt = (void *)output->V[nv][-1][0];
        }
        #ifdef PARALLEL
         fbin = OpenBinaryFile (filename, sz, "w");
         AL_Set_offset(sz, offset);
        #endif
        WriteBinaryArray (Vpt, dsize, sz, fbin, output->stag_var[nv]);
        #ifdef PARALLEL
         offset = AL_Get_offset(sz);
         CloseBinaryFile(fbin, sz);
        #endif
      }
      #ifndef PARALLEL
       CloseBinaryFile(fbin, sz);
      #endif

    }else{              /* -- multiple files -- */

      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;
        sprintf (filename, "%s.%04d.%s", output->var_name[nv], 
                                         output->nfile, output->ext);

        if      (output->stag_var[nv] == -1) {  /* -- cell-centered data -- */
          sz = SZ;
          Vpt = (void *)output->V[nv][0][0];
        } else if (output->stag_var[nv] == 0) { /* -- x-staggered data -- */
          sz  = SZ_stagx;
          Vpt = (void *)(output->V[nv][0][0]-1);
        } else if (output->stag_var[nv] == 1) { /* -- y-staggered data -- */
          sz = SZ_stagy;
          Vpt = (void *)output->V[nv][0][-1];
        } else if (output->stag_var[nv] == 2) { /* -- z-staggered data -- */
           sz = SZ_stagz;
           Vpt = (void *)output->V[nv][-1][0];
        }
        fbin = OpenBinaryFile (filename, sz, "w");
        WriteBinaryArray (Vpt, dsize, sz, fbin, output->stag_var[nv]);
        CloseBinaryFile (fbin, sz);
      }
    }

  } else if (output->type == FLT_OUTPUT) {

  /* ----------------------------------------------------------
                 FLT output for cell-centered data
     ---------------------------------------------------------- */

    single_file = strcmp(output->mode,"single_file") == 0;

    if (single_file){  /* -- single output file -- */

      sprintf (filename, "data.%04d.%s", output->nfile, output->ext);
      fbin = OpenBinaryFile (filename, SZ_float, "w");
      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;
        Vpt = (void *)(Convert_dbl2flt(output->V[nv],0))[0][0];
        WriteBinaryArray (Vpt, sizeof(float), SZ_float, fbin, 
                          output->stag_var[nv]);
      }
      CloseBinaryFile(fbin, SZ_float);

    }else{              /* -- multiple files -- */

      for (nv = 0; nv < output->nvar; nv++) {
        if (!output->dump_var[nv]) continue;
        sprintf (filename, "%s.%04d.%s", output->var_name[nv], 
                                         output->nfile, output->ext);

        fbin = OpenBinaryFile (filename, SZ_float, "w");
        Vpt = (void *)(Convert_dbl2flt(output->V[nv],0))[0][0];
        WriteBinaryArray (Vpt, sizeof(float), SZ_float, fbin, 
                          output->stag_var[nv]);
        CloseBinaryFile (fbin, SZ_float);
      }
    }

  }else if (output->type == DBL_H5_OUTPUT || output->type == FLT_H5_OUTPUT){

  /* ------------------------------------------------------
             HDF5 output (single/double precision)
     ------------------------------------------------------ */

    #ifdef USE_HDF5 
     single_file = YES;
     WriteHDF5 (output, grid);
    #else
     print1 ("! WriteData: HDF5 library not available\n");
     return;
    #endif

  }else if (output->type == VTK_OUTPUT) { 

  /* ------------------------------------------------------------------- */
  /*! - \b VTK output:  
      in order to enable parallel writing, files must be closed and
      opened again for scalars, since the distributed array descriptors 
      used by ArrayLib (Float_Vect) and (float) are different. This is
      done using the AL_Get_offset() and AL_Set_offset() functions.      */
  /* ------------------------------------------------------------------- */
    
    single_file = strcmp(output->mode,"single_file") == 0;
    sprintf (filename, "data.%04d.%s", output->nfile, output->ext);
    if (single_file){  /* -- single output file -- */

      fbin  = OpenBinaryFile(filename, SZ_Float_Vect, "w");
      WriteVTK_Header(fbin, grid);
      for (nv = 0; nv < output->nvar; nv++) {  /* -- write vectors -- */
        if (output->dump_var[nv] != VTK_VECTOR) continue;
        WriteVTK_Vector (fbin, output->V + nv, output->var_name[nv], grid);
      }

      #ifdef PARALLEL
       offset = AL_Get_offset(SZ_Float_Vect);
       CloseBinaryFile(fbin, SZ_Float_Vect);
       fbin  = OpenBinaryFile(filename, SZ_float, "w");
       AL_Set_offset(SZ_float, offset);
      #endif
      
      for (nv = 0; nv < output->nvar; nv++) { /* -- write scalars -- */
        if (output->dump_var[nv] != YES) continue;
        WriteVTK_Scalar (fbin, output->V[nv], output->var_name[nv], grid);
      }
      CloseBinaryFile(fbin, SZ_float);

    }else{          /* -- multiple output files -- */

      for (nv = 0; nv < output->nvar; nv++) { /* -- write vectors -- */
        if (output->dump_var[nv] != VTK_VECTOR) continue;
        if (strcmp(output->var_name[nv],"vx1") == 0) {
          sprintf (filename, "vfield.%04d.%s", output->nfile, output->ext);
        }else if (strcmp(output->var_name[nv],"bx1") == 0) {
          sprintf (filename, "bfield.%04d.%s", output->nfile, output->ext);
        }else{
          print1 ("! WriteData: unknown vector type in VTK output\n"); 
          QUIT_PLUTO(1);
        }

        fbin = OpenBinaryFile(filename, SZ_Float_Vect, "w");
        WriteVTK_Header(fbin, grid);
        WriteVTK_Vector(fbin, output->V + nv, output->var_name[nv], grid);
        CloseBinaryFile(fbin, SZ_Float_Vect);
      }

      for (nv = 0; nv < output->nvar; nv++) {  /* -- write scalars -- */
        if (output->dump_var[nv] != YES) continue;
        sprintf (filename, "%s.%04d.%s", output->var_name[nv], 
                                         output->nfile,
                                         output->ext);
        fbin = OpenBinaryFile(filename, SZ_Float_Vect, "w");
        WriteVTK_Header(fbin, grid);
        #ifdef PARALLEL
         offset = AL_Get_offset(SZ_Float_Vect);
         CloseBinaryFile(fbin, SZ_Float_Vect);
         fbin  = OpenBinaryFile(filename, SZ_float, "w");
         AL_Set_offset(SZ_float, offset);
        #endif
        WriteVTK_Scalar(fbin, output->V[nv], output->var_name[nv], grid);
        CloseBinaryFile (fbin, SZ_float);
      }
    }

  }else if (output->type == TAB_OUTPUT) { 

  /* ------------------------------------------------------
                        TAB OUTPUT 
     ------------------------------------------------------ */

    single_file = YES;
    sprintf (filename,"data.%04d.%s", output->nfile,output->ext);
    WriteTabArray (output, filename, grid);

  }else if (output->type == PPM_OUTPUT) { 

  /* ------------------------------------------------------
                       PPM OUTPUT
     ------------------------------------------------------ */

    single_file = NO;
    for (nv = 0; nv < output->nvar; nv++) {
      if (!output->dump_var[nv]) continue;
      sprintf (filename, "%s.%04d.%s", output->var_name[nv], 
                                       output->nfile,
                                       output->ext);
      WritePPM (output->V[nv], output->var_name[nv], filename, grid);
    }

  }else if (output->type == PNG_OUTPUT) { 

  /* ------------------------------------------------------
                        PNG OUTPUT 
     ------------------------------------------------------ */

    #ifdef USE_PNG
     single_file = NO;
     for (nv = 0; nv < output->nvar; nv++) {
       if (!output->dump_var[nv]) continue;
       sprintf (filename, "%s.%04d.%s", output->var_name[nv], 
                                        output->nfile,
                                        output->ext);
       WritePNG (output->V[nv], output->var_name[nv], filename, grid);
     }
    #else
     print1 ("! PNG library not available\n");
     return;
    #endif

  }

/* -------------------------------------------------------------
           Update corresponding ".out" file
   ------------------------------------------------------------- */

  sprintf (filename,"%s.out",output->ext);

  if (prank == 0) {
    if (output->nfile == 0) {
      fout = fopen (filename, "w");
    }else {
      fout = fopen (filename, "r+");
      for (nv = 0; nv < output->nfile; nv++) fgets (sline, 512, fout);
      fseek (fout, ftell(fout), SEEK_SET);
    }

  /* -- write a multi-column file -- */

    fprintf (fout, "%d %12.6e %12.6e %ld ",
             output->nfile, g_time, g_dt, g_stepNumber);

    if (single_file) fprintf (fout,"single_file ");
    else             fprintf (fout,"multiple_files ");

    if (IsLittleEndian()) fprintf (fout, "little ");
    else                 fprintf (fout, "big ");

    for (nv = 0; nv < output->nvar; nv++) { 
      if (output->dump_var[nv]) fprintf (fout, "%s ", output->var_name[nv]);
    }

    fprintf (fout,"\n");
    fclose (fout);
  }

  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   if (prank == 0){
     time(&tend);
     print1 (" [%5.2f sec]",difftime(tend,tbeg));
   }
  #endif
  print1 ("\n");
}

#ifdef USE_ASYNC_IO
static float ****Vflt;
static int perf_output[16] = {0};
/* ****************************************************************** */
void Async_BegWriteData (const Data *d, Output *output, Grid *grid)
/*!
 *
 * PURPOSE:
 *
 *  Write data to disk using binary format and asyncronous MPI functions: 
 *   dbl, flt.
 *
 *  \author CINECA (g.muscianisi@cineca.it), A. Mignone (mignone@ph.unito.it)
 *
 ******************************************************************** */
{
  int  i, j, k, nv;
  size_t dsize;
  char   filename[128];
  static int last_computed_var = -1;

/* -----------------------------------------------------------
                Increment the file number 
   ----------------------------------------------------------- */

  output->nfile++;
  print1 ("> Writing file #%d (%s) to disk [async: beg]...\n",
             output->nfile, output->ext);

/* --------------------------------------------------------
            Get user var if necessary 
   -------------------------------------------------------- */

  if (last_computed_var != g_stepNumber && d->Vuser != NULL) {
    ComputeUserVar (d, grid);
    last_computed_var = g_stepNumber;
  }

/* ------------------------------------------------------
                  DBL/FLT OUTPUTS 
   ------------------------------------------------------ */

  if (output->type == DBL_OUTPUT) {
    dsize = sizeof(double);
    perf_output[DBL_OUTPUT] = 1;
  } else{
    dsize = sizeof(float);
    perf_output[FLT_OUTPUT] = 1;
  }    
  
  sprintf (filename, "data.%04d.%s", output->nfile, output->ext);

  if (dsize == sizeof(double)) AL_File_open(filename, SZ);
  if (dsize == sizeof(float))  AL_File_open(filename, SZ_float);

  if (dsize == sizeof(double)){
    AL_Write_array_begin ((void *)output->V[0][0][0], SZ, output->stag_var,
                                  output->dump_var, output->nvar);
  }
  if (dsize == sizeof(float)){
    if (Vflt == NULL){
      Vflt = ARRAY_4D(output->nvar, NX3_TOT, NX2_TOT, NX1_TOT, float);
    }
  
    /* similar to CONVERT_TO_FLOAT, with swap_endian disabled */
    
    for (nv = 0; nv < output->nvar; nv++){
      DOM_LOOP(k,j,i) Vflt[nv][k][j][i] = (float)output->V[nv][k][j][i]; 
    }
    AL_Write_array_begin ((void *)Vflt[0][0][0], SZ_float, 
                          output->stag_var, output->dump_var, output->nvar);
  }
}

/* ****************************************************************** */
void Async_EndWriteData (Input *ini)
/*!
 *
 * PURPOSE:
 *
 *  Writing data completition using binary format and asyncronous MPI functions: 
 *   dbl, flt.
 *
 * \author CINECA (g.muscianisi@cineca.it), A. Mignone (mignone@ph.unito.it)
 *
 ******************************************************************** */
{
  char filename[128], sline[512];
  FILE *fout;
  int nv;
  Output *output;

  if (perf_output[DBL_OUTPUT]){   /* asynchronous dbl output */
    output = ini->output + 0;
    print1 ("> Writing file #%d (%s) to disk [async: end]...\n",
             output->nfile, output->ext);
    AL_Write_array_end((void *)output->V[0][0][0], SZ);
    AL_File_close(SZ);

    sprintf (filename,"%s.out",output->ext);   
    if (prank == 0) {
      if (output->nfile == 0) {
        fout = fopen (filename, "w");
      }else{
        fout = fopen (filename, "r+");
        for (nv = 0; nv < output->nfile; nv++) fgets (sline, 512, fout);
           fseek (fout, ftell(fout), SEEK_SET);
      }

      /* -- write a multi-column file -- */
      fprintf (fout, "%d %8.3e %8.3e %ld ", output->nfile, g_time, g_dt, g_stepNumber);
      fprintf (fout,"single_file ");

      if (IsLittleEndian()) fprintf (fout, "little ");
      else                 fprintf (fout, "big ");

      for (nv = 0; nv < output->nvar; nv++) {
         if (output->dump_var[nv]) fprintf (fout, "%s ", output->var_name[nv]);
      }
      fprintf (fout,"\n");
      fclose (fout);
    }
    perf_output[DBL_OUTPUT] = 0;
  }
  if (perf_output[FLT_OUTPUT]){  /* asynchronous flt output */
    output =ini->output + 1;
    print1 ("> Writing file #%d (%s) to disk [async: end]...\n",
             output->nfile, output->ext);

    AL_Write_array_end((void *)Vflt[0][0][0], SZ_float);
    AL_File_close(SZ_float);

    sprintf (filename,"%s.out",output->ext);  
    if (prank == 0) {
      if (output->nfile == 0) {
        fout = fopen (filename, "w");
      }else {
        fout = fopen (filename, "r+");
        for (nv = 0; nv < output->nfile; nv++) fgets (sline, 512, fout);
           fseek (fout, ftell(fout), SEEK_SET);
      }

      /* -- write a multi-column file -- */
      fprintf (fout, "%d %8.3e %8.3e %ld ", 
                output->nfile, g_time, g_dt, g_stepNumber);
      fprintf (fout,"single_file ");

      if (IsLittleEndian()) fprintf (fout, "little ");
      else                 fprintf (fout, "big ");

      for (nv = 0; nv < output->nvar; nv++) {
         if (output->dump_var[nv]) fprintf (fout, "%s ", output->var_name[nv]);
      }
      fprintf (fout,"\n");
      fclose (fout);
    }
    perf_output[FLT_OUTPUT] = 0;
  }

}
#endif /* USE_ASYNC_IO */

/* ********************************************************************* */
void ConsToPrim3D(const Data *d, int where)
/*!
 *  Convert the data array d->Uc into d->Vc inside 
 *  the computational domain.
 *********************************************************************** */
{
  int i,j,k,nv;
  static double **v;
  static unsigned char *flag;
  
  if (v == NULL) {
    v    = ARRAY_2D(NMAX_POINT,NVAR,double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);  
  }
  
  if (where == 0){
    KDOM_LOOP(k){
    JDOM_LOOP(j){
      ConsToPrim(d->Uc[k][j], v, IBEG, IEND, flag);
      IDOM_LOOP(i){
        for (nv = NVAR; nv--;  ) d->Vc[nv][k][j][i] = v[i][nv];
      }
    }} 
  }else if (where == 1){
    KTOT_LOOP(k){
    JTOT_LOOP(j){
      ConsToPrim(d->Uc[k][j], v, 0, NX1_TOT-1, flag);
      ITOT_LOOP(i){
        for (nv = NVAR; nv--;  ) d->Vc[nv][k][j][i] = v[i][nv];
      }
    }} 
  }
}
/* ********************************************************************* */
void PrimToCons3D(const Data *d, int where)
/*!
 *  Convert the data array d->Vc into d->Uc inside 
 *  the computational  domain.
 *********************************************************************** */
{
  int i,j,k,nv;
  static double **v;
  static unsigned char *flag;
  
  if (v == NULL) {
    v    = ARRAY_2D(NMAX_POINT,NVAR,double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);  
  }
  
  if (where == 0){
    KDOM_LOOP(k){
    JDOM_LOOP(j){
      IDOM_LOOP(i){
        for (nv = NVAR; nv--;  ) v[i][nv] = d->Vc[nv][k][j][i];
      }
      PrimToCons(v, d->Uc[k][j], IBEG, IEND);
    }} 
  }else if (where == 1){
    KTOT_LOOP(k){
    JTOT_LOOP(j){
      ITOT_LOOP(i){
        for (nv = NVAR; nv--;  ) v[i][nv] = d->Vc[nv][k][j][i];
      }
      PrimToCons(v, d->Uc[k][j], 0, NX1_TOT-1);
    }} 
  }
}
