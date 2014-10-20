;+
;
; NAME:      H5LOAD
;
; AUTHOR:    C. Zanni
;
; REQUIRES:  HDF5 support.
;
; PURPOSE:   Read PLUTO variables stored in a (pixie) HDF5 data file 
;            (static grid).  It is automatically called by PLOAD when the 
;            /H5 keyword is supplied.
;
; SYNTAX:    H5LOAD, filename
;
; ARGUMENTS
;
;   filename   = the data file name
;
; KEYWORDS 
;
;   /SILENT    Set this keyword to suppress output
;
; LAST MODIFIED
;
;   Sept 25, 2012 by C. Zanni (zanni@oato.inaf.it)
;
;-
PRO H5LOAD, filename, silent=silent

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

 checkfile,filename
 ifile = H5F_OPEN(filename)
 
 igrp  = H5G_OPEN(ifile,'/vars')

 numvars = H5G_GET_NMEMBERS(ifile,"/vars")

 FOR nv=0,numvars-1 DO BEGIN
  vname = H5G_GET_MEMBER_NAME(ifile,"/vars",nv)
  idata = H5D_OPEN(igrp,vname)
  vpt = H5D_READ(idata)
  MATCH_VARNAME, vpt, vname, silent=silent
  H5D_CLOSE,idata 
 ENDFOR

 H5G_CLOSE,igrp

 numgrp = H5G_GET_NMEMBERS(ifile,"/")
 IF (numgrp gt 3) THEN BEGIN

   igrp  = H5G_OPEN(ifile,'/stag_vars')
   numvars = H5G_GET_NMEMBERS(ifile,"/stag_vars")

   FOR nv=0,numvars-1 DO BEGIN
     vname = H5G_GET_MEMBER_NAME(ifile,"/stag_vars",nv)
     idata = H5D_OPEN(igrp,vname)
     vpt = H5D_READ(idata)
     MATCH_VARNAME, vpt, vname, silent=silent
     H5D_CLOSE,idata 
   ENDFOR
   H5G_CLOSE,igrp
 ENDIF

 H5F_CLOSE,ifile
 RETURN
END
