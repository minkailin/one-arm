
c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine hdfall(filename)
c
c  PURPOSE: Makes an hdf dump of all the active field variables.  The
c  set of field variables dumped is problem specific (depends on what
c  physics is defined).  Data is written in the Scientific Data Set
c  format to the file zhzXXNNNNNN.MMM.
c  Note that data must be stored column major and contiguously in order
c  to interface correctly to the C hdf routines.  All variables are
c  dumped as zone centered quantities.
c
c  EXTERNALS: HDF library routines
c
c  LOCALS:
c
c  LAST MODIFIED: by JCH; 3/12/97.
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use field
      use root
      use scratch
      use cons
      use mpiyes
      use mpipar
      use gboundary
c
      implicit NONE
c
      character*15 :: filename
      character*16 :: coordsys
      character*32 :: string
c
      integer      :: i,j,k,indx,kp1
c
c real is correct on PVP Crays, but real*4 is required on T3E.
c This is due to the kind of C floats on the two machines according to
c Albert Cheng of the NCSA HDF group.  The varying integer definitions
c I have decided on by trial and error M-MML 17.5.98
c
C---------------------------------------------------------------------
C      use on non-T3E UNICOS systems
C
C      integer rank,shape(3),ret
C      real data(in*jn*kn),xscale(in),yscale(jn),zscale(kn)
C      integer  dssdims,dssdast,dssdisc,dsadata,dspdata
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C     use on T3E
C
C      integer rank,shape(3),ret
C      real*4 data(in*jn*kn),xscale(in),yscale(jn),zscale(kn)
C      integer  dssdims,dssdast,dssdisc,dsadata,dspdata
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C     use on everything else
C
      integer(kind=4) :: rank,shape(3),ret
      real(rl4)       :: data(in*jn*kn),xscale(in),yscale(jn),
     .                   zscale(kn)
      integer(kind=4) :: dssdims,dssdast,dssdisc,dsadata,dspdata
C---------------------------------------------------------------------
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c
      if(myid_w .eq. 0) then
         print*,'-----------------------------------------'
         print*,'OUTPUT   ',filename
         print*,'-----------------------------------------'
      endif
      if(lgeom .eq. 3) coordsys = 'spherical polar' // char(0)
      if(lgeom .eq. 2) coordsys = 'cylindrical' // char(0)
      if(lgeom .eq. 1) coordsys = 'cartesian' // char(0)
c
      do 10 i=is,ie
        xscale(i-is+1) = real(x1b(i))
10    continue
      do 20 j=js,je
        yscale(j-js+1) = real(x2b(j))
20    continue
      do 30 k=ks,ke
        zscale(k-ks+1) = real(x3b(k))
30    continue
c
      rank     = 3
      shape(1) = nx1z
      shape(2) = nx2z
      shape(3) = nx3z
      ret = dssdims(rank,shape)
      ret = dssdisc(1,shape(1),xscale)
      ret = dssdisc(2,shape(2),yscale)
      ret = dssdisc(3,shape(3),zscale)
c     
c  1-velocity
c
      do 120 k=ks,ke
        do 110 j=js,je
          do 100 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
	    data(indx) = real(0.5*(v1(i,j,k) + v1(i+1,j,k)))
100       continue
110     continue
120   continue
      write(string,"('1-VELOCITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dspdata(filename,rank,shape,data)
c
c  2-velocity
c
      do 220 k=ks,ke
        do 210 j=js,je
          do 200 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
 	    data(indx)  = real(0.5*(v2(i,j,k) + v2(i,j+1,k)))
200       continue
210     continue
220   continue
      write(string,"('2-VELOCITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c
c  3-velocity
c
      do 320 k=ks,ke
        if(ldimen .eq. 3) then
         kp1 = k+1
        else
         kp1 = k
        endif
        do 310 j=js,je
          do 300 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
	    data(indx)  = real(0.5*(v3(i,j,k) + v3(i,j,kp1)))
300       continue
310     continue
320   continue
      write(string,"('3-VELOCITY AT TIME=',1pe8.2,'     ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c
c     this output will  be used for gravitational potential in the future. 
c
c  
c
      do 420 k=ks,ke
        do 410 j=js,je
          do 400 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
             if(xgrav.eq..true.) then
                data(indx)  = -real(gp(i,j,k))
!                data(indx)  = rho_recons(i,j,k)
             else
                data(indx) = 0d0
             endif
400       continue
410     continue
420   continue
      write(string,"('1-MAG FIELD AT TIME=',1pe8.2,'    ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c      if(myid_w .eq. 0) print*, 'magnetic output', ret
      if(xmhd) then
c
c  2-magnetic field
c
      do 520 k=ks,ke
        do 510 j=js,je
          do 500 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
	    data(indx)  = real(0.5*(b2(i,j,k) + b2(i,j+1,k)))
500       continue
510     continue
520   continue
      write(string,"('2-MAG FIELD AT TIME=',1pe8.2,'    ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c
c  3-magnetic field
c
      do 620 k=ks,ke
        if(ldimen .eq. 3) then
         kp1 = k+1
        else
         kp1 = k
        endif
        do 610 j=js,je
          do 600 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
	    data(indx)  = real(0.5*(b3(i,j,k) + b3(i,j,kp1)))
600       continue
610     continue
620   continue
      write(string,"('3-MAG FIELD AT TIME=',1pe8.2,'    ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
      endif ! xmhd
c
c  density
c
      do 720 k=ks,ke
        do 710 j=js,je
          do 700 i=is,ie
             indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
             data(indx) = real(d(i,j,k))
700      continue
710     continue
720   continue
      write(string,"('DENSITY AT TIME=',1pe8.2,'        ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c      if(myid_w .eq. 0) print*, 'density output',ret
c
c  internal energy
c
      do 820 k=ks,ke
        do 810 j=js,je
          do 800 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
            data(indx) = real(e(i,j,k))
800       continue
810     continue
820   continue
      write(string,"('GAS ENERGY AT TIME=',1pe8.2,'   ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
c      if(myid_w .eq. 0) print*, 'energy output',ret
c
c  radiation internal energy
c
      if(lrad .ne. 0) then
      do 1020 k=ks,ke
        do 1010 j=js,je
          do 1000 i=is,ie
            indx = (k-ks)*nx2z*nx1z + (j-js)*nx1z + (i-is) + 1
            data(indx) = real(er(i,j,k))
1000      continue
1010    continue
1020  continue
      write(string,"('RADIATION T AT TIME=',1pe8.2,'    ')") time
      ret = dssdast(string,' ',' ',coordsys)
      ret = dsadata(filename,rank,shape,data)
      endif ! lrad
c
      return
      end
