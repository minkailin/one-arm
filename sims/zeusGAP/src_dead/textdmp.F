c=======================================================================
c/////////////////////////  EXAMPLE USERDUMP  \\\\\\\\\\\\\\\\\\\\\\\\\c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine textdmp
c
c  written by:   Robert Fiedler
c
c  PURPOSE:  Example USERDUMP routine to generate ASCII output
c            files during a ZEUS-MP run.
c
c  USAGE:  To use this routine, include "-u textdmp" on the
c          Make_zeusmp command line.
c
c  DESCRIPTION:  This example USERDUMP routine prints out the
c                coordinates and the density, energy density,
c                and velocity components for a 3-D
c                hydrodynamics simulation.  It must be modified
c                to list magnetic field components, 
c                gravitational potential, radiation energy density, 
c                etc.  
c
c-----------------------------------------------------------------------
      use real_prec
      use param
      use config
      use field
      use grid
      use bndry
      use root
      use cons
      use scratch
      use opac_law
      use mpiyes
      use mpipar
c
      implicit NONE
c
      integer :: i, j, k
      real(rl) :: etotal, ekin, eint, etot_glb, dvb, totdvb,
     .            glb_j, jtot, glb_ek(3),
     .            ek(3), v1scr(in), v2scr(in), b2scr(in), b3scr(in),
     .            v3scr(in), b1scr(in), dscr(in), pscr(in)
c
      real(rl) :: so_eps , tau    , eriibn, timarg
c
      common /soui/ so_eps, tau, eriibn
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c Write out the coordinates and the density for each zone in the
c computational domain to unit 12 (usrfile).
c
12    format(1p5e16.8)
13    format(1p6e16.8)
14    format(1p7e16.8)
15    format(1p8e16.8)
      open(12,file=usrfile)
c
c-----------------------------------------------------------------------
c     if MHD, then check to verify that div(B) = 0 everywhere
c-----------------------------------------------------------------------
c
      if(xmhd) then
       call diverg ( b1, b2, b3, 1, 1, w3da, dvb )
       call MPI_Reduce(dvb, totdvb, 1,
     &                 MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm3d, ierr)
      endif ! xmhd
c
      if(coords(1) .eq. 0) then
       write(12,"('ZEUS-MP field variables at time = ',1pe15.8)") time
       write(12,"(a)")
      endif
c
      if(lrad .eq. 0) then
       if(xmhd) then
        write(12,"('Total div(B) = ',1pd12.4/)")totdvb
       endif ! xmhd
       write(12,"('         x1b             x2b             x3b',
     1            '         density         e               v1 ')")
       if(xtotnrg .eqv. .false.) then
        write(12,13) (((x1b(i), x2b(j), x3b(k),
     .                  d(i,j,k), 
     .                  e(i,j,k),
C     .                  gamm1*mmw*mh*e(i,j,k)/(boltz*d(i,j,k)),
     .                  v1(i,j,k),
     .                  i=is,ie), j=js,je ), k=ks,ke)
       else
        write(12,13) (((x1b(i), x2b(j), x3b(k),
     .                  d(i,j,k), 
     .                  e(i,j,k) - 0.125D0*d(i,j,k)*(
     .                      (v1(i,j,k)+v1(i+1,j,k))**2 
     .                    + (v2(i,j,k)+v2(i,j+1,k))**2 
     .                    + (v3(i,j,k)+v3(i,j,k+1))**2        
     .                                               ),
     .                  v1(i,j,k),
     .                  i=is,ie), j=js,je), k=ks,ke)
       endif
      else ! lrad
       write(12,"('         x1b             x2b             x3b',
     1            '         density         Tgas            Trad ')")
       write(12,13) (((x1b(i), x2b(j), x3b(k), d(i,j,k),
     .                 gamm1*mmw*mh*e(i,j,k)/(everg*d(i,j,k)),
     .                 boltz*sqrt(sqrt(er(i,j,k)/rad_con))/everg,
     .                 i=is,ie), j=js,je ), k=ks,ke)
      endif ! lrad
c
22    format(1p4e16.8)
23    format(1p5e16.8)
      write(12,"(a)")
      write(12,"('                     x1b             x2a        ',
     1 '     x3b          v2')")
      write(12,32) ( ( (i,j,k,x1b(i), x2a(j), x3b(k), v2(i,j,k),
     1                 i=is,ie), j=js,je ), k=ks,ke )
      write(12,"(a)")
      write(12,"('         x1b             x2b             x3a',
     1 '          v3')")
      write(12,42) ( ( ( x1b(i), x2b(j), x3a(k), 
     .                   v3(i,j,k),
     1                 i=is,ie), j=js,je ), k=ks,ke )
c
32    format(3i4,1p4e16.8)
42    format(1p4e16.8)
52    format(1p5e16.8)
c
      close(12)
c
      return
      end
