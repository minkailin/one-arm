c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine gpbv
c     MKL comment: 
c     This is the original boundary potential solver that comes with ZEUSMP2. 
c     It assumes mass distribution is far from boundaries, which is inappropriate for disk
c     problems. So we cannot use this routine. So we calculate boundary potential with spherical
c     harmonics. 
c
c  Written by PSLi (12/4/99)
c
c  PURPOSE: Calculate the gravitational potential at the boundary
c           surfaces contributed from monopole and quadrupole moments
c           of a mass distribution.
c           Boundary potentals are calculated if the flags niis(3),
c           nois(3), nijs(3), nojs(3), niks(3), or noks(3) = 3,
c           respectively.
c
c  EXTERNALS: [none]
c
c  LOCALS:
c  dm(i,j,k)   is mass contained in zone i,j,k
c  lqm         is the local quadrupole moment (qm)
c  lcx,lcy,lcz coords of mass center at local processor (cx,cy,cz)
c  ltm         total mass at local processor (tm)
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use field
      use bndry
      use root
      use scratch
      use mpiyes
      use mpipar
c
      implicit NONE
c
      real(rl) ::     errmax
      integer  :: n
c
      real(rl) ::  dm(in,jn,kn),qm(6,nprocs_w),tm(nprocs_w),
     .             cx(nprocs_w),
     &             cy(nprocs_w),cz(nprocs_w),r1,dx,dy,dz,dx2,dy2,dz2,
     &             xdm,ydm,zdm,ltm,lcx,lcy,lcz,lqm(6),cthe,cphi,sthe,
     .             sphi
c
      integer  i,j,k,m,ism1,jsm1,ksm1,iep1,jep1,kep1
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////////
c=======================================================================
c  This routine calculate the monopole and quadrupole moments of a mass
c  distribution and the gravitational potential at the boundary surfaces
c  for solving Poisson equation.
c
c  Positive potential for the sign convention of ZeusMP.
c
      do i=1,6
        do m=1,nprocs_w
          qm(i,m) = 0.0
        enddo
        lqm(i) = 0.0
      enddo
      do k=1,kn
        do j=1,jn
          gpiib(j,k,1)=0.0
          gpoib(j,k,1)=0.0
        enddo
      enddo
      do k=1,kn
        do i=1,in
          gpijb(i,k,1)=0.0
          gpojb(i,k,1)=0.0
        enddo
      enddo
      do j=1,jn
        do i=1,in
          gpikb(i,j,1)=0.0
          gpokb(i,j,1)=0.0
        enddo
      enddo
      xdm=0.0
      ydm=0.0
      zdm=0.0
      ism1=is-1
      jsm1=js-1
      ksm1=ks-1
      iep1=ie+1
      jep1=je+1
      kep1=ke+1
c
c  Determine the center of mass.
c
      do i=1,nprocs_w
        tm(i)=tiny
        cx(i)=0.0
        cy(i)=0.0
        cz(i)=0.0
      enddo
      ltm=0.0
      lcx=0.0
      lcy=0.0
      lcz=0.0
c
      do k=ks,ke
       if(lgeom .eq. 2 .or. lgeom .eq. 3) then
        cphi=cos(x3b(k))
        sphi=sin(x3b(k))
       endif ! CYL or SPHERE
	do j=js,je
         if(lgeom .eq. 3) then
          cthe=cos(x2b(j))
          sthe=sin(x2b(j))
         endif ! SPHERE
          do i=is,ie
C            dm(i,j,k) = g2b(i)*g31b(i)*g32b(j)*dx1a(i)*
C     &                  dx2a(j)*dx3a(k)*d(i,j,k)
            dm(i,j,k) = d(i,j,k)*dvl1a(i)*dvl2a(j)*dvl3a(k)
            ltm=ltm+dm(i,j,k)
           if(lgeom .eq. 1) then
            xdm=xdm+x1b(i)*dm(i,j,k)
            ydm=ydm+x2b(j)*dm(i,j,k)
            zdm=zdm+x3b(k)*dm(i,j,k)
           endif ! CART
           if(lgeom .eq. 2) then
            xdm=xdm+x2b(j)*cphi*dm(i,j,k)
            ydm=ydm+x2b(j)*sphi*dm(i,j,k)
            zdm=zdm+x1b(i)*dm(i,j,k)
           endif ! CYL
           if(lgeom .eq. 3) then
            xdm=xdm+x1b(i)*sthe*cphi*dm(i,j,k)
            ydm=ydm+x1b(i)*sthe*sphi*dm(i,j,k)
            zdm=zdm+x1b(i)*cthe*dm(i,j,k)
           endif ! SPHERE
          enddo
        enddo
      enddo
c
      lcx=xdm/ltm
      lcy=ydm/ltm
      lcz=zdm/ltm
c
c  Compute the 2nd order terms in the expansion.
c
      do k=ks,ke
       if(lgeom .eq. 2 .or. lgeom .eq. 3) then
        sphi=sin(x3b(k))
        cphi=cos(x3b(k))
       endif ! CYL or SPHERE
	do j=js,je
         if(lgeom .eq. 3) then
          cthe=cos(x2b(j))
          sthe=sin(x2b(j))
         endif ! SPHERE
          do i=is,ie
           if(lgeom .eq. 1) then
            dx=x1b(i)-lcx
            dy=x2b(j)-lcy
            dz=x3b(k)-lcz
           endif ! CART
           if(lgeom .eq. 2) then
            dx=x2b(j)*cphi-lcx
            dy=x2b(j)*sphi-lcy
            dz=x1b(i)-lcz
           endif ! CYL
           if(lgeom .eq. 3) then
            dx=x1b(i)*sthe*cphi-lcx
            dy=x1b(i)*sthe*sphi-lcy
            dz=x1b(i)*cthe-lcz
           endif ! SPHERE
            dx2=dx*dx
            dy2=dy*dy
            dz2=dz*dz
            lqm(1)=lqm(1)+(2.*dx2-dy2-dz2)*dm(i,j,k)
            lqm(2)=lqm(2)+(2.*dy2-dx2-dz2)*dm(i,j,k)
            lqm(3)=lqm(3)+(2.*dz2-dx2-dy2)*dm(i,j,k)
            lqm(4)=lqm(4)+3.*dx*dy*dm(i,j,k)
            lqm(5)=lqm(5)+3.*dx*dz*dm(i,j,k)
            lqm(6)=lqm(6)+3.*dy*dz*dm(i,j,k)
          enddo
        enddo
      enddo
c
c  Sum up qm from all processors.
c
      call MPI_Gather(ltm,1,MPI_DOUBLE_PRECISION,
     &       tm,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Gather(lcx,1,MPI_DOUBLE_PRECISION,
     &       cx,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Gather(lcy,1,MPI_DOUBLE_PRECISION,
     &       cy,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Gather(lcz,1,MPI_DOUBLE_PRECISION,
     &       cz,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Gather(lqm,6,MPI_DOUBLE_PRECISION,
     &       qm,6,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(tm,nprocs_w,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(cx,nprocs_w,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(cy,nprocs_w,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(cz,nprocs_w,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
      call MPI_Bcast(qm,nprocs_w*6,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
c
c  Now compute the potential along each boundary surface, starting
c  with the inner i surface.
c
      if (niis(3) .eq. 3) then
        do m=1,nprocs_w
          do k=ks,ke
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            sphi=sin(x3b(k))
            cphi=cos(x3b(k))
           endif ! CYL OR SPHERE
            do j=js,je
             if(lgeom .eq. 3) then
              cthe=cos(x2b(j))
              sthe=sin(x2b(j))
             endif ! SPHERE
             if(lgeom .eq. 1) then
              dx=x1b(ism1)-cx(m)
              dy=x2b(j)-cy(m)
              dz=x3b(k)-cz(m)
             endif ! CART
             if(lgeom .eq. 2) then
              dx=x2b(j)*cphi-cx(m)
              dy=x2b(j)*sphi-cy(m)
              dz=x1b(ism1)-cz(m)
             endif ! CYL
             if(lgeom .eq. 3) then
              dx=x1b(ism1)*sthe*cphi-cx(m)
              dy=x1b(ism1)*sthe*sphi-cy(m)
              dz=x1b(ism1)*cthe-cz(m)
             endif ! SPHERE
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpiib(j,k,1) = gpiib(j,k,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+
     &                       qm(2,m)*dy2+qm(3,m)*dz2)+
     &                       qm(4,m)*dx*dy+qm(5,m)*dx*dz+
     &                       qm(6,m)*dy*dz)/r1**5
!              print*, (gpiib(j,k,1) - gp(is,j,k))/gp(is,j,k)
            enddo
          enddo
        enddo
      endif
c
c  Initial solution along outer i surface
c
      if (nois(3) .eq. 3) then
        do m=1,nprocs_w
          do k=ks,ke
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            sphi=sin(x3b(k))
            cphi=cos(x3b(k))
           endif
            do j=js,je
             if(lgeom .eq. 3) then
              cthe=cos(x2b(j))
              sthe=sin(x2b(j))
             endif
             if(lgeom .eq. 1) then
              dx=x1b(iep1)-cx(m)
              dy=x2b(j)-cy(m)
              dz=x3b(k)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(j)*cphi-cx(m)
              dy=x2b(j)*sphi-cy(m)
              dz=x1b(iep1)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(iep1)*sthe*cphi-cx(m)
              dy=x1b(iep1)*sthe*sphi-cy(m)
              dz=x1b(iep1)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpoib(j,k,1) = gpoib(j,k,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+
     &                       qm(2,m)*dy2+qm(3,m)*dz2)+
     &                       qm(4,m)*dx*dy+qm(5,m)*dx*dz+
     &                       qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
c
c  Initial solution along inner j surface
c
      if (nijs(3) .eq. 3) then
        do m=1,nprocs_w
         if(lgeom .eq. 3) then
          cthe=cos(x2b(jsm1))
          sthe=sin(x2b(jsm1))
         endif
          do k=ks,ke           
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            sphi=sin(x3b(k))
            cphi=cos(x3b(k))
           endif
            do i=is,ie
             if(lgeom .eq. 1) then
              dx=x1b(i)-cx(m)
              dy=x2b(jsm1)-cy(m)
              dz=x3b(k)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(jsm1)*cphi-cx(m)
              dy=x2b(jsm1)*sphi-cy(m)
              dz=x1b(i)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(i)*sthe*cphi-cx(m)
              dy=x1b(i)*sthe*sphi-cy(m)
              dz=x1b(i)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpijb(i,k,1) = gpijb(i,k,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+
     &                       qm(2,m)*dy2+qm(3,m)*dz2)+
     &                       qm(4,m)*dx*dy+qm(5,m)*dx*dz+
     &                       qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
c
c  Initial solution along outer j surface
c
      if (nojs(3) .eq. 3) then
        do m=1,nprocs_w
         if(lgeom .eq. 3) then
          cthe=cos(x2b(jep1))
          sthe=sin(x2b(jep1))
         endif
          do k=ks,ke           
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            sphi=sin(x3b(k))
            cphi=cos(x3b(k))
           endif
            do i=is,ie
             if(lgeom .eq. 1) then
              dx=x1b(i)-cx(m)
              dy=x2b(jep1)-cy(m)
              dz=x3b(k)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(jep1)*cphi-cx(m)
              dy=x2b(jep1)*sphi-cy(m)
              dz=x1b(i)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(i)*sthe*cphi-cx(m)
              dy=x1b(i)*sthe*sphi-cy(m)
              dz=x1b(i)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpojb(i,k,1) = gpojb(i,k,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+
     &                       qm(2,m)*dy2+qm(3,m)*dz2)+
     &                       qm(4,m)*dx*dy+qm(5,m)*dx*dz+
     &                       qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
c
c  Initial solution along inner k surface
c
      if (niks(3) .eq. 3) then
        do m=1,nprocs_w
         if(lgeom .eq. 3) then
          sphi=sin(x3b(ksm1))
          cphi=cos(x3b(ksm1))
         endif
          do j=js,je           
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            cthe=cos(x2b(j))
            sthe=sin(x2b(j))
           endif
            do i=is,ie
             if(lgeom .eq. 1) then
              dx=x1b(i)-cx(m)
              dy=x2b(j)-cy(m)
              dz=x3b(ksm1)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(j)*cphi-cx(m)
              dy=x2b(j)*sphi-cy(m)
              dz=x1b(i)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(i)*sthe*cphi-cx(m)
              dy=x1b(i)*sthe*sphi-cy(m)
              dz=x1b(i)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpikb(i,j,1) = gpikb(i,j,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+
     &                       qm(2,m)*dy2+qm(3,m)*dz2)+
     &                       qm(4,m)*dx*dy+qm(5,m)*dx*dz+
     &                       qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
c
c  Initial solution along outer k surface
c
      if (noks(3) .eq. 3) then
        do m=1,nprocs_w
         if(lgeom .eq. 3) then
          sphi=sin(x3b(kep1))
          cphi=cos(x3b(kep1))
         endif
          do j=js,je           
           if(lgeom .eq. 2 .or. lgeom .eq. 3) then
            cthe=cos(x2b(j))
            sthe=sin(x2b(j))
           endif
            do i=is,ie
             if(lgeom .eq. 1) then
              dx=x1b(i)-cx(m)
              dy=x2b(j)-cy(m)
              dz=x3b(kep1)-cz(m)
             endif
             if(lgeom .eq. 2) then
              dx=x2b(j)*cphi-cx(m)
              dy=x2b(j)*sphi-cy(m)
              dz=x1b(i)-cz(m)
             endif
             if(lgeom .eq. 3) then
              dx=x1b(i)*sthe*cphi-cx(m)
              dy=x1b(i)*sthe*sphi-cy(m)
              dz=x1b(i)*cthe-cz(m)
             endif
              dx2=dx*dx
              dy2=dy*dy
              dz2=dz*dz
              r1 = sqrt(dx2+dy2+dz2)
              gpokb(i,j,1) = gpokb(i,j,1)+tm(m)/r1+(0.5*(qm(1,m)*dx2+
     &                       qm(2,m)*dy2+qm(3,m)*dz2)+
     &                       qm(4,m)*dx*dy+qm(5,m)*dx*dz+
     &                       qm(6,m)*dy*dz)/r1**5
            enddo
          enddo
        enddo
      endif
c
      return
      end
c--------------------------------------------------------------------------------------------------------------------c
c                                                                                                                    c
c     Functions and routines for calculating the potential at radial boundaries using spherical harmonic expansion.  c  
c     Also some routines for self-gravity in general.                                                                c
c--------------------------------------------------------------------------------------------------------------------c
      real*8 FUNCTION plgndr(l,m,x)
      use param
      implicit none
      INTEGER :: l,m
      REAL*8 :: x
      INTEGER :: i,ll
      REAL*8 :: fact,oldfact,pll,pmm,pmmp1,omx2
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)then
         print*, 'bad arguments in plgndr'
         stop
      endif
      pmm=1d0
      if(m.gt.0) then
         omx2=(1d0-x)*(1d0+x)
         fact=1d0
         do i=1,m
            pmm=pmm*omx2*fact/(fact+1d0)
            fact=fact+2d0
         enddo
      endif
      pmm=sqrt((2*m+1d0)*pmm/(4d0*PI))
      if(mod(m,2).eq.1)pmm=-pmm
      if(l.eq.m) then
         plgndr=pmm
      else
         pmmp1=x*sqrt(2d0*m+3d0)*pmm
         if(l.eq.m+1) then
            plgndr=pmmp1
         else
            oldfact=sqrt(2d0*m+3d0)
            do ll=m+2,l
               fact=sqrt((4d0*ll**2-1d0)/(ll**2-m**2))
               pll=(x*pmmp1-pmm/oldfact)*fact
               oldfact=fact
               pmm=pmmp1
               pmmp1=pll
            enddo
            plgndr=pll
         endif
      endif
      return
      END
      complex*16 function sph_harm(l, m, polar, azimuth)
      implicit none
c
c     Compute the spherical harmonics 
c     = associated legendre polynomials * exp (i.m.phi)
c     
      integer :: l, m, absm
      real*8  :: polar, azimuth, x, mphi, pi ! theta, phi
      real*8, external :: plgndr
      pi = acos(-1d0)
      absm = iabs(m)            ! Get Y_lm for positive m first
      
      x = cos(polar)
      mphi = absm*azimuth 
      sph_harm = plgndr(l,absm,x)*dcmplx(cos(mphi), sin(mphi))
      if(m.lt.0) then
         sph_harm = (-1)**absm*conjg(sph_harm)
      endif
      return
      end
      subroutine get_sph_harm_array
      use gboundary
      use grid
c
c     Get spherical harmonic Ylm(l, m theta, phi) = Plm*exp(i.m.phi)
c     Note                   Yl(-m)= (-1)^m Ylm* , m>0 
c
c     We actually store it as Y_ml(m,l,theta,phi) though. 
c
      implicit none
      integer :: j, k, lcount, mcount, mrange, absm
      real*8  :: theta, phi, x
      real*8  :: plgndr
      complex*16, external :: sph_harm
      
      do k=ks-1, ke+1
         phi = x3b(k)
         do j=js-1, je+1
            theta = x2b(j)
            x = cos(theta)
            
            do lcount = lmin, lmax
               mrange = min(mmax, lcount)
               do mcount=-mrange, mrange
                  absm = iabs(mcount)
                  Y_ml(mcount, lcount, j, k) 
     &                 = sph_harm(lcount, mcount, theta, phi) ! This is the full spherical harmonic 
c
c     Below: we need Y_ml* to do integrate (rho.sin(theta).dtheta.dphi.Y_ml*)
c     Y_lm*    = Plm. exp (-i.m.phi)     for m > 0
c     Y_l(-m)* = (-1)^m.Plm.exp(i.m.phi) for m > 0 
c     or
c     Y_lm* = (-1)^|m|.Pl|m|.exp(-i.m.phi) for m < 0
c
c     So that in both cases (m>0, m<0) the exponential factor is exp(-i.m.phi)
c
                  Y_ml_star(mcount, lcount, j) ! This is is just the Legendre polynomial part of Y_ml* for POSITIVE m
     &                 = plgndr(lcount, absm, x)
                  if(mcount.lt.0) Y_ml_star(mcount,lcount,j) ! Insert factor of (-1)^|m| if m<0
     &                 =(-1)**(absm)*Y_ml_star(mcount,lcount,j) 
               enddo
            enddo
            
         enddo
      enddo
      print*, 'got the 4D spherical harmonic array'
      return
      end
      subroutine get_pot_ml(m, l)
c     
c     Get the pot_ml(r) for a given (m,l). 
c     Do the phi-theta integration locally. Then store it in appropriate subarray of the global rho_ml ( which 
c     covers the entire radial range). Then sum everyone's rho_ml, each CPU only contributes to this global sum 
c     at the appropriate radial cells. 
      
      use root
      use gboundary
      use grid
      use field
      use mpiyes
      use param
      use mpiyes
      use mpipar
      
      implicit none
      
      integer, parameter :: itmax=20
      real*8, parameter :: tol=1d-3, relax = 1.5
      integer :: i, j, k, m, l, ibeg, info, count,
     &     errest
      real*8  :: res, resmax, error
      real*8  :: thetaint(1:ntab_j),phiint(1:ntab_k),
     &     result, dphi, phiint_real(1:ntab_k), 
     &     phiint_imag(1:ntab_k), sintheta(1:ntab_j),
     &     result2, res_real, res_imag, dblem,
     &     ywork_j(3,ntab_j), ework_j(ntab_j), work_j(ntab_j),
     &     ywork_k(3,ntab_k), ework_k(ntab_k), work_k(ntab_k),
     &     hwork(2*(ntab_j-1)), hwork2(2*(ntab_k-1)),
     &     res_real2, res_imag2, corr
      complex*16 :: dd(0:nrad+1), du(0:nrad), dl(1:nrad+1), 
     &     trial(0:nrad+1), trialold(0:nrad+1), rhsinit(0:nrad+1)
      complex*16 :: contribution
      
      dblem = dble(m)
      sintheta(1:ntab_j) = sin(x2b(js:je))
      
c----------------------------------------------------------c
c     Step 1: get the rho_ml's by integrating              c
c                                                          c
c     rho.sin(theta).Ylm*.dtheta.dphi                      c
c                                                          c
c     over the local angular range. Then do global sum.    c  
c----------------------------------------------------------c
c     
c     Zero the local copy of the global rho_ml
c     
      rho_ml_1d_glob = dcmplx(0d0, 0d0)
      
c     
c     First do local phi-theta integration. Active density zones only.
c     Since we are only using this procedure to get BOUNDARY POTENTIAL,
c     we don't do accurate interpolations (which costs more) and simply do 
c     
c     sum (integrand.dx) 
c
c     I've tried those quadrature routines which uses interpolants but not sure
c     if that's a good idea because our problem will be parallelised, interpolation scheme
c     has no knowledge of other CPU's sub-grid. There may be issues at boundaries. 
      do i=is, ie               ! Go through each radial grid cell
         do k=ks, ke            ! Fix phi, do theta integration. No phi dependence here
            
            thetaint(1:ntab_j) = d(i,js:je,k)*Y_ml_star(m,l,js:je) ! This is the part of the integrand that depends on theta. 
            
            result = sum(thetaint(1:ntab_j)*dvl2a(js:je)) ! dvl2a = d(-cos theta) = sin(theta)*dtheta [Oct '12: changed to a-grid spacing]
                                                          ! Gives 0.011 error (relative difference between reconstructed density and actual density) 
            phiint(k - ks + 1)= result
         enddo
            
c     
c     Now do phi integration, do real part (cosine) and imaginary part( sine) separately. 
c     
         if(m.eq.0) then        ! No phi dependence. 
            result2 = sum(phiint(1:ntab_k)*dx3a(ks:ke)) ![Oct '12:changed to a-grid spacing]
            contribution = dcmplx(result2)            
         else                   ! There is phi dependence. Insert factor exp(-i.m.phi) = cos(m.phi) - i.sin(m.phi)
               
c          Nov '12 correctoin factor for use later, ASSUME UNIFORM PHI
c          SPACING. note there are no end-point corrections because periodicity
c          in phi
       
            corr = dblem*dx3a(ks)/2d0     
            corr = (sin(corr)/corr)**2
            phiint_real(1:ntab_k) =  phiint(1:ntab_k)
     &           *cos(dblem*phiaxis(1:ntab_k))
            phiint_imag(1:ntab_k) = -phiint(1:ntab_k)
     &           *sin(dblem*phiaxis(1:ntab_k))
            
            res_real = sum(phiint_real(1:ntab_k)*dx3a(ks:ke)) ! Actually looks good! [Oct '12: changed to a-grid spacing]
            res_imag = sum(phiint_imag(1:ntab_k)*dx3a(ks:ke)) ! [Oct '12: changed to a-grid spacing]
c           Nov '12 correction factors [assume function constant in cell, do
c           trig integral by hand]
            
            res_real = res_real*corr
            res_imag = res_imag*corr
            contribution = dcmplx(res_real,res_imag)  
         endif
         if(sgzsymm.eq..true.) contribution = contribution*2d0 ! If assume symmetric about z=0, need to double up to account for lower half disk. 
         rho_ml_1d_loc(i-is+1) = contribution*4d0*pi*x1b(i)**2 ! Effective RHS of 1D Poisson in log grid. Assumed G = 1 here. 
      enddo
      
c     
c     Put the local 1D array above into the correct subarray of the global version.
c     Only the appropriate subarray is filled, radii corresponding to other CPUs are left zero. 
c      
      ibeg = coords(1)*nrad_loc + 1
      rho_ml_1d_glob(ibeg:ibeg+nrad_loc-1) = rho_ml_1d_loc(1:nrad_loc)
c     
c     Now do global sum. rho_ml is the receiving array. 
c     
      call mpi_allreduce(rho_ml_1d_glob, rho_ml, nrad, 
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, comm3d,ierr)
c    Save the rho_lm for diagonstic checks if desired. ***DISABLED***
c
c      rho_ml_3d(m,l,is:ie) = rho_ml(ibeg:ibeg+nrad_loc-1)
c     &     /(4d0*pi*x1b(is:ie)**2)
c----------------------------------------------------------c
c                                                          c
c     Step 2: matrix problem                               c
c     Everyone has the global distribution of rho_ml(r)    c
c     Everyone solves the global problem. Since this is 1D,c
c     shouldn't cost too much.                             c
c----------------------------------------------------------c
c     
c     Construct RHS interior points first. This result from active density cells. 
c     
      rhs(1:nrad) = rho_ml(1:nrad)
c     
c     Boundary condition for zero mass outside boundaries or zero gradient in potential
c
       rhs(0)     = 0d0 
       rhs(nrad+1)= 0d0 
c     
c     Construct elements of the tridiagonal matrix. Note that zgtsv overwrites these 1D arrays so they are
c     no good for the next l
c
c     Interior points, for active zones
      du(1:nrad) = udiag(1:nrad)
      dd(1:nrad) = -2d0*dlogri**2 - l*(l+1d0)
      dl(1:nrad) = ldiag(1:nrad) 
c     IBC: du(0) and dd(0)
c     For zero mass inside boundary:
      du(0)      = dlogri - l/2d0
      dd(0)      = -(dlogri + l/2d0)
c     OBC: dd(nrad+1) and dl(nrad+1) controls OBC
c     For zero mass outside boundary
      dd(nrad+1) = dlogri + (l+1d0)/2d0
      dl(nrad+1) = -dlogri + (l+1d0)/2d0 
c     
c     Solve matrix problem
c
      call zgtsv(nrad+2, 1, dl, dd, du, rhs, nrad+2, info)
      if(info.ne.0) then
         print*, 'did not solve poisson successfully, info=', info
         stop
      endif
      
c
c     rhs is now the phi_ml(r) 1D vector that covers the global disk. 
c     Note that phi_ml has size nrad + 2. phi_ml(0)      is potential at inner radial ghost cell
c                                         phi_ml(nrad+1) is potential at outer radial ghost cell 
c     Copy appropriate subsection into local array for potential calculation. Including boundary/ghost cells. 
      pot_ml_3d_loc(m,l,is-1:ie+1) = rhs(ibeg-1:ibeg+nrad_loc)
c     Save a copy of the global solution.
c     pot_ml_3d_glob(m,l,0:nrad+1) = rhs(0:nrad+1)
c
      return
      end
      subroutine get_potBC
c
c     Given the local pot_ml for all (l,m) of interest, get the potential at r and theta BOUNDARIES
c     Sum is performed over (l,m)
c
c     Solution method is for standard convention (negative potential) 
c     but ZEUS convection is positive potential, so add minus sign
c
      use bndry
      use grid
      use gboundary
      use field
      implicit none
      integer :: i, j, k, l, m
      integer :: mrange
      complex*16 :: rgpin, rgpout, tgpin, tgpout
      complex*16 :: contribution      
c
c     Everyone calculates their boundary/ghost cell potential. But only 
c     the CPUs that are responsible for physical boundaries will use the result. 
      do k=ks-1, ke+1           ! Go through each azimuth
c
c     Radial boundaries, at is-1 and ie+1, do all j
c
         do j=js-1, je+1
            
c            rgpin = sum(pot_ml_3d_loc(-mmax:mmax, lmin:lmax, is-1)
c     &           *Y_ml(-mmax:mmax, lmin:lmax, j, k))             
c            
c            rgpout= sum(pot_ml_3d_loc(-mmax:mmax,lmin:lmax,  ie+1)
c     &           *Y_ml(-mmax:mmax,lmin:lmax , j, k))
            
c   Nov' 12. only have coefficients for positive m's account for negative m's
c   using the fact that negative coefficients are related to positive ones
            
             rgpin = sum(pot_ml_3d_loc(0, lmin:lmax, is-1)
     &           *Y_ml(0, lmin:lmax, j, k))
             rgpin = rgpin 
     &       + 2d0*sum(pot_ml_3d_loc(1:mmax, 1:lmax, is-1)
     &                 *Y_ml(1:mmax, 1:lmax, j, k))            
             rgpout= sum(pot_ml_3d_loc(0,lmin:lmax,  ie+1)
     &           *Y_ml(0,lmin:lmax , j, k))
             rgpout = rgpout 
     &       + 2d0*sum(pot_ml_3d_loc(1:mmax,1:lmax,ie+1)
     &                 *Y_ml(1:mmax,1:lmax , j, k))
            gpiib(j,k,1) = -dble(rgpin)
            gpoib(j,k,1) = -dble(rgpout)
         enddo
c     
c     theta boundaries, at js-1 and je+1, do all i
c
         do i=is-1, ie+1
c            tgpin = sum(pot_ml_3d_loc(-mmax:mmax, lmin:lmax, i)
c     &           *Y_ml(-mmax:mmax, lmin:lmax, js-1, k))             
c            
c            tgpout= sum(pot_ml_3d_loc(-mmax:mmax,lmin:lmax,  i)
c     &           *Y_ml(-mmax:mmax,lmin:lmax , je+1, k))
           
c   Nov' 12. only have coefficients for positive m's account for negative m's
c   using the fact that negative coefficients are related to positive ones
            tgpin = sum( pot_ml_3d_loc(0, lmin:lmax, i)
     &                  *Y_ml(0, lmin:lmax, js-1, k))
            tgpin =tgpin 
     &      + 2d0*sum(pot_ml_3d_loc(1:mmax, 1:lmax, i)
     &                     *Y_ml(1:mmax, 1:lmax, js-1, k))                      
 
            tgpout= sum( pot_ml_3d_loc(0,lmin:lmax,  i)
     &                  *Y_ml(0,lmin:lmax , je+1, k))
            tgpout= tgpout 
     &      + 2d0*sum(pot_ml_3d_loc(1:mmax,1:lmax,i)
     &                *Y_ml(1:mmax,1:lmax , je+1, k))
            gpijb(i,k,1) = -dble(tgpin)
            gpojb(i,k,1) = -dble(tgpout)
         enddo
 
      enddo
      return
      end
     
    
      real*8 function sgacc_cut(tt)
c     Tapering SG acceleration. If vertical structure not setup to include SG then
c     switch on theta-SG acceleration slowly. To switch it off, put sg_on to be large negative,
c     so full theta acceleration applied from start 
      use gboundary
      use param
      implicit none
      real*8 :: tt, temp
      if(tt .lt. tsg) then
         sgacc_cut = 0d0
      else if((tt .ge. tsg).and.
     &        (tt .lt. (tsg + sg_on))) then
         temp = (tt - tsg)/(sg_on)
         sgacc_cut = (sin(temp*pi/2d0))**2
      else 
         sgacc_cut = 1d0
      endif
      return
      end        
  
      subroutine sgBCboss
c
c     Driver routine to get potential at boundaries. 
c
      use gboundary
      use mpiyes
      use mpipar
      use field
      use bndry
      implicit none
      integer :: l, m
      integer :: mrange
      pot_ml_3d_loc = dcmplx(0d0, 0d0) ! Initialize
      do l=lmin, lmax
         mrange = min(mmax, l)
c         do m=-mrange, mrange
          do m=0, mrange ![Nov '12: just get positive m]
            if(mod(l+m,2).eq.0) call get_pot_ml(m, l) ! Only do l+m = even modes
         enddo
      enddo
      
      call mpi_barrier(comm3d, ierr)
      call get_potBC    
      return
      end
