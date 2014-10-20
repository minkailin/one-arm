c$$$      subroutine pot_den_iteration
c$$$!
c$$$!     iteration scheme to get potential-density pair consistent with
c$$$!     vertical structure
c$$$!     solve for Phi and rho according to:
c$$$!     
c$$$!     0 = c^2 dlog rho/dz - dPhi_*/dz + dPhi/dz
c$$$!     we actually solve for beta, such that rho = beta*rho_NSG 
c$$$!     because rho_NSG was constructed to balance central star. 
c$$$!     the z derivatives are defined at j-faces (theta faces)
c$$$!     positive sign with SG potenital is just zeus convention
c$$$
c$$$      use field
c$$$      use grid
c$$$      use mpipar
c$$$      use mpiyes
c$$$      use planet
c$$$      use bndry
c$$$      
c$$$      implicit none
c$$$
c$$$      integer, parameter :: itmax=200, itmin=50
c$$$      integer :: bvstat_save(8,nbvar) !we're going to use boundary communication for density field
c$$$                                       !which over rides bvstat, use dummy to restore to value before
c$$$                                       !calling this routine
c$$$      real*8, parameter :: tolerence = 1d-4, relax = 1.0
c$$$      integer :: i,j,k, count, ii
c$$$      integer*4 :: timeArray(3)
c$$$      real :: rand 
c$$$      real*8 :: dold(in,jn,kn), beta(in,jn,kn), logbeta(in,jn,kn),
c$$$     . beta_old(in,jn,kn)
c$$$      real*8 :: resmax_local, resmax_global, dnew, dnsg,
c$$$     .     rad, theta, csq, bigR, mass, discmass,
c$$$     .     dphidr, dphidtheta, dphidz, discmass_old,
c$$$     .     dlogbetadr, logbeta_new, temp, taper,
c$$$     .     dphistar_dz, zz
c$$$      real*8, external :: initialdens
c$$$
c$$$      bvstat_save = bvstat
c$$$      if(myid.eq.0) then
c$$$         print*, '-------------------------------------------------'
c$$$         print*, 'parallel SOR to get vertical structure with SG...'
c$$$         print*, '-------------------------------------------------'
c$$$      endif
c$$$      
c$$$      beta = d
c$$$!
c$$$!     iteration scheme
c$$$!
c$$$      count = 0
c$$$      resmax_global = 1d0
c$$$      resmax_local  = 1d0
c$$$      discmass      = 1d0
c$$$
c$$$      do while((resmax_global.gt.tolerence).or.
c$$$     &     (count.lt.itmin))
c$$$         logbeta = log(d)
c$$$         beta_old = d
c$$$         discmass_old = discmass
c$$$!
c$$$!     get potential of old density field
c$$$!
c$$$         call gravity
c$$$!
c$$$!     communcate boundaries, need i, j BC
c$$$!            
c$$$         nreq = 0 
c$$$         nsub = nsub + 1
c$$$         call bvald(1,1,1,1,0,0,logbeta)
c$$$         if(nreq.ne.0) call MPI_WAITALL(nreq, req, stat, ierr)
c$$$
c$$$         mass = 0d0
c$$$         do k=ks, ke
c$$$            do j=js, je
c$$$               theta = x2a(j)
c$$$               do i=is, ie
c$$$                  rad = x1b(i)
c$$$!
c$$$!     work out quantities at the j-face(e.g. a-mesh for theta above)
c$$$!
c$$$                  bigR= rad*sin(theta)
c$$$                  zz  = rad*cos(theta)
c$$$                  csq = smallh**2/bigR
c$$$                  
c$$$                  dphistar_dz = -zz/rad**3
c$$$
c$$$                  dphidtheta  = (gp(i,j,k)-gp(i,j-1,k))*dx2bi(j-1)
c$$$                  dphidr = (0.5*(gp(i+1,j,k)+gp(i+1,j-1,k)) - 
c$$$     &                 0.5*(gp(i-1,j,k)+gp(i-1,j-1,k)))
c$$$     &                 /(x1b(i+1) - x1b(i-1))
c$$$
c$$$                  dphidz = 0.0 !cos(theta)*dphidr - sin(theta)*dphidtheta/rad
c$$$                  dphidz = (dphistar_dz + dphidz)/csq !this maybe very large
c$$$
c$$$
c$$$
c$$$
c$$$!
c$$$!     do not have full SG from the start, taper it
c$$$!
c$$$!                  if(count .lt. itmin) then
c$$$!                     temp = dble(count)/dble(itmin)
c$$$!                     taper = (sin(temp*pi/2d0))**2
c$$$!                  else
c$$$!                     taper = 1d0
c$$$!                  endif
c$$$!                  dphidz = dphidz*taper
c$$$
c$$$
c$$$                  dlogbetadr=(0.5*(logbeta(i+1,j,k)+logbeta(i+1,j-1,k))- 
c$$$     &                 0.5*(logbeta(i-1,j,k)+logbeta(i-1,j-1,k)))
c$$$     &                 /(x1b(i+1) - x1b(i-1))
c$$$                  dlogbetadr= dlogbetadr*cos(theta)
c$$$                  
c$$$                  logbeta_new = dphidz - dlogbetadr
c$$$                  
c$$$
c$$$
c$$$                  logbeta_new =-logbeta_new*rad*dx2b(j-1)/sin(theta)
c$$$                  logbeta_new = logbeta_new + logbeta(i,j-1,k)
c$$$                  logbeta_new = (1d0 - relax)*logbeta(i,j,k) 
c$$$     &                 +relax*logbeta_new
c$$$                  
c$$$
c$$$
c$$$                  beta(i,j,k) = exp(logbeta_new)
c$$$!                  if((beta(i,j,k).lt.0.8).or.
c$$$!     &                 (beta(i,j,k).gt.1.2)) then
c$$$!                     beta(i,j,k) = beta_old(i,j,k)
c$$$!                  endif
c$$$
c$$$                  d(i,j,k) = beta(i,j,k)
c$$$
c$$$                  mass = mass + d(i,j,k)
c$$$     &                 *rad**2*sin(theta)*dx1a(i)
c$$$     &                 *dx2a(j)*dx3a(k)
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$!
c$$$!     we updated interior points. so boundaries are outdated
c$$$!
c$$$
c$$$         bvstat = 0
c$$$
c$$$         resmax_local = maxval(abs(
c$$$     &        (beta(is:ie,js:je,ks:ke) - beta_old(is:ie,js:je,ks:ke)) 
c$$$     &        /beta_old(is:ie,js:je,ks:ke)))
c$$$!         print*, 'resmax and mass local=', resmax_local, mass
c$$$         call MPI_allREDUCE(resmax_local, 
c$$$     &        resmax_global, 1, MPI_DOUBLE_PRECISION,
c$$$     &        MPI_MAX, comm3d, ierr)
c$$$         call MPI_allREDUCE(mass, 
c$$$     &        discmass, 1, MPI_DOUBLE_PRECISION,
c$$$     &        MPI_SUM, comm3d, ierr)
c$$$         
c$$$!         if(discmass.gt.discmass_old) then 
c$$$!            d = d*(discmass_old/discmass)
c$$$!         endif
c$$$         if(myid.eq.0) then
c$$$            write(*,"('iter, resmax, mdisc, den0 = ',i4,3e22.15)")
c$$$     .              count, resmax_global, discmass, den0
c$$$            if(discmass.gt.0.1) then
c$$$               print*, 'disc too massive!'
c$$$               stop
c$$$            endif
c$$$         endif
c$$$
c$$$         count = count + 1
c$$$         if(count.gt.itmax) then
c$$$            print*, 'could not find sg eqm solution'
c$$$            stop
c$$$         endif
c$$$      enddo
c$$$      call MPI_barrier(comm3d, ierr)
c$$$!
c$$$!     restore bvstat
c$$$!
c$$$      bvstat = bvstat_save
c$$$
c$$$!      if(myid.eq.0) print*,'SG vertical: new mass =', 
c$$$!     &     discmass
c$$$      return
c$$$      end


c$$$        i=2  
c$$$        do  k=ks-1, ke+1
c$$$           do  j=js-1, je+1
c$$$              theta = x2b(j)
c$$$
c$$$              rad   = x1b(ie + i - 1)
c$$$              radp1 = x1b(ie + i)
c$$$              rho0    = initialdens(rad, theta, sig, den0,
c$$$     &                   innerhole)
c$$$              vphi0 = sqrt(subkepvel(rad,theta,sig,innerhole))
c$$$              cs0 = sqrt(smallh**2/(rad*sin(theta)))
c$$$              cs1 = sqrt(smallh**2/(radp1*sin(theta)))
c$$$
c$$$              dW  = godon_doib(j,k,1) - rho0
c$$$              dU =  godon_doib(j,k,1)
c$$$     &              *(godon_v1oib(j,k,1) + godon_v1oib(j,k,2))/2d0
c$$$
c$$$              coeff1 = (cs0*dW-dU)/(2d0*cs0)
c$$$              coeff3 = (cs0*dW+dU)/(2d0*cs0)
c$$$
c$$$              if((godon_v1oib(j,k,1)+cs0).lt. 0.0) coeff3 = 0d0
c$$$              if((godon_v1oib(j,k,1)-cs0).lt. 0.0) coeff1 = 0d0
c$$$
c$$$              godon_doib(j,k,i) = rho0 + (coeff1 +coeff3)
c$$$              godon_v1oib(j,k,i)= (-cs1*coeff1 + cs1*coeff3)
c$$$     &             /(rho0 + (coeff1 +coeff3))
c$$$
c$$$              dW = 0.5*(godon_doib(j,k,1)+godon_doib(j,k-1,1)) - rho0
c$$$              dU = 0.5*(godon_doib(j,k,1)+godon_doib(j,k-1,1))
c$$$     &             *(godon_v1oib(j,k,1) + godon_v1oib(j,k,2)
c$$$     &                    +godon_v1oib(j,k-1,1))/3d0
c$$$              dV = godon_v3oib(j,k,1)*godon_doib(j,k-1,1)-rho0*vphi0
c$$$
c$$$              coeff1 = (cs0*dW-dU)/(2d0*cs0)
c$$$              coeff2 = -vphi0*dW+dV
c$$$              coeff3 = (cs0*dW+dU)/(2d0*cs0)
c$$$
c$$$              if((godon_v1oib(j,k,1)+cs0).lt.0.0) coeff3 = 0d0
c$$$              if(godon_v1oib(j,k,1).lt. 0.0)      coeff2 = 0d0
c$$$              if((godon_v1oib(j,k,1)-cs0).lt. 0.0)coeff1 = 0d0
c$$$
c$$$
c$$$              godon_v3oib(j,k,i)=(rho0*vphi0+vphi0*coeff1+vphi0*coeff3+
c$$$     &             coeff2)/(rho0 + (coeff1 +coeff3))
c$$$
c$$$           enddo
c$$$        enddo


c$$$      dphi = 2d0*pi/128d0
c$$$      dtheta = pi/64d0
c$$$
c$$$      do j=1, 64
c$$$         theta(j) = dtheta*j
c$$$      enddo
c$$$      do k=1, 128
c$$$         azi(k) = dphi*k
c$$$      enddo
c$$$
c$$$      l = 4
c$$$      m = 2
c$$$      
c$$$      res = dcmplx(0d0, 0d0)
c$$$      do k=1, 128
c$$$         do j=1, 64
c$$$            res = res +
c$$$     &           sph_harm(l,m,theta(j), azi(k))
c$$$     &           *conjg(sph_harm(l,m,theta(j), azi(k)))
c$$$     &           *sin(theta(j))*dphi*dtheta
c$$$         enddo
c$$$      enddo
c$$$      print*, 'l,m the same, res=', res
c$$$
c$$$      lprime = 4
c$$$      mprime = 3
c$$$      
c$$$      res = dcmplx(0d0, 0d0)
c$$$      do k=1, 128
c$$$         do j=1, 64
c$$$            res = res +
c$$$     &           sph_harm(l,m,theta(j), azi(k))
c$$$     &           *conjg(sph_harm(lprime,mprime,theta(j), azi(k)))
c$$$     &           *sin(theta(j))*dphi*dtheta
c$$$         enddo
c$$$      enddo
c$$$      print*, 'l,m different, res=', res
c$$$      stop

c$$$      subroutine get_Ilm_and_Qlm
c$$$      use gboundary
c$$$      use grid
c$$$      use field
c$$$      use param
c$$$      use mpiyes
c$$$      use mpipar
c$$$
c$$$      implicit none
c$$$
c$$$      integer :: i, j, k, l, m, count
c$$$      real*8 :: rad, theta, phi
c$$$      real*8 :: dm
c$$$      complex*16 :: element, element2
c$$$
c$$$      do l=0, lmax         
c$$$         do m= -l, l
c$$$            
c$$$            element = dcmplx(0d0, 0d0)
c$$$            element2= dcmplx(0d0, 0d0)
c$$$
c$$$            do k=ks, ke
c$$$               phi = x3b(k)
c$$$               do j=js, je
c$$$                  theta = x2b(j)
c$$$                  do i=is, ie 
c$$$                     rad = x1b(i)
c$$$
c$$$                     dm = rad**2*sin(theta)*dx1a(i)*dx2a(j)*dx3a(k)
c$$$     &                    *conjg(Y_lm(l, m, j, k))
c$$$                     element = element + dm/(rad**(l+1)) !for Ilm
c$$$                     element2 = element2 + dm*rad**l !for Qlm
c$$$
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$           
c$$$            locI_lm(l,m) = element*4d0*pi/(2d0*l + 1d0)
c$$$            locQ_lm(l,m) = element2*4d0*pi/(2d0*l + 1d0)
c$$$
c$$$         enddo
c$$$      enddo
c$$$
c$$$!      print*, 'got the local Qlm'
c$$$
c$$$
c$$$c     in theory, got local processor contribution to Qlm. now sum over all processors
c$$$
c$$$!      do l=0, lmax
c$$$!         count = 2*l + 1 !number of elements to sum
c$$$         
c$$$!         call mpi_allreduce(locI_lm(l, -l), globI_lm(l, -l), count, !!fortran will start increasing the first index when going through memory!!! but we want to fix the
c$$$!     &        MPI_DOUBLE_COMPLEX, MPI_SUM, comm3d,ierr)             !! first index(the l index) 
c$$$!
c$$$!         call mpi_allreduce(locQ_lm(l, -l), globQ_lm(l, -l), count, 
c$$$!     &        MPI_DOUBLE_COMPLEX, MPI_SUM, comm3d,ierr)
c$$$!      enddo
c$$$
c$$$
c$$$! i think the stuff above is VERY WRONG because of how fortran stores arrays. below just sum everything
c$$$
c$$$      call mpi_allreduce(locI_lm, globI_lm, (lmax+1)*(2*lmax+1), 
c$$$     &        MPI_DOUBLE_COMPLEX, MPI_SUM, comm3d,ierr)
c$$$
c$$$      call mpi_allreduce(locQ_lm, globQ_lm,(lmax+1)*(2*lmax+1), 
c$$$     &        MPI_DOUBLE_COMPLEX, MPI_SUM, comm3d,ierr)
c$$$
c$$$
c$$$!      print*, 'got the global Qlm'
c$$$
c$$$      return
c$$$      end
c$$$
c$$$      subroutine get_gpiib
c$$$      use grid
c$$$      use gboundary
c$$$      use bndry
c$$$      use field
c$$$
c$$$      implicit none
c$$$
c$$$      integer :: j, k, l, m 
c$$$      real*8 :: innerrad
c$$$      complex*16 :: contribution
c$$$
c$$$      innerrad = x1b(is-1)
c$$$
c$$$
c$$$      do k=ks, ke
c$$$         do j=js, je
c$$$           
c$$$            contribution = dcmplx(0d0, 0d0)
c$$$            
c$$$            do l=0, lmax
c$$$               do m=-l, l
c$$$                  
c$$$                  contribution = contribution
c$$$     &                 + innerrad**l*Y_lm(l,m,j,k)*globI_lm(l,m)
c$$$               enddo
c$$$            enddo
c$$$           
c$$$            gpiib(j,k,1) = dble(contribution)
c$$$            if(gpiib(j,k,1).le.0d0) print*, 'negative phi iib'
c$$$!            print*, (gpiib(j,k,1) - gp(is,j,k))/gp(is,j,k)
c$$$         enddo
c$$$      enddo
c$$$
c$$$!      print*, 'got inner boundary potential'
c$$$
c$$$      return
c$$$      end
c$$$
c$$$      subroutine get_gpoib
c$$$      use grid
c$$$      use gboundary
c$$$      use bndry
c$$$      use field
c$$$
c$$$      implicit none
c$$$
c$$$      integer :: j, k, l, m 
c$$$      real*8 :: outerrad
c$$$      complex*16 :: contribution
c$$$
c$$$      outerrad = x1b(ie+1)
c$$$
c$$$
c$$$      do k=ks, ke
c$$$         do j=js, je
c$$$           
c$$$            contribution = dcmplx(0d0, 0d0)
c$$$            
c$$$            do l=0, lmax
c$$$               do m=-l, l
c$$$                  
c$$$                  contribution = contribution
c$$$     &                 + Y_lm(l,m,j,k)*globQ_lm(l,m)/(outerrad**(l+1))
c$$$               enddo
c$$$            enddo
c$$$           
c$$$            gpoib(j,k,1) = dble(contribution)
c$$$            if(gpoib(j,k,1).le.0d0) print*, 'negative phi oib'
c$$$!            print*, gpoib(j,k,1) - gp(ie,j,k)
c$$$         enddo
c$$$      enddo
c$$$
c$$$!      print*, 'got inner boundary potential'
c$$$
c$$$      return
c$$$      end


c$$$      REAL*8 FUNCTION plgndr(l,m,x)
c$$$
c$$$      implicit none
c$$$
c$$$      INTEGER :: l,m
c$$$      REAL*8 :: x
c$$$      INTEGER :: i,ll
c$$$      REAL*8 :: fact,pll,pmm,pmmp1,somx2
c$$$      if(m.lt.0.or.m.gt.l.or.dabs(x).gt.1d0) then
c$$$         print*,'bad arguments in plgndr'
c$$$         stop
c$$$      endif
c$$$      pmm=1d0
c$$$      if(m.gt.0) then
c$$$        somx2=sqrt((1d0-x)*(1d0+x))
c$$$        fact=1d0
c$$$        do 11 i=1,m
c$$$          pmm=-pmm*fact*somx2
c$$$          fact=fact+2d0
c$$$11      continue
c$$$      endif
c$$$      if(l.eq.m) then
c$$$        plgndr=pmm
c$$$      else
c$$$        pmmp1=x*(2*m+1)*pmm
c$$$        if(l.eq.m+1) then
c$$$          plgndr=pmmp1
c$$$        else
c$$$          do 12 ll=m+2,l
c$$$            pll=(x*(2*ll-1d0)*pmmp1-(ll+m-1d0)*pmm)/(ll-m)
c$$$            pmm=pmmp1
c$$$            pmmp1=pll
c$$$12        continue
c$$$          plgndr=pll
c$$$        endif
c$$$      endif
c$$$      return
c$$$      END

      real*8 function factorial(l,m)
c     calculate (l-m)! / (l+m)!
      
      implicit none

      integer :: l, m, diff, lpm, j
      real*8 :: result
      if(m.lt.0.or.m.gt.l) then
         print*,'bad arguments in factorial'
         stop
      endif
      
      if(m.eq.0) then
         factorial = 1d0
         return
      endif

      lpm = l + m
      if(lpm.eq.0) then
         factorial = 1d0
         return
      endif

      diff = 2*m !number of terms in multiplication
      
      result = dble(1d0)

      do j = 1, diff
         result = result*(lpm - (j-1))
      enddo

      factorial = 1d0/result
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c  To get dgamma,  "send dgamma from fnlib".
c  To get d1mach, mail netlib
c       send d1mach from core
c
      subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
c
c           this set of routines computes the nodes t(j) and weights
c        w(j) for gaussian-type quadrature rules with pre-assigned
c        nodes.  these are used when one wishes to approximate
c
c                 integral (from a to b)  f(x) w(x) dx
c
c                              n
c        by                   sum w  f(t )
c                             j=1  j    j
c
c        (note w(x) and w(j) have no connection with each other.)
c        here w(x) is one of six possible non-negative weight
c        functions (listed below), and f(x) is the
c        function to be integrated.  gaussian quadrature is particularly
c        useful on infinite intervals (with appropriate weight
c        functions), since then other techniques often fail.
c
c           associated with each weight function w(x) is a set of
c        orthogonal polynomials.  the nodes t(j) are just the zeroes
c        of the proper n-th degree polynomial.
c
c     input parameters (all real numbers are in double precision)
c
c        kind     an integer between 1 and 6 giving the type of
c                 quadrature rule:
c
c        kind = 1:  legendre quadrature, w(x) = 1 on (-1, 1)
c        kind = 2:  chebyshev quadrature of the first kind
c                   w(x) = 1/sqrt(1 - x*x) on (-1, +1)
c        kind = 3:  chebyshev quadrature of the second kind
c                   w(x) = sqrt(1 - x*x) on (-1, 1)
c        kind = 4:  hermite quadrature, w(x) = exp(-x*x) on
c                   (-infinity, +infinity)
c        kind = 5:  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
c                   beta on (-1, 1), alpha, beta .gt. -1.
c                   note: kind=2 and 3 are a special case of this.
c        kind = 6:  generalized laguerre quadrature, w(x) = exp(-x)*
c                   x**alpha on (0, +infinity), alpha .gt. -1
c
c        n        the number of points used for the quadrature rule
c        alpha    real parameter used only for gauss-jacobi and gauss-
c                 laguerre quadrature (otherwise use 0.d0).
c        beta     real parameter used only for gauss-jacobi quadrature--
c                 (otherwise use 0.d0)
c        kpts     (integer) normally 0, unless the left or right end-
c                 point (or both) of the interval is required to be a
c                 node (this is called gauss-radau or gauss-lobatto
c                 quadrature).  then kpts is the number of fixed
c                 endpoints (1 or 2).
c        endpts   real array of length 2.  contains the values of
c                 any fixed endpoints, if kpts = 1 or 2.
c        b        real scratch array of length n
c
c     output parameters (both double precision arrays of length n)
c
c        t        will contain the desired nodes.
c        w        will contain the desired weights w(j).
c
c     underflow may sometimes occur, but is harmless.
c
c     references
c        1.  golub, g. h., and welsch, j. h., "calculation of gaussian
c            quadrature rules," mathematics of computation 23 (april,
c            1969), pp. 221-230.
c        2.  golub, g. h., "some modified matrix eigenvalue problems,"
c            siam review 15 (april, 1973), pp. 318-334 (section 7).
c        3.  stroud and secrest, gaussian quadrature formulas, prentice-
c            hall, englewood cliffs, n.j., 1966.
c
c        original version 20 jan 1975 from stanford
c        modified 21 dec 1983 by eric grosse
c          imtql2 => gausq2
c          hex constant => d1mach (from core library)
c          compute pi using datan
c          removed accuracy claims, description of method
c          added single precision version
c
      double precision b(n), t(n), w(n), endpts(2), muzero, t1,
     &     gam, solve, dsqrt, alpha, beta
c     
      call class (kind, n, alpha, beta, b, t, muzero)
c     
c     the matrix of coefficients is assumed to be symmetric.
c     the array t contains the diagonal elements, the array
c     b the off-diagonal elements.
c     make appropriate changes in the lower right 2 by 2
c     submatrix.
c     
      if (kpts.eq.0)  go to 100
      if (kpts.eq.2)  go to  50
c     
c           if kpts=1, only t(n) must be changed
c
      t(n) = solve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      go to 100
c
c     if kpts=2, t(n) and b(n-1) must be recomputed
c     
 50   gam = solve(endpts(1), n, t, b)
      t1 = ((endpts(1) - endpts(2))/(solve(endpts(2), n, t, b) - gam))
      b(n-1) = dsqrt(t1)
      t(n) = endpts(1) + gam*t1
c     
c     note that the indices of the elements of b run from 1 to n-1
c     and thus the value of b(n) is arbitrary.
c     now compute the eigenvalues of the symmetric tridiagonal
c     matrix, which has been modified as necessary.
c     the method used is a ql-type method with origin shifting
c     
 100  w(1) = 1.0d0
      do 105 i = 2, n
 105     w(i) = 0.0d0
c     
         call gausq2 (n, t, b, w, ierr)
         do 110 i = 1, n
 110        w(i) = muzero * w(i) * w(i)
c     
            return
            end
c     
c     
c     
      double precision function solve(shift, n, a, b)
c     
c     this procedure performs elimination to solve for the
c     n-th component of the solution delta to the equation
c     
c     (jn - shift*identity) * delta  = en,
c     
c     where en is the vector of all zeroes except for 1 in
c     the n-th position.
c     
c     the matrix jn is symmetric tridiagonal, with diagonal
c     elements a(i), off-diagonal elements b(i).  this equation
c     must be solved to obtain the appropriate changes in the lower
c     2 by 2 submatrix of coefficients for orthogonal polynomials.
c     
c     
      double precision shift, a(n), b(n), alpha
c     
      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
 10      alpha = a(i) - shift - b(i-1)**2/alpha
         solve = 1.0d0/alpha
         return
         end
c     
c     
c     
      subroutine class(kind, n, alpha, beta, b, a, muzero)
c     
c     this procedure supplies the coefficients a(j), b(j) of the
c     recurrence relation
c     
c     b p (x) = (x - a ) p   (x) - b   p   (x)
c     j j            j   j-1       j-1 j-2
c     
c     for the various classical (normalized) orthogonal polynomials,
c     and the zero-th moment
c     
c     muzero = integral w(x) dx
c     
c     of the given polynomial's weight function w(x).  since the
c     polynomials are orthonormalized, the tridiagonal matrix is
c     guaranteed to be symmetric.
c     
c     the input parameter alpha is used only for laguerre and
c     jacobi polynomials, and the parameter beta is used only for
c     jacobi polynomials.  the laguerre and jacobi polynomials
c     require the gamma function.
c     
      double precision a(n), b(n), muzero, alpha, beta
      double precision abi, a2b2, dgamma, pi, dsqrt, ab
c     
      pi = 4.0d0 * datan(1.0d0)
      nm1 = n - 1
      go to (10, 20, 30, 40, 50, 60), kind
c     
c     kind = 1:  legendre polynomials p(x)
c     on (-1, +1), w(x) = 1.
c     
 10   muzero = 2.0d0
      do 11 i = 1, nm1
         a(i) = 0.0d0
         abi = i
 11      b(i) = abi/dsqrt(4*abi*abi - 1.0d0)
         a(n) = 0.0d0
         return
c     
c     kind = 2:  chebyshev polynomials of the first kind t(x)
c     on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
c     
 20      muzero = pi
         do 21 i = 1, nm1
            a(i) = 0.0d0
 21         b(i) = 0.5d0
            b(1) = dsqrt(0.5d0)
            a(n) = 0.0d0
            return
c     
c     kind = 3:  chebyshev polynomials of the second kind u(x)
c     on (-1, +1), w(x) = sqrt(1 - x*x)
c     
 30         muzero = pi/2.0d0
            do 31 i = 1, nm1
               a(i) = 0.0d0
 31            b(i) = 0.5d0
               a(n) = 0.0d0
               return
c     
c     kind = 4:  hermite polynomials h(x) on (-infinity,
c     +infinity), w(x) = exp(-x**2)
c     
 40            muzero = dsqrt(pi)
               do 41 i = 1, nm1
                  a(i) = 0.0d0
 41               b(i) = dsqrt(i/2.0d0)
                  a(n) = 0.0d0
                  return
c     
c     kind = 5:  jacobi polynomials p(alpha, beta)(x) on
c     (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
c     beta greater than -1
c     
 50               ab = alpha + beta
                  abi = 2.0d0 + ab
                  muzero = 2.0d0 ** (ab + 1.0d0) * dgamma(alpha + 1.0d0) * dgamma(
     x                 beta + 1.0d0) / dgamma(abi)
                  a(1) = (beta - alpha)/abi
                  b(1) = dsqrt(4.0d0*(1.0d0 + alpha)*(1.0d0 + beta)/((abi + 1.0d0)*
     1                 abi*abi))
                  a2b2 = beta*beta - alpha*alpha
                  do 51 i = 2, nm1
                     abi = 2.0d0*i + ab
                     a(i) = a2b2/((abi - 2.0d0)*abi)
 51                  b(i) = dsqrt (4.0d0*i*(i + alpha)*(i + beta)*(i + ab)/
     1                    ((abi*abi - 1)*abi*abi))
                     abi = 2.0d0*n + ab
                     a(n) = a2b2/((abi - 2.0d0)*abi)
                     return
c     
c     kind = 6:  laguerre polynomials l(alpha)(x) on
c     (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
c     than -1.
c     
 60                  muzero = dgamma(alpha + 1.0d0)
                     do 61 i = 1, nm1
                        a(i) = 2.0d0*i - 1.0d0 + alpha
 61                     b(i) = dsqrt(i*(i + alpha))
                        a(n) = 2.0d0*n - 1 + alpha
                        return
                        end
c     
c     
      subroutine gausq2(n, d, e, z, ierr)
c     
c     this subroutine is a translation of an algol procedure,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c     this is a modified version of the 'eispack' routine imtql2.
c     
c     this subroutine finds the eigenvalues and first components of the
c     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
c     method.
c     
c     on input:
c     
c        n is the order of the matrix;
c
c     d contains the diagonal elements of the input matrix;
c
c        e contains the subdiagonal elements of the input matrix
c          in its first n-1 positions.  e(n) is arbitrary;
c
c        z contains the first row of the identity matrix.
c
c      on output:
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1, 2, ..., ierr-1;
c
c        e has been destroyed;
c
c        z contains the first components of the orthonormal eigenvectors
c          of the symmetric tridiagonal matrix.  if an error exit is
c          made, z contains the eigenvectors associated with the stored
c          eigenvalues;
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     ------------------------------------------------------------------
c
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8 d(n), e(n), z(n), b, c, f, g, p, r, s, machep
      real*8 dsqrt, dabs, dsign, d1mach
c     
      machep=d1mach(4)
c     
      ierr = 0
      if (n .eq. 1) go to 1001
c     
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
c     :::::::::: look for small sub-diagonal element ::::::::::
 105     do 110 m = l, n
            if (m .eq. n) go to 120
            if (dabs(e(m)) .le. machep * (dabs(d(m)) + dabs(d(m+1))))
     x           go to 120
 110     continue
c     
 120     p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     :::::::::: form shift ::::::::::
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = dsqrt(g*g+1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     
c     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (dabs(f) .lt. dabs(g)) go to 150
            c = g / f
            r = dsqrt(c*c+1.0d0)
            e(i+1) = f * r
            s = 1.0d0 / r
            c = c * s
            go to 160
 150        s = f / g
            r = dsqrt(s*s+1.0d0)
            e(i+1) = g * r
            c = 1.0d0 / r
            s = s * c
 160        g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     :::::::::: form first component of vector ::::::::::
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
 200        z(i) = c * z(i) - s * f
c     
            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0d0
            go to 105
 240     continue
c     
c     :::::::::: order eigenvalues and eigenvectors ::::::::::
         do 300 ii = 2, n
            i = ii - 1
            k = i
            p = d(i)
c     
            do 260 j = ii, n
               if (d(j) .ge. p) go to 260
               k = j
               p = d(j)
 260        continue
c     
            if (k .eq. i) go to 300
            d(k) = d(i)
            d(i) = p
            p = z(i)
            z(i) = z(k)
            z(k) = p
 300     continue
c     
         go to 1001
c     :::::::::: set error -- no convergence to an
c     eigenvalue after 30 iterations ::::::::::
 1000    ierr = l
 1001    return
c     :::::::::: last card of gausq2 ::::::::::
         end




!        contribution = dcmplx(0d0, 0d0)
!         do k=ks, ke
!            do j=js, je
!               contribution = contribution 
!     &              +d(i,j,k)*conjg(Y_ml( m,l, j, k))
!     &              *dx3a(k)*dvl2a(j)               
!            enddo
!         enddo





c$$$      contribution = sum(d(i,js:je, ks:ke)
c$$$     &        *Y_ml_star(m,l,js:je,ks:ke))


!              contribution = dcmplx(0d0, 0d0)               
!               do l=lmin, lmax
!                  mrange = min(mmax, l)
!                   do m=-mrange, mrange
!                     
!                     contribution = contribution 
!     &                     +pot_ml_3d_loc(m,l,i)*Y_ml(m,l ,j,k)
!                  enddo
!               enddo


!               contribution = dcmplx(0d0, 0d0)               
!               do l=lmin, lmax
!                  mrange = min(mmax, l)
!                  do m=-mrange, mrange
!                     
!                     contribution = contribution 
!     &                     +rho_ml_3d(m,l,i)*Y_ml(m, l ,j,k)
!                  enddo
!               enddo



c$$$               rgpin  = dcmplx(0d0, 0d0)
c$$$               rgpout = dcmplx(0d0, 0d0)
c$$$               do l=lmin, lmax
c$$$                  mrange = min(mmax, l)
c$$$                   do m=-mrange, mrange
c$$$                     rgpin = rgpin
c$$$     &                     +pot_ml_3d_loc(m,l,is-1)*Y_ml(m,l ,j,k)
c$$$                     rgpout = rgpout
c$$$     &                     +pot_ml_3d_loc(m,l,ie+1)*Y_ml(m,l ,j,k) 
c$$$                  enddo
c$$$               enddo


c$$$               tgpin  = dcmplx(0d0, 0d0)
c$$$               tgpout = dcmplx(0d0, 0d0)
c$$$               do l=lmin, lmax
c$$$                  mrange = min(mmax, l)
c$$$                   do m=-mrange, mrange
c$$$                     tgpin = tgpin
c$$$     &                     +pot_ml_3d_loc(m,l,i)*Y_ml(m,l ,js-1, k)
c$$$                     tgpout = tgpout
c$$$     &                     +pot_ml_3d_loc(m,l,i)*Y_ml(m,l ,js+1, k)
c$$$                  enddo
c$$$               enddo




      subroutine filon_sin ( ftab, a, b, ntab, t, result )
      use grid
!
!***********************************************************************
!
!! FILON_SIN uses Filon's method on integrals with a sine factor.
!
!
!  Discussion:
!
!    The integral to be approximated has the form
!
!      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
!
!    where T is user specified.
!
      implicit none
!     
      integer ntab
!     
      real*8 :: a, alpha, b, beta,
     &     cost, gamma, h, result,
     &     s2n, s2nm1, sint, t, theta
      
      real*8 :: ftab(ntab), xtab(ntab)
!     
      if ( a == b ) then
         result = 0.0d0
         return
      end if
      
      if ( ntab <= 1 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
         write ( *, '(a)' ) '  NTAB < 2'
         write ( *, '(a,i6)' ) '  NTAB = ',ntab
         stop
      end if
      
      if ( mod ( ntab, 2 ) /= 1 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
         write ( *, '(a)' ) '  NTAB must be odd.'
         write ( *, '(a,i6)' ) '  NTAB = ',ntab
         stop
      end if
!     
!     Set up a vector of the NTAB X values.
!     
      call rvec_even ( a, b, ntab, xtab )
      h = ( b - a ) / dble ( ntab - 1d0 )
!       xtab(1:ntab) = x3b(ks:ke)
!       h = dx3b(ks)


      theta = t * h
      sint = sin ( theta )
      cost = cos ( theta )
      
      alpha = ( theta**2 + theta * sint * cost 
     &     - 2d0 * sint**2 ) / theta**3
      
      beta = ( 2d0 * theta + 2d0 * theta * cost**2 
     &     - 4d0 * sint * cost ) / theta**3
      
      gamma = 4.0d0* ( sint - theta * cost ) / theta**3

     
      s2n = 
     &     sum ( ftab(1:ntab:2) * sin ( t * xtab(1:ntab:2) ) ) 
     &     - 0.5d0 * ( ftab(ntab) * sin ( t * xtab(ntab) )    
     &     + ftab(1) * sin ( t * xtab(1) ) )
      
      s2nm1 = sum ( ftab(2:ntab-1:2) * sin ( t * xtab(2:ntab-1:2) ) )
      
      result = h * ( 
     &     alpha * ( ftab(1) * cos ( t * xtab(1) ) 
     &     - ftab(ntab) * cos ( t * xtab(ntab) ) ) 
     &     + beta * s2n 
     &     + gamma * s2nm1 )
      
      return
      end



!     &     +(h/2d0)* ( ftab_cos(ntab) * cos ( t * xtab(ntab) )     
!     &     + ftab_cos(1) * cos ( t * xtab(1) ) )

!     &      +(h/2d0) * ( ftab_sin(ntab) * sin ( t * xtab(ntab) )    
!     &      + ftab_sin(1) * sin ( t * xtab(1) ) ) 


!            call cubint(phiint_real, phiaxis, ntab_k, 1, !shit
!     &           ntab_k, res_real, error)
!            call cubint(phiint_imag, phiaxis, ntab_k, 1,
!     &           ntab_k, res_imag, error)

!            call cspint (phiint_real, phiaxis, ntab_k, x3a(ks), !shit
!     &              x2a(ke+1), ywork_k, ework_k, work_k, res_real)
!            call cspint (phiint_imag, phiaxis, ntab_k, x3a(ks),
!     &              x2a(ke+1), ywork_k, ework_k, work_k, res_imag)

!            call simpsn((x3a(ke+1)-x3a(ks))/ntab_k, !ah shit
!     &           phiint_real, ntab_k, res_real)
!            call simpsn((x3a(ke+1)-x3a(ks))/ntab_k,
!     &           phiint_imag, ntab_k, res_imag)

!            call cspint (phiint_real, phiaxis, ntab_k, x3a(ks), !few percent error at high m, but non-axisymmetric even when it should be
!     &           x3a(ke+1), ywork_k, ework_k, work_k, res_real) 
!            call cspint (phiint_imag, phiaxis, ntab_k, x3a(ks),
!     &           x3a(ke+1), ywork_k, ework_k, work_k, res_imag) 
            
!            call plint (phiint_real, phiaxis, ntab_k, x3a(ks) , !shit
!     &           x3a(ke+1), res_real)
!            call plint (phiint_imag, phiaxis, ntab_k, x3a(ks) ,
!     &           x3a(ke+1), res_imag)

!            call hiordq (ntab_k, phiint_real, dx3a(ks)!shit
!     &           , hwork2, res_real )
!            call hiordq (ntab_k, phiint_imag, dx3a(ks)
!     &           , hwork2, res_imag )



!            call eightP_newton_cotes(phiint, ntab_k, !massive increase in error
!     &           (x3a(ke+1)-x3a(ks))/ntab_k, result2)

!            call fourP_newton_cotes(phiint, ntab_k,!massive increase in error
!     &           (x3a(ke+1)-x3a(ks))/ntab_k, result2) 

!            call cubint(phiint, phiaxis, ntab_k, 1, !increases error
!     &           ntab_k, result2, error)

             
!            call simpsn((x3a(ke+1)-x3a(ks))/ntab_k,!increases error
!     &              phiint, ntab_k, result2)

!               call eightP_newton_cotes(thetaint*sintheta, ntab_j, !error 0.11 ! can paralelise
!     &              (x2a(je+1)-x2a(js))/ntab_j, result)            
               
!               call eightP_newton_cotes_v2(thetaint*sintheta, ntab_j,!error 0.11!
!     &              (x2a(je+1)-x2a(js))/ntab_j, result) 
                    
!               call fourP_newton_cotes(thetaint*sintheta, ntab_j,  0.24 error
!     &              (x2a(je+1)-x2a(js))/ntab_j, result) 
               
!               call fourP_newton_cotes_open(thetaint*sintheta, ntab_j, !error 0.25!
!     &              (x2a(je+1)-x2a(js))/dble(ntab_j), result)
               



!
!     iterate the density field so it's consistent with self-gravity potential in the vertical direction
!     but try to keep total disc mass the same. newton raphson iteration. 
! 
       if(vertSG.eq..true.) then
          i=0
          beta_out    = 1.3
          dbeta_out   = 0.05
          diff        = 1d0
          do while(abs(diff).gt.tol)
             dsg = d            !work on dsg array
             call verticalSG(beta_out, dsg, sgdiscmass)
             diff = (sgdiscmass - discmass)/discmass
             
             dsg = d 
             call verticalSG(beta_out + dbeta_out, dsg, sgdiscmass)
             diffnew = (sgdiscmass - discmass)/discmass

             beta_out = beta_out - diff*dbeta_out/(diffnew - diff)

             i = i + 1
             if(i.gt.100) then
                if(myid.eq.0) then
                   print*, 'could not tweak vertical structure , abort'
                   stop
                endif
             endif
          enddo 
          call mpi_barrier(comm3d, ierr)
          d = dsg               !acceptable total mass after tweaking     
       
          
       endif


subroutine verticalSG(beta_out, dsg, sgdiscmass)
!
!     get the modification to vertical density structure due to self-gravity.
!     assume vertical gradient is largest in poisson's equation, so that
!     d^2Phi/dz^2 = 4 pi G rho
!     combine with vertical equilibrium:
!     0 = -(1/rho)dP/dz - dPhi_*/dz - dPhi/dz 
!     
!     assume rho = rho_NSG * beta. such that rho_NSG balances the star potential term
!     need to solve d^2 u/dz^2 = -4 pi G rho_NSG/c^2 exp(u)
!     where u = log beta
!
!     BC: u=du/dz=0 at z=0


      use planet
      use domain
      use grid
      use mpiyes
      use mpipar

      implicit none
      real*8, parameter :: beta_in = 1.0
      integer:: i, j, k, nok, nbad
      real*8, intent(in) :: beta_out
      real*8, intent(out):: sgdiscmass
      real*8, intent(inout) :: dsg(in,jn,kn)
      real*8 :: u(2)
      real*8 :: zmax, dz, dr, bigR, beta, mass,
     .     theta, dtheta, rad, zfinal, deltar
      external :: rkqs, dudz
      common /vertical/ bigR
    
      deltar = rout - rin
!      beta_out = 1.5 !initial guess for outer radius over-density 

      if(myid.eq.0) then
         print*, '--------------------------------------'
         print*, 'tweaking vertical structure to suit SG'
         print*, '--------------------------------------'
      endif

      mass = 0d0
      do k=ks, ke
         do j=js, je
            theta = x2b(j)
            dtheta= dx2a(j)
            do i=is, ie
               rad = x1b(i)
               dz = rad*dtheta
               bigR   = rad*sin(theta)
               zfinal = abs(rad*cos(theta))
!
!     use ode int to integrate from z=0 to |z|
!
               u(1) =log((bigR-rin)*beta_out/deltar
     &                  +(rout-bigR)*beta_in/deltar)
               u(2) = 0d0
               call odeint(u ,2, 0d0, zfinal, 1d-8,dz/20d0,
     &              0d0, nok,nbad, dudz, rkqs)
              
               beta = exp(u(1))
               dsg(i,j,k) = dsg(i,j,k)*beta
               
               mass = mass 
     &        + dsg(i,j,k)*rad**2*sin(theta)*dx1a(i)*dtheta*dx3a(k)
            enddo
         enddo
      enddo
      
      call MPI_allREDUCE(mass, 
     &     sgdiscmass, 1, MPI_DOUBLE_PRECISION,
     &     MPI_SUM, comm3d, ierr)
      
      if(myid.eq.0) then
         print*,'new disc mass, beta_out =', sgdiscmass, beta_out
      endif
    
!      if(discmass.gt.(1d0+tol)*nsgdiscmass) then !total mass now too large
!          beta_out = beta_out*(1d0 - tol)
!          goto 11
!      endif
!      if(discmass.lt.(1d0-tol)*nsgdiscmass) then !total mass now too small
!          beta_out = beta_out*(1d0 + tol)
!          goto 11
!      endif
      return
      end

      subroutine dudz(z, u, du)
      use planet
      use param
      
      implicit none
      
      real*8 :: z, u(2), du(2)
      real*8 :: bigR, dnsg, rad, theta, csq
      real*8, external :: initialdens
      common /vertical/ bigR

      du(1) = u(2)
      
      rad = sqrt(bigR**2 + z**2)
      theta = asin(bigR/rad)
      
      csq = smallh**2/bigR
      dnsg = initialdens(rad, theta, sig, den0, 
     &     innerhole)

      du(2) = -(4d0*pi*dnsg/csq)*exp(u(1))

      return
      end



subroutine density_pert
!
!    random density perturbation
!
      use field
      use grid
      use mpipar
      use mpiyes

      implicit none

      integer :: i,j,k,ii 
      integer*4 :: timeArray(3)
      real :: rand 
      real*8 :: mass, rad, theta, pert, discmass

      mass = 0d0
      call itime(timeArray) 
      ii = rand(timeArray(1)+timeArray(2)+timeArray(3))
      do 290 k=ks,ke
         do 280 j=js,je
            do 270 i=is,ie
               rad   = x1b(i)
               theta = x2b(j)
               pert = 1d0 + (2d0*rand(0) - 1d0)*0.1
               d(i,j,k) = pert*d(i,j,k)
               mass = mass + d(i,j,k)*rad**2*sin(theta)*dx1a(i)
     &              *dx2a(j)*dx3a(k)
 270        continue
 280     continue
 290  continue

      call MPI_REDUCE(mass, discmass, 1, MPI_DOUBLE_PRECISION,
     &     MPI_SUM, 0, comm3d, ierr )
      
      if(myid .eq. 0) print*,'disc mass after random pert =',discmass
      return
      end



      subroutine get_pot_ml(m, l)
!     
!     get the pot_ml's. do the phi-theta integration locally. then store it in appropriate subarray of the global rho_ml ( which 
!     covers the entire radial range   
      
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
     &     res_real2, res_imag2
      complex*16 :: dd(0:nrad+1), du(0:nrad), dl(1:nrad+1), 
     &     trial(0:nrad+1), trialold(0:nrad+1), rhsinit(0:nrad+1)
      complex*16 :: contribution
      
      dblem = dble(m)
      sintheta(1:ntab_j) = sin(x2b(js:je))
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     part 1: get the rho_ml's                             c
c                                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     
!     zero the local copy of the global rho_ml
!     
      rho_ml_1d_glob = dcmplx(0d0, 0d0)
      
!     
!     first do local phi-theta integration. active density zones only
!     
      do i=is, ie
         
         do k=ks, ke !fix phi, do theta integration. no m-phi dependence here
                     !note that some routines assume uniform spacing in theta, which is NOT
                     !      the same as uniform spacing in -cos(theta) \eqv vol2 axis!
                     !if the integrator does interpolations, then not good idea to parallelise in that direction
                     ! because we don't know end point interpolants are consistent between CPU's

               thetaint(1:ntab_j) = d(i,js:je,k)*Y_ml_star(m,l,js:je)
c
c     UNEVEN spaced integrators, so int  thetaint * d(-cos theta), thetaaxis holds -cos(theta) data points
c

                result = sum(thetaint(1:ntab_j)*dvl2b(js:je)) ! the standard (constant in each cell) gives 0.011 error
               
!               call cubint(thetaint, thetaaxis, ntab_j, 1,                  !***seems good, less than 0.003 error. thetaaxis is -cos(theta)
!     &              ntab_j, result, error)                                  !but should NOT PARALLELIZE in theta. cubic splined.
               
               
!     call cspint (thetaint, thetaaxis, ntab_j, vol2a(js), !0.011 error. splined. seems OK to parallelise
!     &              vol2a(je+1), ywork_j, ework_j, work_j, result)
!               
!     call davint(thetaaxis, thetaint, ntab_j, vol2a(js), !about 0.011. OK to parallelise
!     &              vol2a(je+1), result, errest)        
!               
!     call simpne ( thetaaxis, thetaint, ntab_j, result) !0.003 error. DO NOT PARALLELIZE
!     
!     call plint (thetaint, thetaaxis, ntab_j, vol2a(js) , !0.011 error. OK to parallel
!     &             vol2a(je+1), result)               
c     
c     EVEN spacing integrators, so int thetaint*sin(theta) dtheta 
c     
!     call hiordq (ntab_j, thetaint*sintheta,(x2a(je+1)-x2a(js)!0.05 (parallel) or 0.003 if single CPU in theta
!     &              )/ntab_j, hwork, result )
!     
!     call simpsn((x2a(je+1)-x2a(js))/ntab_j, !error 0.0025 but do not parallelise (error is 5% if so)
!     &              thetaint*sintheta, ntab_j, result)
!              
               phiint(k - ks + 1)= result
            enddo
            
!     
!     now do phi integration, do real part (cosine) and imaginary part( sine) separately. 
!     
            if(m.eq.0) then     ! no m-phi dependence, use cubic interpolation
               result2 = sum(phiint(1:ntab_k)*dx3b(ks:ke))
               
!               call davint(phiaxis, phiint, ntab_k, x3a(ks), !maintains 0.003 error. no change with CPU number
!     &              x3a(ke+1), result2, errest) 
               
               
!     call cspint (phiint, phiaxis, ntab_k, x3a(ks), !maintains 0.003 error. but increases a little when more cpu's put in phi 
!     &              x3a(ke+1), ywork_k, ework_k, work_k, result2)!fortunetly not obvious problem at tile boundary
!     call hiordq (ntab_k, phiint,(x3a(ke+1)-x3a(ks))/ntab_k!increases error. error doubles when cpu numbers double
!     &           , hwork2, result2 )
!     call plint (phiint, phiaxis, ntab_k, x3a(ks) , !maintains 0.003 error. no change with cpu number
!     &           x3a(ke+1), result2)
               
               contribution = dcmplx(result2)            
            else                !there is phi dependence
               
               phiint_real(1:ntab_k) =  phiint(1:ntab_k)
     &           *cos(dblem*phiaxis(1:ntab_k))
               phiint_imag(1:ntab_k) = -phiint(1:ntab_k)
     &           *sin(dblem*phiaxis(1:ntab_k))
               
               res_real = sum(phiint_real(1:ntab_k)*dx3b(ks:ke)) !actually looks good!
               res_imag = sum(phiint_imag(1:ntab_k)*dx3b(ks:ke))
               
!     call filon(phiint_real(1:ntab_k),
!     &           phiint_imag(1:ntab_k),                            
!     &           ntab_k, dblem, res_real, res_imag)               
               
!               call flinn(phiint_real(1:ntab_k),
!     &              phiint_imag(1:ntab_k),                          
!     &              ntab_k, dblem, res_real, res_imag) 
               
!     call davint(phiaxis, phiint_real, ntab_k, x3a(ks), !seems ok! goes crazy for high m (though error is few per cent). but nonaxisymmetric error
!     &           x3a(ke+1), res_real, errest)  
!     call davint(phiaxis, phiint_imag, ntab_k, x3a(ks), 
!     &           x3a(ke+1), res_imag, errest)  
!     
!     call eightP_newton_cotes(phiint_real, ntab_k, !you need at lease 16 cells per wavelength of the mmax mode. e.g. mmax=16 and Nphi=256 then
!     &              (x3a(ke+1)-x3a(ks))/ntab_k, res_real) !have 256/16=16, that's two newton-cotes intervals per wavelength so you're ok
!     call eightP_newton_cotes(phiint_imag, ntab_k, !if you try mmax = 32, then only one NC interval and the result is shit . same goes for 4P NC method
!     &              (x3a(ke+1)-x3a(ks))/ntab_k, res_imag)
!     
!     call fourP_newton_cotes(phiint_real, ntab_k, !problem if you try 4 cells per mode. e.g. Nphi=128 and mmax=32 then 
!     &           dx3a(ks), res_real)                     !for m=32 it's 4 grids per wave-length=one newton-cotes. not good enough
!     call fourP_newton_cotes(phiint_imag, ntab_k,
!     &           dx3a(ks), res_imag) 
               contribution = dcmplx(res_real,res_imag)  
            endif
            if(sgzsymm.eq..true.) contribution = contribution*2d0 !if assume symmetric about z=0, need to double up to account for lower half disc
            rho_ml_1d_loc(i-is+1) = contribution*4d0*pi*x1b(i)**2 
         enddo
     
!      print*, 'got the 1D rho_lm (local)'
!     
!     put the local 1D array above into the correct subarray of the global version. first set zero it though
!     
      
      ibeg = coords(1)*nrad_loc + 1
!      print*, 'r-coord, ibeg', coords(1), ibeg, ie-is+1, nrad_loc, nrad
      rho_ml_1d_glob(ibeg:ibeg+nrad_loc-1) = rho_ml_1d_loc(1:nrad_loc)

!     
!     now do global sum.
!     
      call mpi_allreduce(rho_ml_1d_glob, rho_ml, nrad, 
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, comm3d,ierr)
!      print*, 'did global sum'


!    save the rho_lm for diagonstic checks if desired
!
      rho_ml_3d(m,l,is:ie) = rho_ml(ibeg:ibeg+nrad_loc-1)
     &     /(4d0*pi*x1b(is:ie)**2)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                          c
c part II: matrix problem                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     
!     matrix problem. construct rhs. interior points first.
!     
      rhs(1:nrad) = rho_ml(1:nrad)

!     
!     boundary condition for zero mass outside boundaries or zero gradient in potential
!
       rhs(0)     = 0d0 
       rhs(nrad+1)= 0d0 

!
!     boundary condition for ghost density being same as last active zones
!     AND d^2Phi/dr^2 = 0 at ghost/active boundary. basically evaluate the 
!     ODE for Phi_ml here and put Phi" = 0. don't think it's working for OBC, though

!       rhs(0)     = rho_ml(1)*exp(-1d0/dlogri)
!      rhs(nrad+1)= rho_ml(nrad)*exp(1d0/dlogri)


!     
!  construct elements of the tridiagonal matrix. note that zgtsv overwrites these 1D arrays so they are
!  no good for the next l
!

!!! interior points, for active zone
      du(1:nrad) = udiag(1:nrad)
      dd(1:nrad) = -2d0*dlogri**2 - l*(l+1d0)
      dl(1:nrad) = ldiag(1:nrad) 

!!! IBC: du(0) and dd(0)
!! for zero mass inside boundary:
         du(0)      = dlogri - l/2d0
         dd(0)      = -(dlogri + l/2d0)
!! for zero radial potential gradient at inner boundary
!      du(0)      = dlogri 
!      dd(0)      = -dlogri
!! for zero second radial potential gradient at inner boundary (use original ODE to relate Phi and Phi')
!       du(0)      = dlogri - (l+1d0)*l/2d0
!       dd(0)      = -(dlogri + l*(l+1d0)/2d0)

!!! OBC: dd(nrad+1) and dl(nrad+1) controls OBC
!!  for zero mass outside boundary
         dd(nrad+1) = dlogri + (l+1d0)/2d0
         dl(nrad+1) = -dlogri + (l+1d0)/2d0 
!! for zero radial potential gradient at outer boundary
!       dd(nrad+1) =  dlogri
!       dl(nrad+1) = -dlogri
!! for zero second radial potential gradient at outer boundary (doesn't work - crashes)
!      dd(nrad+1) = dlogri - l*(l+1d0)/2d0
!      dl(nrad+1) = -(dlogri + l*(l+1d0)/2d0)


!     
!     solve matrix problem
!

!         if((nhy.eq.0).or.(nfail.gt.itmax/2)) then !compute exact solution at very first step, or if 
                                                 !we were using SOR but kept failing to converge witin itmax for the previous m,l's
!            nfail = 0
            call zgtsv(nrad+2, 1, dl, dd, du, rhs, nrad+2, info)
            if(info.ne.0) then
               print*, 'did not solve poisson successfully, info=', info
               stop
            endif
!         else
!!           
!!    we are beyond step 0, so have earlier solution, use that in an SOR algorithm
!!     
!            trial(0:nrad+1) = pot_ml_3d_glob(m,l,0:nrad+1)
!            rhsinit = rhs
!            count = 0
!            resmax= 1d0 
!            do while(resmax.gt.tol)  
!               if(count.gt.itmax) then !too many iterations, just call exact solver
!                   nfail = nfail + 1
!                   call zgtsv(nrad+2, 1, dl, dd, du, rhsinit, nrad+2, 
!     &                  info)
!                   rhs = rhsinit
!                   exit
!                endif 
!               
!               trialold = trial
!               trial(0) = (1d0-relax)*trial(0)
!     &              + relax*(rhsinit(0) - du(0)*trial(1))/dd(0)
!               do i=1, nrad
!                  trial(i) = (1d0 - relax)*trial(i)
!     &                 + relax*(rhsinit(i) - dl(i-1)*trial(i-1) 
!     &                 -   du(i)*trial(i+1))/dd(i)
!               enddo
!               trial(nrad+1) = (1d0-relax)*trial(nrad+1)
!     &              + relax*(rhsinit(nrad+1) - dl(nrad+1)*trial(nrad))
!     &              /dd(nrad+1)
!                resmax = maxval(abs((trial - trialold)/trialold))
!                count = count + 1
!                rhs(0:nrad+1) = trial(0:nrad+1) ! replace rhs with the new interated solution
!!                if(myid.eq.0) print*, count
!             enddo
!         endif

!
!     rhs is now the phi_ml(r) 1D vector that covers the global disk. copy appropriate subsection into local array
!     for potential calculation. save a copy of the global solution. 
      pot_ml_3d_loc(m,l,is-1:ie+1) = rhs(ibeg-1:ibeg+nrad_loc)
!      pot_ml_3d_glob(m,l,0:nrad+1) = rhs(0:nrad+1)

      return
      end



      subroutine get_pot
!
!     get the potential. do ghost zones too, so we don't need to communicate
!
      use field
      use grid
      use gboundary

      implicit none

      integer :: i, j, k, l, m
      integer :: mrange
      complex*16 :: contribution
      
      do k=ks-1, ke+1
         do j=js-1, je+1
            do i=is-1, ie+1               


               contribution = sum(pot_ml_3d_loc(-mmax:mmax, lmin:lmax,i)
     &              *Y_ml(-mmax:mmax, lmin:lmax, j,k))

      
               gp(i,j,k) = -dble(contribution)
!!
!! solution method is for standard convention (negative potential) but zeus standard is positive potential
!! 
            enddo
         enddo
      enddo


!      do k=ks, ke !get the reconstructed density. for diagonistic and testing runs only. comment out otherwise
!         do j=js, je
!            do i=is, ie               
               
               
!                contribution = sum(rho_ml_3d(-mmax:mmax, lmin:lmax, i)
!     &              *Y_ml(-mmax:mmax, lmin:lmax, j, k))


               
!              rho_recons(i,j,k)=(dble(contribution)-d(i,j,k))/d(i,j,k)
               
!            enddo
!         enddo
!      enddo

      return
      end




 subroutine sgboss
!
!     compute the self-potential everywhere using spherical harmonic expansion
!
      use gboundary
      use mpiyes
      use mpipar

      implicit none
      
      integer :: j, k, l,m, mprime, lprime
      integer :: mrange
      real*8 :: azi(128), theta(64)
      real*8 :: dphi, dtheta
      complex*16 :: sph_harm, res
      
      do l=lmin, lmax
         mrange = min(mmax, l)
         do m=-mrange, mrange
            if(mod(l+m,2).eq.0) call get_pot_ml(m, l) ! only do l+m = even modes

 !           call get_pot_ml(m, l)
         enddo
      enddo

      call mpi_barrier(comm3d, ierr)
      call get_pot
      
      return
      end
