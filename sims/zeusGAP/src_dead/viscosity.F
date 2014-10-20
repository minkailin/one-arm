      subroutine get_stress_tensor(oldv1, oldv2, oldv3)
      use param
      use grid
      use field
      use mpiyes
      use mpipar
      implicit none
      integer :: ibeg,iend,jbeg,jend,kbeg,kend
      real*8 :: oldv1(in,jn,kn), oldv2(in,jn,kn), oldv3(in,jn,kn)
!!     boundary communications (some are already done in srcstep) (label these "skip")
       nreq = 0
       nsub = nsub + 1
!!     additional communication to get T_rp
!      call bvald  (0,1,0,0,0,0,d)
!      call bvalv3 (0,1,0,0,0,0,oldv3)
!      call bvalv1 (0,1,0,0,1,1,oldv1)
      
!!     additional communication to get T_tp
!      call bvalv2 (0,0,0,0,1,1,oldv2)
!      call bvalv3 (0,0,0,1,0,0,oldv3)
!!     additional communication to get T_rt
!      call bvalv1 (0,0,1,1,0,0,oldv1)
!      call bvalv2 (0,1,0,0,0,0,oldv2)
!!     additional commmuncation to get diag cpt of stress tensor. these are cell centred, as well as
!!     velocity divergences required to compute it. need velocities at ALL first boundaries
!      call bvalv1 (1,0,0,0,0,0,oldv1)
!      call bvalv2 (0,0,1,0,0,0,oldv1)
!!    combination of the above ("skipcomb"). we need to delay resetting bvstat to after viscous update. 
!      call bvald  (0,1,0,0,0,0,d)
!      call bvalv1 (1,1,1,1,1,1,oldv1)
!      call bvalv2 (0,1,1,0,1,1,oldv2)
!      call bvalv3 (0,1,0,1,0,0,oldv3)
!!     request all boundary information("reqALL"). we need not skip the bvstat reset  just after forces_d (do it again after viscous update.)
       call bvald   (1,1,1,1,1,1, d)
       call bvalv1  (1,1,1,1,1,1,oldv1) 
       call bvalv2  (1,1,1,1,1,1,oldv2)
       call bvalv3  (1,1,1,1,1,1,oldv3)
       if(nreq.ne.0) call MPI_WAITALL(nreq, req, stat, ierr)
       call get_T_components(oldv1, oldv2, oldv3)
!      print*, 'got stress tensor'
      return
      end
      subroutine get_T_components(oldv1, oldv2, oldv3)
      use param
      use grid
      use field
c
c     get the edge-centred T_rphi stress tensor. we need to finite difference it in the phi direction for updating
c     vrad, and finite difference in the radial direction for updating vphi
c
c     require velocity field AFTER all boundaries have been communicated
      implicit none
      real*8, parameter :: two_thirds = 2d0/3d0
      integer:: ibeg, iend, jbeg, jend, kbeg, kend, i, j, k,
     & ip1, jp1, kp1
      real*8 :: oldv1(in,jn,kn), oldv2(in,jn,kn), oldv3(in,jn,kn)
      real*8 :: bigR, den, mu, domega_dr, dvr_dphi, dvtheta_dphi, 
     &     dvphi_dtheta, dvr_dtheta, dvtheta_dr, div_v,
     &     dvr_dr, vr_cell, vtheta_cell
      real*8, external :: kinematic_visc
!      print*, 'getting Trphi'
      do k=ks-1, ke + 1
         kp1 = k + 1
         do j=js-1, je+1
            jp1 = j + 1
            do i=is-1, ie + 1
               ip1 = i + 1
!     cell centered stuff
               bigR = x1b(i)*sin(x2b(j))
               mu = kinematic_visc(bigR)*d(i,j,k)
               vr_cell = 0.5*(oldv1(i,j,k) + oldv1(i+1,j,k))
               vtheta_cell = 0.5*(oldv2(i,j,k) + oldv2(i,j+1,k))
!     velocity divergence of this cell
               div_v = ( g2a(ip1) * g31a(ip1) * oldv1(ip1,j,k)
     1              - g2a(i  ) * g31a(i  ) * oldv1(i  ,j,k) )
     2              *                         dvl1ai(i)
     3              + ( g32a(jp1) * oldv2(i,jp1,k)
     4              - g32a(j  ) * oldv2(i,j  ,k) )
     5              *   g2bi(i)             * dvl2ai(j)
     6              + ( oldv3(i,j,kp1) - oldv3(i,j,k) )
     7              *   g31bi(i) * g32bi(j) * dvl3ai(k)
               div_v = div_v*two_thirds*mu !always have this term for diag cpts 
!     T_rr 
               T_rr(i,j,k) = (oldv1(i+1,j,k) 
     &              - oldv1(i,j,k))*dx1ai(i)*2d0*mu
     &              -div_v
!     T_tt
               T_tt(i,j,k) = 2d0*mu*((oldv2(i,j+1,k)
     &              -oldv2(i,j,k))*dx2ai(j) 
     &              + vr_cell)/x1b(i)
     &              -div_v
!     T_pp     
               T_pp(i,j,k) = (oldv3(i,j,k+1)-oldv3(i,j,k))*dx3ai(k)/bigR
     &              + vr_cell/x1b(i) + vtheta_cell*cos(x2b(j))/bigR
               T_pp(i,j,k) = 2d0*mu*T_pp(i,j,k) -div_v
               
!     T_rphi
               bigR = x1a(i)*sin(x2b(j))
               den = (d(i, j, k) + d(i-1, j, k)  
     &              +d(i-1,j, k-1) + d(i, j, k-1))*0.25d0
               mu = kinematic_visc(bigR)*den
               domega_dr = (oldv3(i ,j ,k)/x1b(i) 
     &              - oldv3(i-1,j,k)/x1b(i-1))*dx1bi(i-1)
               dvr_dphi = (oldv1(i,j,k) - oldv1(i,j,k-1))*dx3bi(k-1)
               T_rp(i,j,k) = mu*(x1a(i)*domega_dr + dvr_dphi/bigR)
!     T_tp
               bigR = x1b(i)*sin(x2a(j))
               den  = (d(i, j, k) + d(i, j-1, k)  
     &              +d(i,j-1, k-1) + d(i, j, k-1))*0.25d0
               mu   = kinematic_visc(bigR)*den
               
               dvtheta_dphi = (oldv2(i,j,k) - oldv2(i,j,k-1))*dx3bi(k-1)
               dvphi_dtheta = (oldv3(i,j,k)/sin(x2b(j)) 
     &              - oldv3(i,j-1,k)/sin(x2b(j-1)))*dx2bi(j-1) !d(vphi/sin(theta))/dtheta
               T_tp(i,j,k) = mu*(dvtheta_dphi/bigR 
     &              + sin(x2a(j))*dvphi_dtheta/x1b(i))
!     T_rt
               bigR = x1a(i)*sin(x2a(j))
               den = (d(i, j, k) + d(i-1, j, k)  
     &              +d(i-1,j-1, k) + d(i, j-1, k))*0.25d0
               mu   = kinematic_visc(bigR)*den
               
               dvr_dtheta = (oldv1(i,j,k) - oldv1(i,j-1,k))*dx2bi(j-1)
               dvtheta_dr = (oldv2(i,j,k)/x1b(i) 
     &              - oldv2(i-1,j,k)/x1b(i-1))*dx1bi(i-1)
               T_rt(i,j,k) = mu*(dvr_dtheta/x1a(i) + x1a(i)*dvtheta_dr)
            enddo
         enddo
      enddo
!      print*, 'got Trphi'
      return
      end
      subroutine viscosity_update(vrad, vtheta, vphi)
      use param
      use field
      use grid
      use root
      use bndry
c
c     update velocity field due to stress tensor. only T_rphi implemented for now. can only update vrad and vtheta (see tassoul)
c
      implicit none
      
      integer :: i, j ,k 
      real*8 :: vrad(in,jn,kn), vtheta(in,jn,kn), vphi(in,jn,kn)
      real*8 :: dvrad, dvtheta, dvphi, den, bigR, contrib, partial 
!    print*, 'performing viscosity update'
      do k=ks, ke      
         do j=js, je
            do i=is, ie
!        do j=js+1, je
!           do i=is+1, ie
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     vrad  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
               den  = 0.5*(d(i,j,k) + d(i-1, j, k)) !just use arithmic averge. not totlly accurate because of non-uniform radial spacing
               bigR = x1a(i)*sin(x2b(j))
!     due to T_rp
               contrib = (T_rp(i,j,k+1) - T_rp(i,j,k))*dx3ai(k)/bigR
               partial =  contrib
!     due to T_rt
               contrib = (T_rt(i,j+1,k)*sin(x2a(j+1))
     &              -T_rt(i,j,k)*sin(x2a(j)))*dvl2ai(j)/x1a(i)
               partial = partial + contrib
!     due to T_rr
!               contrib = (T_rr(i,j,k)*x1b(i)**2 
!     &              - T_rr(i-1,j,k)*x1b(i-1)**2)*dvl1b(i-1)
               contrib = (T_rr(i,j,k)-T_rr(i-1,j,k))*dx1bi(i-1)
     &                    +(T_rr(i,j,k)+T_rr(i-1,j,k))/x1a(i)
               partial = partial + contrib
!     due to T_tt
               contrib = -0.5*(T_tt(i,j,k) + T_tt(i-1,j,k))/x1a(i)
               partial = partial + contrib
!     due to T_pp
               contrib = -0.5*(T_pp(i,j,k) + T_pp(i-1,j,k))/x1a(i)
               partial = partial + contrib
               dvrad = dt*partial/den
               vrad(i,j,k) = vrad(i,j,k) + dvrad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    vtheta  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               den  = 0.5*(d(i,j,k) + d(i,j-1,k))
               bigR = x1b(i)*sin(x2a(j))
!     due to T_tp
               contrib = (T_tp(i,j,k+1) - T_tp(i,j,k))*dx3ai(k)
     &              /bigR
               partial = contrib
!     due to T_rt
               contrib = (T_rt(i+1,j,k)*x1a(i+1)**3
     &              -T_rt(i,j,k)*x1a(i)**3)*dvl1ai(i)/x1b(i)
               partial = partial + contrib
!     due to T_tt
               contrib = (T_tt(i,j,k)*sin(x2b(j))
     &              -T_tt(i,j-1,k)*sin(x2b(j-1)))*dvl2bi(j-1)/x1b(i)
               partial = partial + contrib
!     due to T_pp
               contrib = -0.5*(T_pp(i,j,k) + T_pp(i,j-1,k))
     &              *cos(x2a(j))/bigR
               partial = partial + contrib
               dvtheta = dt*partial/den
               vtheta(i,j,k) = vtheta(i,j,k) + dvtheta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     vphi  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               den   = 0.5*(d(i,j,k) + d(i, j, k-1)) !uniform phi spacing
               bigR  = x1b(i)*sin(x2b(j))
!     due to T_rp               
               contrib = (T_rp(i+1,j,k)*x1a(i+1)**3 
     &              - T_rp(i,j,k)*x1a(i)**3)*dvl1ai(i)/x1b(i)
               partial = contrib
!     due to T_tp
               contrib = (T_tp(i,j+1,k)*sin(x2a(j+1))**2 
     &              - T_tp(i,j,k)*sin(x2a(j))**2)*dvl2ai(j)
     &              /bigR
               partial = partial + contrib
               
!     due to T_pp
               contrib = (T_pp(i,j,k) - T_pp(i,j,k-1))*dx3bi(k-1)/bigR
               partial = partial + contrib
               dvphi = dt*partial/den
               vphi(i,j,k) = vphi(i,j,k) + dvphi
            enddo
         enddo
      enddo
!     
!     make the above loop over i=is+1 to ie, similar for j. these are active for sure. k is periodic so ks to ke are all active
!     if(niis(1).eq.0) update vrad for i=is, js to je
!     if(nijs(1.eq.0) update vtheta for i=is to ie and j=js 
! 
!move below to after forces_d
      do 10 i = 1,6
         bvstat(i,3) = 0        !  v1
         bvstat(i,4) = 0        !  v2
         bvstat(i,5) = 0        !  v3
 10   continue
      
!      print*, 'updated with viscosity'
      return
      end
      real*8 function kinematic_visc(cylind_rad)
      use planet
      use domain     
      implicit none
      real*8, parameter :: width= 2d0
      real*8  :: cylind_rad, cs, bigH
      real*8  :: cylind_min, factor, factor2, dr, dr2
      if((nu.gt.0d0).and.(alpha.gt.0d0)) then
         print*, 'error: cannot define both nu and alpha'
         stop
      endif
      if(nu.gt.0d0) kinematic_visc = nu
      if(alpha.gt.0d0) then
         cs = smallh/sqrt(cylind_rad)
         bigH = smallh*cylind_rad
         kinematic_visc = alpha*cs*bigH
      endif
c     taper the viscosity so it only applies near inner bc and outer bc
c      cylind_min = rin*sin(thetamin)
c      dr         = width*smallh*cylind_min
c      dr2        = width*smallh*rout
c      factor = exp(-0.5*(cylind_rad - cylind_min)**2/dr**2)
c      factor2= exp(-0.5*(cylind_rad - rout)**2/dr2**2) 
c      kinematic_visc = kinematic_visc*(factor+factor2)
c      viscosity jump using tanh 
c      dr = softeps*planetrad
c      factor =0.5d0*( 1d0 - tanh( (cylind_rad - planetrad)/dr ) )      
c      kinematic_visc = kinematic_visc*factor
      return
      end
         
      
      
