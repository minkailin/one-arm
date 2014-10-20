	module vortensity_global
  	use planet
	use domain
  	implicit none
        integer, parameter :: nx=512
  	real*8, parameter :: mstar = 1.0, bigG=1.0, pi=2d0*asin(1d0)
  	real*8 :: sigma_out, omega_out, valpha, beta, sigma_in, qin 
        real*8 :: radius(nx), dr, omega_out_trial
	end module vortensity_global
	
	
	subroutine vortensity(cylind_rad, surfD, dlogsurfD, get_omega_out)
	use vortensity_global
	implicit none
        logical, intent(in) :: get_omega_out 
        integer, parameter :: n=1, lwa=16
	real*8, parameter :: rtol=1d-15 
	real*8, intent(in) :: cylind_rad
	real*8, intent(out) :: surfD, dlogsurfD
	integer :: i, info, iflag
	real*8 :: surfD_floor, den_final, omega_final, rin_cylind
	real*8, external :: sound_speed_sq, 
     & 	 dsound_speed_sq, vortensity_profile, omega_k
	real*8 :: y(2), wa(lwa), x(n), fvec(n)
	external :: funcv1
	common /result/ den_final, omega_final 
        valpha    = sig - 3d0/2d0
        sigma_out = den0
   
        dr = (rout-cylind_rad)/(nx-1d0)
        do i=1, nx
          radius(i) = cylind_rad + dr*(i-1d0)
        enddo
        if(get_omega_out.eqv..true.) then ! first ever call, figure out what omega_out should be
!        qin  = qout*(rout/rin)**(2d0 - sig)
        qin  = innerhole
        !the actual minimum inner cylindrical radius is rin*sin(thetamin)
	rin_cylind = rin*sin(thetamin)
        sigma_in  = sqrt(sound_speed_sq(rin_cylind))*omega_k(rin_cylind)
     &             /(pi*bigG*qin)
      
	beta            = sqrt(1d0 - smallh**2*(2d0+sig)) 
        omega_out_trial = beta*omega_k(rout) !initial guess
        omega_out = omega_out_trial
        x(1) = omega_out
        call hybrd1(funcv1, n, x, fvec, rtol, info, wa, lwa)
        omega_out = x(1) !updated omega_out for later use
       ! print*, omega_out
        else
	x(1) = omega_out 
        call funcv1(n,x,fvec, iflag)
      
	surfD = den_final
        
!! dlogsurfD = dlog (Sigma/H)/dr
	dlogsurfD = cylind_rad*omega_final**2 -dsound_speed_sq(cylind_rad) - 
     & 	bigG*mstar/cylind_rad**2
	dlogsurfD = dlogsurfD/sound_speed_sq(cylind_rad)
        endif
        return
	end subroutine vortensity
  
        subroutine funcv1(n,x,fvec, iflag)
 	use vortensity_global
	implicit none
 	integer, parameter :: neq = 2, itol=1, iopt=0, mf=22, itask=1
 	real*8, parameter :: rtol=1d-15, atol=1d-15
 	integer :: i, istate, lrw, liw, n, iflag
 	integer, allocatable :: iwork(:)
 	real*8 :: t, tout, x(n), fvec(n), y(2), den_final, omega_final 
 	real*8, allocatable :: rwork(:)
 	external :: vdiff, vjac
        common /result/ den_final, omega_final 
 	lrw = 22 +  9*neq + neq**2
 	liw = 20 + neq
 	allocate(rwork(lrw))
 	allocate(iwork(liw))
 	y(1) = sigma_out
        y(2) = x(1)
 	
        do i=nx-1, 1, -1
        t = radius(i+1)
        tout = radius(i)
        istate = 1
        call lsode (vdiff, neq, y, t, tout, itol, rtol,
     &    atol, itask, istate, iopt, rwork, lrw, iwork, liw, vjac, mf)
        enddo
        den_final   = y(1)
        omega_final = y(2) 
        fvec(1) = log( y(1)/sigma_in )
       
        !print*, fvec(1) 
        return
        end subroutine funcv1
	subroutine vjac  
	return
	end subroutine vjac
	subroutine vdiff(neq, r, y, ydot)
	use vortensity_global
	implicit none
	integer :: neq
	real*8 :: y(2), ydot(2), den, omg, r, dden
	real*8, external :: vortensity_profile, sound_speed_sq, 
     &	dsound_speed_sq 
	den = y(1)
	omg = y(2)
	ydot(1) = r*omg**2 - dsound_speed_sq(r) - bigG*mstar/r**2 
     &	+ sound_speed_sq(r)/r
	ydot(1) = ydot(1)*den/sound_speed_sq(r)
	ydot(2) = den*vortensity_profile(r)/r - 2d0*omg/r
 
	return
	end subroutine vdiff
	real*8 function sound_speed_sq(rad)
	use vortensity_global
	implicit none
	real*8 :: rad
	sound_speed_sq = smallh**2*bigG*mstar/rad
	return
	end function sound_speed_sq
	real*8 function dsound_speed_sq(rad)
	use vortensity_global
	implicit none
	real*8 :: rad
	dsound_speed_sq = -smallh**2*bigG*mstar/rad**2
	return
	end function dsound_speed_sq
	real*8 function omega_k(rad)
	use vortensity_global
	implicit none
	real*8 :: rad
  
	omega_k = sqrt(bigG*mstar/rad**3)
	return
	end function omega_k
	real*8 function vortensity_profile(rad)
	use vortensity_global
	implicit none
	real*8 :: rad
	real*8, external :: inner_hole, omega_k
  
	vortensity_profile = omega_out_trial/(2d0*sigma_out)
	vortensity_profile = vortensity_profile*(rad/rout)**valpha 
     &	*inner_hole(rad)
	return
	end function vortensity_profile
