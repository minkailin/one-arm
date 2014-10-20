c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                    E O S                  \\\\\\\\\c
c                            Developed by
c                Laboratory of Computational Astrophysics
c                 University of California at San Diego
c
c=======================================================================
      subroutine eos(ibeg,iend,jbeg,jend,kbeg,kend)
c
      use real_prec
      use config
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use cons
      use eos_par
      use mpiyes
      use mpipar
      use planet
c
      implicit NONE
c
      integer  :: i, j, k, ibeg, iend, jbeg, jend, 
     .            kbeg, kend
      real*8   :: cylind_rad, cs, rad, theta, custcs,
     .            distp, x, y, z, omegap_sq, omegas_sq,
     .            dcs
      real*8   :: csq, bigH
c
      do k = kbeg, kend
       do j = jbeg, jend
        do i = ibeg, iend
         if(leos .eq. 1) then
          if(.not. xiso) then
           if(.not. xtotnrg) then
            p(i,j,k) = gamm1*e(i,j,k)
           else
            p(i,j,k) = gamm1*(
     .                 e(i,j,k) - 0.125D0*d(i,j,k)*(
     .                  (v1(i  ,j  ,k  )+v1(i+1,j  ,k  ))**2
     .                 +(v2(i  ,j  ,k  )+v2(i  ,j+1,k  ))**2
     .                 +(v3(i  ,j  ,k  )+v3(i  ,j  ,k+1))**2 )
     .                        )
           endif
          else
c
c     Locally isothermal equation of state, now with arbitrary power-law
c     With option to do Peplinski (2008)'s EOS, which heats up the gas
c     near the planet. Enabled if parameter hp > 0. 
             rad = x1b(i)
             theta = x2b(j)
             cylind_rad = rad*sin(theta)
c             cs = smallh/sqrt(cylind_rad)
             cs = sqrt(csq(cylind_rad))
             if(hplanet.gt.0d0) then
                x = rad*sin(theta)*cos(x3b(k))
                y = rad*sin(theta)*sin(x3b(k))
                z = rad*cos(theta)
                distp = (x-xp)**2 + (y-yp)**2 +(z-zp)**2
     &               + planetsoft**2
                distp = sqrt(distp)
                omegas_sq = 1d0/cylind_rad**3 ! Assumed G = M_* = 1
                omegap_sq = planetmass/distp**3
                custcs = bigH(cylind_rad)*hplanet*distp
                custcs = custcs*sqrt(omegas_sq + omegap_sq)
                custcs = custcs/(bigH(cylind_rad)**(7d0/2d0)
     &               +(hplanet*distp)**(7d0/2d0))**(2d0/7d0)
                dcs = masstaper*(custcs - cs)
                cs = cs + dcs
             endif
             p(i,j,k) = cs*cs*d(i,j,k)
          endif
         endif ! IDEAL GAS
        enddo ! i
       enddo ! j
      enddo ! k
c
      return
      end
