cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                  c
c     integration routines                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine cubint ( ftab, xtab, ntab, ia, ib, result, error )
!     INTEGRAL (XTAB(IB) to XTAB(IA)) F(X) DX
!     
!     Input, real FTAB(NTAB), contains the tabulated function
!     values, FTAB(I) = F(XTAB(I)).
!     
!     Input, real XTAB(NTAB), contains the points at which the
!     function was tabulated.  XTAB should contain distinct
!     values, given in ascending order.
!     
!     Input, integer NTAB, the number of tabulated points.
!     NTAB must be at least 4.
!     
!     Input, integer IA, the entry of XTAB at which integration
!     is to begin.  IA must be no less than 1 and no greater
!     than NTAB.
!     
!     Input, integer IB, the entry of XTAB at which integration
!     is to end.  IB must be no less than 1 and no greater than
!     NTAB.
!     
!     Output, real RESULT, the approximate value of the
!     integral from XTAB(IA) to XTAB(IB) of the function.
!     
!     Output, real ERROR, an estimate of the error in
!     integration.
!     
      implicit none
!     
      integer ntab, i, ia, ib, ind,
     &     it, j, k
!     
      real*8 :: c, d1, d2, d3, error, 
     &           h1, h2, h3, h4,
     &           r1, r2, r3, r4, result,
     &            s, term
      real*8 :: ftab(ntab), xtab(ntab)
!     
      result = 0.0d0
      error  = 0.0d0
 
      if ( ia == ib ) then
         return
      end if
 
      if ( ntab < 4 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'CUBINT - Fatal error!'
         write ( *, '(a,i6)' ) 
     &        '  NTAB must be at least 4, but input NTAB = ',ntab
         stop
      end if
 
      if ( ia < 1 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'CUBINT - Fatal error!'
         write ( *, '(a,i6)' ) 
     &        '  IA must be at least 1, but input IA = ',ia
         stop
      end if
 
      if ( ia > ntab ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'CUBINT - Fatal error!'
         write ( *, '(a,i6)' ) 
     &        '  IA must be <= NTAB, but input IA=',ia
         stop
      end if
 
      if ( ib < 1 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'CUBINT - Fatal error!'
         write ( *, '(a,i6)' ) 
     &        '  IB must be at least 1, but input IB = ',ib
         stop
      end if
 
      if ( ib > ntab ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'CUBINT - Fatal error!'
       write ( *, '(a,i6)' ) 
     &      '  IB must be <= NTAB, but input IB=',ib
       stop
      end if
!
!  Temporarily switch IA and IB, and store minus sign in IND
!  so that, while integration is carried out from low X's
!  to high ones, the sense of the integral is preserved.
!
      if ( ia > ib ) then
         ind = -1
         it = ib
         ib = ia
         ia = it
      else
         ind = 1
      end if
      
      s = 0.0d+00
      c = 0.0d+00
      r4 = 0.0d+00
      j = ntab-2
      if ( ia < ntab-1 .or. ntab == 4 ) then
         j=max(3,ia)
      end if
      
      k = 4
      if ( ib > 2 .or. ntab == 4 ) then
         k=min(ntab,ib+2)-1
      end if
      
      do i = j, k
         
         if ( i <= j ) then
            
            h2 = xtab(j-1)-xtab(j-2)
            d3 = (ftab(j-1)-ftab(j-2)) / h2
            h3 = xtab(j)-xtab(j-1)
            d1 = (ftab(j)-ftab(j-1)) / h3
            h1 = h2+h3
            d2 = (d1-d3)/h1
            h4 = xtab(j+1)-xtab(j)
            r1 = (ftab(j+1)-ftab(j)) / h4
            r2 = (r1-d1) / (h4+h3)
            h1 = h1+h4
            r3 = (r2-d2) / h1
            
            if ( ia <= 1 ) then
               result = h2 * (ftab(1)+h2*(0.5d0*d3-h2*(d2/6d0-(h2+h3+h3)
     &              *r3/12d0)))
               s = -h2**3 * (h2*(3d0*h2+5d0*h4)+10d0*h3*h1)/60d0
            end if
            
         else
            
            h4 = xtab(i+1)-xtab(i)
            r1 = (ftab(i+1)-ftab(i))/h4
            r4 = h4+h3
            r2 = (r1-d1)/r4
            r4 = r4+h2
            r3 = (r2-d2)/r4
            r4 = (r3-d3)/(r4+h1)
            
         end if
         
         if ( i > ia .and. i <= ib ) then
            
            term = h3*((ftab(i)+ftab(i-1))*0.5d0-h3*h3*(d2+r2+(h2-h4)*r3
     &           )/12d0)
            result = result+term
            c = h3**3*(2d0*h3*h3+5d0*(h3*(h4+h2)+2d0*h2*h4))/120d0
            error = error+(c+s)*r4
            
            if ( i /= j ) then
               s = c
            else
               s = s+c+c
            end if
            
         else
            
            error = error+r4*s
            
         end if
         
         if ( i >= k ) then
            
            if ( ib >= ntab ) then
               term = h4*(ftab(ntab) - h4*(0.5d0*r1+h4*(r2/6d0 +(h3+h3
     &              + h4)*r3/12d0)))
               result = result + term
               error = error - h4**3 * r4 * 
     &              ( h4 * ( 3d0 * h4 + 5d0* h2 ) 
     &              + 10d0* h3 * ( h2 + h3 + h4 ) ) / 60.0d0
            end if
            
            if ( ib >= ntab-1 ) error=error+s*r4
         else
            h1 = h2
            h2 = h3
            h3 = h4
            d1 = r1
            d2 = r2
            d3 = r3
         end if
         
      end do
!     
!     Restore original values of IA and IB, reverse signs
!     of RESULT and ERROR, to account for integration
!     that proceeded from high X to low X.
!     
      if ( ind /= 1 ) then
         it = ib
         ib = ia
         ia = it
         result = -result
         error = -error
      end if
      
      return
      end
      
      subroutine filon ( ftab_cos,ftab_sin, ntab, t, result_cos,
     &     result_sin )
      use grid
      use mpiyes
      use mpipar
!
!***********************************************************************
!
!! FILON_COS uses Filon's method on integrals with a cosine factor.
!
!
!  Discussion:
!
!    The integral to be approximated has the form:
!
!      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
!
!    where T is user specified.
!    Output, real RESULT, the approximate value of the integral.
!
      implicit none
!
      integer :: ntab, ntabp1
!
      real*8, parameter :: tol=1d0/6d0
      real*8, parameter :: a3 = 2d0/45d0, a5=-2d0/315d0,
     &     a7=2d0/4725d0, a9= -8d0/467775d0, 
     &     a11=4d0/8513505d0,
     &     b0 = 2d0/3d0, b2=2d0/15d0, b4=-4d0/105d0,
     &     b6 = 2d0/567d0, b8=-4d0/22275d0,
     &     b10 = 4d0/675675d0, b12=-8d0/58046625d0,
     &     g0=4d0/3d0, g2=-2d0/15d0, g4=1d0/210d0,
     &     g6=-1d0/11340d0, g8 = 1d0/997920d0,
     &     g10=-1d0/129729600d0, g12=1d0/23351328000d0
      real*8 :: a, alpha, b, beta, c2n,
     &          c2nm1, cost, gamma, h, s2n, s2nm1,
     &          result_cos, result_sin, sint, t, theta,
     &          df, f0, fnp1
      real*8 :: ftab_cos(ntab), ftab_sin(ntab), xtab(ntab),
     &     ftab_ce(ntab+1), ftab_se(ntab+1), xtab_e(ntab+1)
!
!      if ( a == b ) then
!         result_cos = 0.0d0
!         result_sin = 0.0d0
!         return
!      end if
!      if ( ntab <= 1 ) then
!         write ( *, '(a)' ) ' '
!         write ( *, '(a)' ) 'FILON_COS - Fatal error!'
!         write ( *, '(a)' ) '  NTAB < 2'
!         write ( *, '(a,i6)' ) '  NTAB = ', ntab
!         stop
!      end if
!      if ( mod ( ntab, 2 ) /= 1 ) then
!         write ( *, '(a)' ) ' '
!         write ( *, '(a)' ) 'FILON_COS - Fatal error!'
!         write ( *, '(a)' ) '  NTAB must be odd.'
!         write ( *, '(a,i6)' ) '  NTAB = ', ntab
!         stop
!      end if
      if ( mod ( ntab, 2 ).ne. 0 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'FILON_COS - Fatal error!'
         write ( *, '(a)' ) 'NTAB must be even.'
         write ( *, '(a,i6)' ) '  NTAB = ', ntab
         stop
      end if
!      
!     given ntab EVEN zones, work out ntab+1 ODD zones for use in Filon's method
!
      ntabp1 = ntab + 1
      xtab_e(1:ntabp1) = x3a(ks:ke+1)
      h = dx3a(ks)
!
!     Filon parameters
!
      theta = t * h
      sint = sin ( theta )
      cost = cos ( theta )
      
      if(theta.gt.tol) then
         alpha = ( theta**2 + theta * sint * cost
     &        - 2d0 * sint**2 ) / theta**3
         
         beta = ( 2d0 * theta + 2d0 *theta* cost**2 
     &        -4d0* sint * cost ) / theta**3
         
         gamma = 4d0 * ( sint - theta * cost ) /theta**3
      else
         alpha = a3*theta**3 + a5*theta**5 + a7*theta**7
     
         beta = b0 + b2*theta**2 + b4*theta**4 
     &        + b6*theta**6 + b8*theta**8
         
         gamma = g0  + g2*theta**2 + g4*theta**4
     &        + g6*theta**6 
 
      endif   
!
!     cosine integral. first extrapolate to get function at innermost and outermost cell EDGE.
!     interior averages are arithmetic averages of the cell centered stuff
!
!      df = -3d0*ftab_cos(5) + 16d0*ftab_cos(4)-36d0*ftab_cos(3)
!     &     +48d0*ftab_cos(2) - 25d0*ftab_cos(1)
!      df = df/(12d0*h)
      df = -ftab_cos(3)
     &     + 4d0*ftab_cos(2) - 3d0*ftab_cos(1)
      df = df/(2d0*h)
      f0  = ftab_cos(1) - (h/2d0)*df
      
!      df = 3d0*ftab_cos(ntab-4) - 16d0*ftab_cos(ntab-3)
!     &     +36d0*ftab_cos(ntab-2)
!     &     -48d0*ftab_cos(ntab-1) + 25d0*ftab_cos(ntab)
!      df = df/(12d0*h)
      df = 
     &     ftab_cos(ntab-2)
     &     -4d0*ftab_cos(ntab-1) + 3d0*ftab_cos(ntab)
      df = df/(2d0*h)
      fnp1= ftab_cos(ntab) + (h/2d0)*df
      
      ftab_ce(1)    =  f0
      ftab_ce(2:ntab) = (ftab_cos(1:ntab-1) + ftab_cos(2:ntab))/2d0
      ftab_ce(ntabp1)=  fnp1
      
      c2n = 
     &     sum ( ftab_ce(1:ntabp1:2) * cos ( t * xtab_e(1:ntabp1:2) ) )
     &     -0.5d0 * ( ftab_ce(ntabp1) * cos ( t * xtab_e(ntabp1) )     
     &     + ftab_ce(1) * cos ( t * xtab_e(1) ) )
      
!      c2nm1 = sum ( ftab_ce(2:ntabp1-1:2) *cos (t*xtab2(2:ntab2-1:2)))
      c2nm1 = sum ( ftab_ce(2:ntab:2) *cos (t*xtab_e(2:ntab:2)))
      
      result_cos =  h*( 
     &     alpha * ( ftab_ce(ntabp1) * sin ( t * xtab_e(ntabp1) ) 
     &     -ftab_ce(1) * sin ( t * xtab_e(1) ) ) 
     &     + beta * c2n 
     &     + gamma * c2nm1 )
!
!     repeat for sine integral
!
!      df = -3d0*ftab_sin(5) + 16d0*ftab_sin(4)-36d0*ftab_sin(3)
!     &     +48d0*ftab_sin(2) - 25d0*ftab_sin(1)
!      df = df/(12d0*h)
      df = -ftab_sin(3)
     &     + 4d0*ftab_sin(2) - 3d0*ftab_sin(1)
      df = df/(2d0*h)
      f0  = ftab_sin(1) - (h/2d0)*df
      
!      df = 3d0*ftab_sin(ntab-4) - 16d0*ftab_sin(ntab-3)
!     &     +36d0*ftab_sin(ntab-2)
!     &     -48d0*ftab_sin(ntab-1) + 25d0*ftab_sin(ntab)
!      df = df/(12d0*h)
      df = 
     &     ftab_sin(ntab-2)
     &     -4d0*ftab_sin(ntab-1) + 3d0*ftab_sin(ntab)
      df = df/(2d0*h)
      fnp1= ftab_sin(ntab) + (h/2d0)*df
      ftab_se(1)    =  f0
      ftab_se(2:ntab) = (ftab_sin(1:ntab-1) + ftab_sin(2:ntab))/2d0
      ftab_se(ntabp1)=  fnp1
      s2n = 
     &     sum ( ftab_se(1:ntabp1:2) * sin ( t * xtab_e(1:ntabp1:2) ) )
     &     - 0.5d0 * ( ftab_se(ntabp1) * sin ( t * xtab_e(ntabp1) )    
     &     + ftab_se(1) * sin ( t * xtab_e(1) ) )
      
!      s2nm1 = sum ( ftab_se(2:ntabp1-1:2) * 
!     &     sin ( t * xtab2(2:ntabp1-1:2)))
      
      s2nm1 = sum ( ftab_se(2:ntab:2) * 
     &     sin ( t * xtab_e(2:ntab:2)))
    
      result_sin = h * ( 
     &     alpha * ( ftab_se(1) * cos ( t * xtab_e(1) ) 
     &     - ftab_se(ntabp1) * cos ( t * xtab_e(ntabp1) ) ) 
     &     + beta * s2n 
     &     + gamma * s2nm1 )
      return
      end
      
      subroutine flinn(ftab_cos,ftab_sin,ntab,t,result_cos,
     &     result_sin )
      use grid
      use mpiyes
      use mpipar
!
!***********************************************************************
!
!!Flinn's method (1960) is a modified Filon scheme
!
!
!  Discussion:
!
!    The integral to be approximated has the form:
!
!      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
!
!    where T is user specified.
!    Output, real RESULT, the approximate value of the integral.
!
      implicit none
!
      integer ntab, i, ntabp1
!
      real*8, parameter :: tol=0.9d0
      real*8, parameter :: N0=16d0/15d0, N2=-8d0/105d0,
     &       N4=2d0/945d0, N6=-1d0/31185d0, N8=1d0/3243240d0,
     &       R0=14d0/15d0, R2=-16d0/105d0, R4=22d0/945d0,
     &       R6=-608d0/311850d0, R8=268d0/2837835d0,
     &       S1=19d0/105d0, S3=-2d0/63d0, S5=1d0/275d0,
     &       S7=-2d0/8775d0, S9=34d0/3869775d0
      real*8 :: a, alpha, b, beta, c2n,
     &     c2nm1, cost, gamma, h, s2n, s2nm1,
     &     result_cos, result_sin, sint, t, theta,
     &     N, R, S, f0, fnp1, df,
     &     t2, t4, t6, t8
      real*8 :: ftab_cos(ntab), ftab_sin(ntab), xtab(ntab),
     &     ftab_ce(ntab+1), ftab_se(ntab+1), xtab_e(ntab+1)
!
      if ( a == b ) then
         result_cos = 0.0d0
         result_sin = 0.0d0 
         return
      end if
 
!      if ( ntab <= 1 ) then
!         write ( *, '(a)' ) ' '
!         write ( *, '(a)' ) 'FILON_COS - Fatal error!'
!         write ( *, '(a)' ) '  NTAB < 2'
!         write ( *, '(a,i6)' ) '  NTAB = ', ntab
!         stop
!      end if
 
!      if ( mod ( ntab, 2 ) /= 1 ) then
!         write ( *, '(a)' ) ' '
!         write ( *, '(a)' ) 'FILON_COS - Fatal error!'
!         write ( *, '(a)' ) '  NTAB must be odd.'
!         write ( *, '(a,i6)' ) '  NTAB = ', ntab
!         stop
!      end if
      if ( mod ( ntab, 2 ).ne. 0 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'FILON_COS - Fatal error!'
         write ( *, '(a)' ) '  NTAB must be odd.'
         write ( *, '(a,i6)' ) '  NTAB = ', ntab
         stop
      end if
!      
!     given ntab EVEN zones, work out ntab+1 ODD zones for use in Flinn's method
!
      ntabp1 = ntab + 1
      xtab_e(1:ntabp1) = x3a(ks:ke+1)
      h = dx3a(ks)
!
!     Flinn parameters. we're going to ignore terms proportional to function gradients. 
!     these have a h^2 dependence
      theta = t * h
      if(theta.gt.tol) then
         print*, 'theta too large for expansion'
         stop
      endif
      sint = sin ( theta )
      cost = cos ( theta )
      t2 = theta**2
      t4 = theta**4
      t6 = theta**6
      t8 = theta**8
      N = N0 + N2*t2 + N4*t4
     &     +N6*t6 + N8*t8
      
      R = R0 + R2*t2 + R4*t4 +R6*t6
     &     + R8*t8
      S = S1*theta + S3*theta**3 + S5*theta**5
     &     +S7*theta**7 + S9*theta**9
      
!
!     cosine integral. first extrapolate to get function at innermost and outermost cell EDGE.
!     interior averages are arithmetic averages of the cell centered stuff. this approximation saves communication costs
!
!      df = -3d0*ftab_cos(5) + 16d0*ftab_cos(4)-36d0*ftab_cos(3)
!     &     +48d0*ftab_cos(2) - 25d0*ftab_cos(1)
!      df = df/(12d0*h)
      df = -ftab_cos(3)
     &     + 4d0*ftab_cos(2) - 3d0*ftab_cos(1)
      df = df/(2d0*h)
      f0  = ftab_cos(1) - (h/2d0)*df
      
!      df = 3d0*ftab_cos(ntab-4) - 16d0*ftab_cos(ntab-3)
!     &     +36d0*ftab_cos(ntab-2)
!     &     -48d0*ftab_cos(ntab-1) + 25d0*ftab_cos(ntab)
!      df = df/(12d0*h)
      df = 
     &     ftab_cos(ntab-2)
     &     -4d0*ftab_cos(ntab-1) + 3d0*ftab_cos(ntab)
      df = df/(2d0*h)
      fnp1= ftab_cos(ntab) + (h/2d0)*df
      ftab_ce(1)    =  f0
      ftab_ce(2:ntab) = (ftab_cos(1:ntab-1) + ftab_cos(2:ntab))/2d0
      ftab_ce(ntabp1)=  fnp1
      c2n = 
     &     sum ( ftab_ce(1:ntabp1:2) * cos ( t * xtab_e(1:ntabp1:2) ) ) 
     &     -0.5d0 * ( ftab_ce(ntabp1) * cos ( t * xtab_e(ntabp1) )     
     &     + ftab_ce(1) * cos ( t * xtab_e(1) ) )
      
      c2nm1 = sum ( ftab_ce(2:ntab:2) *cos (t*xtab_e(2:ntab:2)))
      
      result_cos =  h*( 
     &     S * ( ftab_ce(ntabp1) * sin ( t * xtab_e(ntabp1) ) 
     &     -ftab_ce(1) * sin ( t * xtab_e(1) ) ) 
     &     + R * c2n 
     &     + N * c2nm1 )
!
!     repeat for sine integral
!
      
!      df = -3d0*ftab_sin(5) + 16d0*ftab_sin(4)-36d0*ftab_sin(3)
!     &     +48d0*ftab_sin(2) - 25d0*ftab_sin(1)
!      df = df/(12d0*h)
      df = -ftab_sin(3)
     &     + 4d0*ftab_sin(2) - 3d0*ftab_sin(1)
      df = df/(2d0*h)
      f0  = ftab_sin(1) - (h/2d0)*df
      
!      df = 3d0*ftab_sin(ntab-4) - 16d0*ftab_sin(ntab-3)
!     &     +36d0*ftab_sin(ntab-2)
!     &     -48d0*ftab_sin(ntab-1) + 25d0*ftab_sin(ntab)
!      df = df/(12d0*h)
      df = 
     &     ftab_sin(ntab-2)
     &     -4d0*ftab_sin(ntab-1) + 3d0*ftab_sin(ntab)
      df = df/(2d0*h)
      fnp1= ftab_sin(ntab) + (h/2d0)*df
      ftab_se(1)    =  f0
      ftab_se(2:ntab) = (ftab_sin(1:ntab-1) + ftab_sin(2:ntab))/2d0
      ftab_se(ntabp1)=  fnp1
      s2n = 
     &     sum ( ftab_se(1:ntabp1:2) * sin ( t * xtab_e(1:ntabp1:2) ) ) 
     &     - 0.5d0 * ( ftab_se(ntabp1) * sin ( t * xtab_e(ntabp1) )    
     &     + ftab_se(1) * sin ( t * xtab_e(1) ) )
      
      s2nm1 = sum ( ftab_se(2:ntab:2) * sin ( t * xtab_e(2:ntab:2)))
      
      result_sin = h * ( 
     &     S * ( ftab_se(1) * cos ( t * xtab_e(1) ) 
     &     - ftab_se(ntabp1) * cos ( t * xtab_e(ntabp1) ) ) 
     &     + R * s2n 
     &     + N * s2nm1 )
  
      return
      end
      
      
      
      subroutine rvec_even ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Output, real A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
!
      implicit none
!
      integer n
!
      real*8 :: a(n)
      real*8 :: ahi
      real*8 :: alo
      integer i
!     
      if ( n == 1 ) then
         
         a(1) = 0.5E+00 * ( alo + ahi )
         
      else
         
         do i = 1, n
            a(i) = ( dble ( n - i ) * alo + dble ( i - 1 ) * ahi ) 
     &           / dble( n - 1 )
         end do
         
      end if
      
      return
      end
      subroutine eightP_newton_cotes(integrand, size, dx, tot)
      
      implicit none
      real*8, parameter :: c1=751.d0, c2=3577.d0, c3=1323.d0, 
     .     c4=2989.d0, coeff= 7.d0/17280.d0
      real*8 :: caxis(8)
      integer :: size, nseg ! size must be divisible by 8 
      integer :: i, j, ibeg
      real*8  :: integrand(size), sample(8)
      real*8  ::dx, tot
      if(mod(size,8).ne.0) then
         print*, '8-point newton cotes need multiples of 8 zones'
         stop
      endif
      caxis(1:8) = (/c1, c2, c3, c4, c4, c3, c2, c1/)*coeff
      nseg = size/8
      
      tot = 0.d0
      do i=1, nseg
         ibeg = (i-1)*8 + 1
         sample(1:8) = integrand(ibeg:ibeg+7)
         tot = tot +
     &        sum(sample(1:8)*caxis(1:8))*dx
      enddo
      
      return
      end
      
      subroutine eightP_newton_cotes_v2(integrand, size, dx, tot)
      
      implicit none
      real*8, parameter :: c1=3d0/8d0, c2=9d0/8d0, c3=9d0/8d0, 
     .     c4=17d0/24d0, c5=4d0/3d0, c6=2d0/3d0, c7=4d0/3d0,
     .     c8=1d0/3d0
      real*8 :: caxis(8)
      integer :: size, nseg ! size must be divisible by 8 
      integer :: i, j, ibeg
      real*8  :: integrand(size), sample(8)
      real*8  ::dx, tot
      if(mod(size,8).ne.0) then
         print*, '8-point newton cotes need multiples of 8 zones'
         stop
      endif
      caxis(1:8) = (/c1, c2, c3, c4, c5, c6, c7, c8/)*dx
      nseg = size/8
      
      do i=1, nseg
         ibeg = (i-1)*8 + 1
         integrand(ibeg:ibeg+7)
     &        = integrand(ibeg:ibeg+7)*caxis(1:8)
      enddo
      
      tot = sum(integrand)
      return
      end
      subroutine fourP_newton_cotes(integrand, size, dx, tot)
      
      implicit none
      real*8, parameter :: c1=1.0d0, c2=3.0d0, c3=3.d0, 
     .     c4=1.d0, coeff= 3.d0/8.d0
      real*8 :: caxis(4)
      integer :: size, nseg ! size must be divisible by 4 
      integer :: i, j, ibeg
      real*8  :: integrand(size)
      real*8  ::dx, tot
      if(mod(size,4).ne.0) then
         print*, '4-point newton cotes need multiples of 8 zones'
         stop
      endif
      caxis(1:4) = (/c1, c2, c3, c4/)*coeff*dx
      nseg = size/4
      
      do i=1, nseg
         ibeg = (i-1)*4 + 1
         integrand(ibeg:ibeg+3) = integrand(ibeg:ibeg+3)*caxis(1:4)
      enddo
      
      tot = sum(integrand)
      return
      end
      subroutine fourP_newton_cotes_open(integrand, size, dx, tot)
      
      implicit none
      real*8, parameter :: c1=11.d0, c2=1.0d0, c3=1.d0, 
     .     c4=11.d0, coeff= 5.d0/24.d0
      real*8 :: caxis(4)
      integer :: size, nseg ! size must be divisible by 4 
      integer :: i, j, ibeg
      real*8  :: integrand(size)
      real*8  ::dx, tot
      if(mod(size,4).ne.0) then
         print*, '4-point newton cotes need multiples of 8 zones'
         stop
      endif
      caxis(1:4) = (/c1, c2, c3, c4/)*coeff*dx
      nseg = size/4
      
      do i=1, nseg
         ibeg = (i-1)*4 + 1
         integrand(ibeg:ibeg+3) = integrand(ibeg:ibeg+3)*caxis(1:4)
      enddo
      
      tot = sum(integrand)
      return
      end
      
 
      subroutine cspint ( ftab, xtab, ntab, a, b, y, e, work, result )
!
!***********************************************************************
!
!! CSPINT estimates the integral of a tabulated function.
!
!
!
      implicit none
!     
      integer ntab
!     
      real*8 ::  a
      real*8 ::  b
      real*8 ::  e(ntab)
      real*8 ::  ftab(ntab)
      integer i
      integer j
      real*8 ::  r
      real*8 ::  result
      real*8 ::  s
      real*8 ::  term
      real*8 ::  u
      real*8 ::  work(ntab)
      real*8 ::  xtab(ntab)
      real*8 ::  y(3,ntab)
!     
      if ( ntab < 3 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'CSPINT - Fatal error!'
         write ( *, '(a,i6)' ) 
     &        '  NTAB must be at least 3, but input NTAB = ',ntab
         stop
      endif
      
      do i = 1, ntab-1
         if ( xtab(i+1) <= xtab(i) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'CSPINT - Fatal error!'
            write ( *, '(a)' ) '  Nodes not in strict increasing order.'
            write ( *, '(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
            write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
            write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
            stop
         end if
      enddo
      
      s = 0.d0
      do i = 1, ntab-1
         r = ( ftab(i+1) - ftab(i) ) / ( xtab(i+1) - xtab(i) )
         y(2,i) = r - s
         s = r
      end do
      
      result = 0.0d0
      s = 0.0d0
      r = 0.0d0
      y(2,1) = 0.0d0
      y(2,ntab) = 0.0d0
      
      do i = 2, ntab-1
         y(2,i) = y(2,i)+r*y(2,i-1)
         work(i) = 2.0d0 * ( xtab(i-1) - xtab(i+1) ) - r * s
         s = xtab(i+1) - xtab(i)
         r = s / work(i)
      end do
      
      do j = 2, ntab-1
         i = ntab+1-j
         y(2,i) = ((xtab(i+1)-xtab(i))*y(2,i+1)-y(2,i)) / work(i)
      end do
      
      do i = 1, ntab-1
         s = xtab(i+1)-xtab(i)
         r = y(2,i+1)-y(2,i)
         y(3,i) = r / s
         y(2,i) = 3.0d0 * y(2,i)
         y(1,i) = (ftab(i+1)-ftab(i)) / s-(y(2,i)+r)*s
      end do
      
      e(1) = 0.0d0
      do i = 1, ntab-1
         s = xtab(i+1)-xtab(i)
         term = (((y(3,i)* 0.25d0 *s+y(2,i) / 3d0 ) *s+y(1,i)* 0.5d0 )*s
     &        +ftab(i))*s
         e(i+1) = e(i) + term
      end do
      
                                !
!     Determine where the endpoints A and B lie in the mesh of XTAB's.
!     
      r = a
      u = 1.0d0
      
      do j = 1, 2
!     
!     The endpoint is less than or equal to XTAB(1).
!     
         if ( r <= xtab(1) ) then
            result = result-u*((r-xtab(1))*y(1,1)*0.5d0 +ftab(1))
     &           *(r-xtab(1))
!     
!     The endpoint is greater than or equal to XTAB(NTAB).
!     
         else if ( r >= xtab(ntab) ) then
            
            result = result-u*(e(ntab)+(r-xtab(ntab))*(ftab(ntab)+ 
     &       0.5d0 *(ftab(ntab-1)+(xtab(ntab)-xtab(ntab-1))*y(1,ntab-1)) 
     &           *(r-xtab(ntab))))
!     
!     The endpoint is strictly between XTAB(1) and XTAB(NTAB).
!     
         else
            do i = 1,ntab-1
               if ( r <= xtab(i+1) ) then
                  r = r-xtab(i)
                  result = result-u*(e(i)+(((y(3,i)*0.25d0*r
     &                 +y(2,i)/3d0)*r 
     &                 +y(1,i)*0.5d0 )*r+ftab(i))*r)
                  go to 120
               end if
            end do
         end if
 120     continue
         
         u = -1.0d0
         r = b
         
      end do
      
      return
      end
      subroutine hiordq ( n, y, delt, work, result )
!
!***********************************************************************
!
!! HIORDQ approximates the integral of a function using equally spaced data.
!
!
!  Discussion:
!
!    The method applies the trapezoidal rule to various subsets of the
!    data, and then applies Richardson extrapolation.
!
!
      implicit none
!     
      integer n
!     
      real*8 ::  delt
      real*8 ::  fac
      integer i
      integer j
      integer jbak
      integer jj
      integer k
      real*8 ::  result
      real*8 ::  sum2
      real*8 ::  sum1
      real*8 ::  work(2*(n-1))
      real*8 ::  y(n)
!     
!     Determine initial trapezoidal rule
!     
      sum1 = ( y(1) + y(n) ) / 2.0d0
      j = -1
      
      do k = 1, n-1
!     
!     Check if K divides N-1
!     
         if ( ((n-1)/k)*k == n-1 ) then
!     
!     Determine the K-point trapezoidal rule.
!     
            sum2 = -sum1
            do i = 1, n, (n-1)/k
               sum2 = sum2 + y(i)
            end do
            
            j = j + 2
            work(j) = delt * sum2 * dble( ( n - 1 ) / k )
            work(j+1) = dble( ((n-1)/k)**2 )
!     
!     Apply Richardson extrapolation.
!     
            if ( k /= 1 ) then
               
               do jj = 3, j, 2
                  jbak = j+1-jj
                  fac = work(j+1) / ( work(j+1) - work(jbak+1) )
                  work(jbak) = work(jbak+2) + fac * ( work(jbak) 
     &                 - work(jbak+2) )
               end do
               
            end if
            
         end if
         
      end do
      
      result = work(1)
      
      return
      end
      subroutine simpsn ( h, y, num, result )
!
!***********************************************************************
!     
!     ! SIMPSN approximates the integral of evenly spaced data.
!     
!     
!     Discussion:
!     
!     Simpson's rule is used.
!     
!     
      implicit none
!     
      integer num
!     
      real*8 ::  del(3)
      real*8 ::  f
      real*8 ::  g(3)
      real*8 ::  h
      integer i
      integer n
      real*8 ::  pii(3)
      real*8 ::  result
      real*8 ::  sum1
      real*8 ::  y(num)
!     
      result = 0.0d0
      
      if ( num <= 2 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SIMPSN - Fatal error!'
         write ( *, '(a,i6)' ) '  NUM < 2, NUM = ', num
         stop
      end if
      
      if ( mod ( num, 2 ) == 0 ) then
         n = num-1
      else
         n = num
      end if
      
      result = y(1) + y(n) + 4.0d0 * y(n-1)
      do i = 2, n-2, 2
         result = result + 4.0d0 * y(i) + 2.0d0 * y(i+1)
      end do
      result = h * result / 3.0d0
      
      if ( mod(num,2) == 1 ) then
         return
      end if
      
      f = h*h*h
      del(1) = h
      del(2) = -2.0d0 * h
      del(3) = h
      g(1) = h
      g(2) = 0.0d0
      g(3) = -h
      pii(1) = 0.0d0
      pii(2) = -h*h
      pii(3) = 0.0d0
      n = n-1
      
      sum1 = 0.0d0
      do i = 1, 3
         sum1 = sum1 + y(n-1+i) * del(i) * 
     &        ( f / 3.0d0 - g(i) * 0.5d0 * h * h + pii(i) * h )
      enddo
      
      result = result + 0.5d0 * sum1 / h**3
      
      return
      end
C      *DECK DAVINT
      SUBROUTINE DAVINT (X, Y, N, XLO, XUP, ANS, IERR)
C***BEGIN PROLOGUE  DAVINT
C***PURPOSE  Integrate a function tabulated at arbitrarily spaced
C            abscissas using overlapping parabolas.
C***LIBRARY   SLATEC
C***CATEGORY  H2A1B2
C***TYPE      DOUBLE PRECISION (AVINT-S, DAVINT-D)
C***KEYWORDS  INTEGRATION, QUADRATURE, TABULATED DATA
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C         DAVINT integrates a function tabulated at arbitrarily spaced
C         abscissas.  The limits of integration need not coincide
C         with the tabulated abscissas.
C
C         A method of overlapping parabolas fitted to the data is used
C         provided that there are at least 3 abscissas between the
C         limits of integration.  DAVINT also handles two special cases.
C         If the limits of integration are equal, DAVINT returns a
C         result of zero regardless of the number of tabulated values.
C         If there are only two function values, DAVINT uses the
C         trapezoid rule.
C
C     Description of Parameters
C         The user must dimension all arrays appearing in the call list
C              X(N), Y(N)
C
C         Input--
C      X    - DOUBLE PRECISION array of abscissas, which must be in
C             increasing order.
C      Y    - DOUBLE PRECISION array of function values. i.e.,
C                Y(I)=FUNC(X(I))
C      N    - The integer number of function values supplied.
C                N .GE. 2 unless XLO = XUP.
C      XLO  - DOUBLE PRECISION lower limit of integration
C      XUP  - DOUBLE PRECISION upper limit of integration.  Must have
C              XLO.LE.XUP
C
C         Output--
C      ANS  - Double Precision computed approximate value of integral
C      IERR - A status code
C           --Normal Code
C                =1 Means the requested integration was performed.
C           --Abnormal Codes
C                =2 Means XUP was less than XLO.
C                =3 Means the number of X(I) between XLO and XUP
C                   (inclusive) was less than 3 and neither of the two
C                   special cases described in the abstract occurred.
C                   No integration was performed.
C                =4 Means the restriction X(I+1).GT.X(I) was violated.
C                =5 Means the number N of function values was .lt. 2.
C                   ANS is set to zero if IERR=2,3,4,or 5.
C
C    DAVINT is documented completely in SC-M-69-335
C    Original program from *Numerical Integration* by Davis & Rabinowitz
C    Adaptation and modifications by Rondall E Jones.
C
C***REFERENCES  R. E. Jones, Approximate integrator of functions
C                 tabulated at arbitrarily spaced abscissas,
C                 Report SC-M-69-335, Sandia Laboratories, 1969.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   690901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAVINT
C
      INTEGER I, IERR, INLFT, INRT, ISTART, ISTOP, N
      DOUBLE PRECISION A, ANS, B, C, CA, CB, CC, FL, FR, R3, RP5,
     1     SLOPE, SUM, SYL, SYL2, SYL3, SYU, SYU2, SYU3, TERM1, TERM2,
     2     TERM3, X, X1, X12, X13, X2, X23, X3, XLO, XUP, Y
      DIMENSION X(*),Y(*)
C     BEGIN BLOCK PERMITTING ...EXITS TO 190
C     BEGIN BLOCK PERMITTING ...EXITS TO 180
C***  FIRST EXECUTABLE STATEMENT  DAVINT
      IERR = 1
      ANS = 0.0D0
      IF (XLO .GT. XUP) GO TO 160
      IF (XLO .EQ. XUP) GO TO 150
      IF (N .GE. 2) GO TO 10
      IERR = 5
      print*, 'SLATEC', 'DAVINT',
     &     'LESS THAN TWO FUNCTION VALUES WERE SUPPLIED.'
      stop
C     ...............EXIT
      GO TO 190
 10   CONTINUE
      DO 20 I = 2, N
C     ............EXIT
         IF (X(I) .LE. X(I-1)) GO TO 180
C     ...EXIT
         IF (X(I) .GT. XUP) GO TO 30
 20   CONTINUE
 30   CONTINUE
      IF (N .GE. 3) GO TO 40
C     
C     SPECIAL N=2 CASE
      SLOPE = (Y(2) - Y(1))/(X(2) - X(1))
      FL = Y(1) + SLOPE*(XLO - X(1))
      FR = Y(2) + SLOPE*(XUP - X(2))
      ANS = 0.5D0*(FL + FR)*(XUP - XLO)
C     ...............EXIT
      GO TO 190
 40   CONTINUE
      IF (X(N-2) .GE. XLO) GO TO 50
      IERR = 3
      print*, 'SLATEC', 'DAVINT',
     &     'THERE WERE LESS THAN THREE FUNCTION VALUES ',
     &     'BETWEEN THE LIMITS OF INTEGRATION.'
      stop
C     ...............EXIT
      GO TO 190
 50   CONTINUE
      IF (X(3) .LE. XUP) GO TO 60
      IERR = 3
      print*, 'SLATEC', 'DAVINT',
     &     'THERE WERE LESS THAN THREE FUNCTION VALUES ',
     &     'BETWEEN THE LIMITS OF INTEGRATION.'
      stop
C     ...............EXIT
      GO TO 190
 60   CONTINUE
      I = 1
 70   IF (X(I) .GE. XLO) GO TO 80
      I = I + 1
      GO TO 70
 80   CONTINUE
      INLFT = I
      I = N
 90   IF (X(I) .LE. XUP) GO TO 100
      I = I - 1
      GO TO 90
 100  CONTINUE
      INRT = I
      IF ((INRT - INLFT) .GE. 2) GO TO 110
      IERR = 3
      print*, 'SLATEC', 'DAVINT',
     &     'THERE WERE LESS THAN THREE FUNCTION VALUES ',
     &     'BETWEEN THE LIMITS OF INTEGRATION.'
      stop
C     ...............EXIT
      GO TO 190
 110  CONTINUE
      ISTART = INLFT
      IF (INLFT .EQ. 1) ISTART = 2
      ISTOP = INRT
      IF (INRT .EQ. N) ISTOP = N - 1
C     
      R3 = 3.0D0
      RP5 = 0.5D0
      SUM = 0.0D0
      SYL = XLO
      SYL2 = SYL*SYL
      SYL3 = SYL2*SYL
C     
      DO 140 I = ISTART, ISTOP
         X1 = X(I-1)
         X2 = X(I)
         X3 = X(I+1)
         X12 = X1 - X2
         X13 = X1 - X3
         X23 = X2 - X3
         TERM1 = Y(I-1)/(X12*X13)
         TERM2 = -Y(I)/(X12*X23)
         TERM3 = Y(I+1)/(X13*X23)
         A = TERM1 + TERM2 + TERM3
         B = -(X2 + X3)*TERM1 - (X1 + X3)*TERM2
     1        - (X1 + X2)*TERM3
         C = X2*X3*TERM1 + X1*X3*TERM2 + X1*X2*TERM3
         IF (I .GT. ISTART) GO TO 120
         CA = A
         CB = B
         CC = C
         GO TO 130
 120     CONTINUE
         CA = 0.5D0*(A + CA)
         CB = 0.5D0*(B + CB)
         CC = 0.5D0*(C + CC)
 130     CONTINUE
         SYU = X2
         SYU2 = SYU*SYU
         SYU3 = SYU2*SYU
         SUM = SUM + CA*(SYU3 - SYL3)/R3
     1        + CB*RP5*(SYU2 - SYL2) + CC*(SYU - SYL)
         CA = A
         CB = B
         CC = C
         SYL = SYU
         SYL2 = SYU2
         SYL3 = SYU3
 140  CONTINUE
      SYU = XUP
      ANS = SUM + CA*(SYU**3 - SYL3)/R3
     1     + CB*RP5*(SYU**2 - SYL2) + CC*(SYU - SYL)
 150  CONTINUE
      GO TO 170
 160  CONTINUE
      IERR = 2
      print*, 'SLATEC', 'DAVINT',
     &     'THE UPPER LIMIT OF INTEGRATION WAS NOT GREATER ',
     &     'THAN THE LOWER LIMIT.'
      stop
 170  CONTINUE
C     ......EXIT
      GO TO 190
 180  CONTINUE
      IERR = 4
      print*, 'SLATEC', 'DAVINT',
     &     'THE ABSCISSAS WERE NOT STRICTLY INCREASING.'
      stop
 190  CONTINUE
      RETURN
      END
      
      subroutine plint ( ftab, xtab, ntab, a, b, result )
!     
!***********************************************************************
!     
!     ! PLINT approximates the integral of unequally spaced data.
!     
!     
!     Discussion:
!     
!     The method uses piecewise linear interpolation.
!     
!     Reference:
!     
!     Philip Davis and Philip Rabinowitz,
!     Methods of Numerical Integration,
!     Blaisdell Publishing, 1967.
!     
!     Modified:
!     
!     30 October 2000
!     
!     Parameters:
!     
!    Input, real FTAB(NTAB), the function values, FTAB(I) = F(XTAB(I)).
!     
!     Input, real XTAB(NTAB), the abscissas at which the
!     function values are given.  The XTAB's must be distinct
!     and in ascending order.
!     
!     Input, integer NTAB, the number of entries in FTAB and
!     XTAB.  NTAB must be at least 2.
!     
!     Input, real A, the lower limit of integration.  A should
!     be, but need not be, near one endpoint of the interval
!     (X(1), X(NTAB)).
!     
!     Input, real B, the upper limit of integration.  B should
!     be, but need not be, near one endpoint of the interval
!     (X(1), X(NTAB)).
!     
!     Output, real RESULT, the approximate value of the integral.
!     
      implicit none
!     
      integer ntab
!     
      real*8 :: a
      real*8 :: b
      real*8 :: fa
      real*8 :: fb
      real*8 :: ftab(ntab)
      integer i
      integer ihi
      integer ilo
      integer ind
      real*8 :: result
      real*8 :: slope
      real*8 :: syl
      real*8 :: xtab(ntab)
!     
      if ( a == b ) then
         result = 0.0d0
         return
      end if
      
      if ( ntab < 2 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'PLINT - Fatal error!'
         write ( *, '(a,i6)' ) '  NTAB < 2, NTAB = ',ntab
         stop
      end if
      
      do i = 2, ntab
         if ( xtab(i) <= xtab(i-1) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PLINT - Fatal error!'
            write ( *, '(a)' ) '  Nodes not in strict increasing order.'
            write ( *, '(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
            write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
            write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
            stop
         end if
      end do
!     
!     If A > B, temporarily switch A and B, and store sign.
!     
      if ( a > b ) then
         syl = b
         b = a
         a = syl
         ind = -1
      else
         syl = a
         ind = 1
      end if
!     
!     Find ILO and IHI so that A <= XTAB(ILO) <= XTAB(IHI) <= B
!     with the possible exception that A and B may be in the same
!     interval, or completely to the right or left of the XTAB's.
!     
      ilo = 1
      ihi = ntab
      do i = 1, ntab
         if ( a <= xtab(i) ) then
            exit
         end if
         ilo = ilo+1
      end do
      
      do i = 1, ntab
         if ( b >= xtab(i) ) then
            exit
         end if
         ihi = ihi-1
      end do
!     
!     Treat special cases where A, B lie both to left or both to right
!     of XTAB interval, or inbetween same pair of XTAB's.
!     
      if ( ihi == 0 ) then
         slope = (ftab(2)-ftab(1))/(xtab(2)-xtab(1))
         fa = ftab(1) + slope*(a-xtab(1))
         fb = ftab(1) + slope*(b-xtab(1))
         result = 0.5d0 * (b-a) * (fa+fb)
         go to 110
      else if ( ilo == ntab+1 ) then
         slope = (ftab(ntab)-ftab(ntab-1))/(xtab(ntab)-xtab(ntab-1))
         fa = ftab(ntab-1)+slope*(a-xtab(ntab-1))
         fb = ftab(ntab-1)+slope*(b-xtab(ntab-1))
         result = 0.5d0 * (b-a) * (fa+fb)
         go to 110
      else if ( ihi+1 == ilo ) then
         slope = (ftab(ilo)-ftab(ihi))/(xtab(ilo)-xtab(ihi))
         fa = ftab(ihi)+slope*(a-xtab(ihi))
         fb = ftab(ihi)+slope*(b-xtab(ihi))
         result = 0.5d0 * (b-a) * (fa+fb)
         go to 110
      end if
!     
!     Carry out approximate integration.  We know that ILO is no greater
!     than IHI-1, but equality is possible; A and B may be on either side
!     of a single XTAB(I).  That's OK, then the loop below won't be executed
!     at all.
!     
      result = 0.0d0
      do i = ilo, ihi-1
         result = result + 0.5d0*(xtab(i+1)-xtab(i))*(ftab(i)+ftab(i+1))
      end do
!     
!     Add contribution from A-ILO and IHI-B.
!     Still have to watch out if ILO = 1 or IHI=NTAB...
!     
      if ( ilo == 1 ) then
         slope = (ftab(2)-ftab(1)) / (xtab(2)-xtab(1))
         fa = ftab(1) + slope*(a-xtab(1))
         result = result + 0.5d0 * (xtab(ilo)-a)*(fa+ftab(ilo))
      else
         slope = (ftab(ilo)-ftab(ilo-1)) / (xtab(ilo)-xtab(ilo-1))
         fa = ftab(ilo-1) + slope*(a-xtab(ilo-1))
         result = result + 0.5d0 * (xtab(ilo)-a)*(fa+ftab(ilo))
      end if
      
      if ( ihi == ntab ) then
         slope = (ftab(ntab)-ftab(ntab-1)) / (xtab(ntab)-xtab(ntab-1))
         fb = ftab(ntab-1) + slope*(b-xtab(ntab-1))
         result = result + 0.5d0*(b-xtab(ntab))*(fb+ftab(ntab))
      else
         slope = (ftab(ihi+1)-ftab(ihi)) / (xtab(ihi+1)-xtab(ihi))
         fb = ftab(ihi) + slope*(b-xtab(ihi))
         result = result + 0.5d0*(b-xtab(ihi))*(fb+ftab(ihi))
      end if
!     
!     Restore original values of A and B, reverse sign of integral
!     because of earlier switch.
!     
 110  continue
      
      if ( ind /= 1 ) then
         ind = 1
         syl = b
         b = a
         a = syl
         result = -result
      end if
      
      return
      end
  
      subroutine simpne ( x, y, num, result )
!     
!***********************************************************************
!     
!     ! SIMPNE approximates the integral of unevenly spaced data.
!     
!     
!     Discussion:
!     
!     The routine repeatedly interpolates a 3-point Lagrangian polynomial 
!     to the data and integrates that exactly.
!     
!     Reference:
!     
!     Philip Davis and Philip Rabinowitz,
!     Methods of Numerical Integration,
!     Blaisdell Publishing, 1967.
!     
!     Modified:
!     
!     30 October 2000
!     
!     Parameters:
!     
!     Input, real*8 :: X(NUM), contains the X values of the data, in order.
!     
!    Input, real*8 :: Y(NUM), contains the Y values of the data.
!
!    Input, integer NUM, number of data points.  NUM must be at least 3.
!
!    Output, real*8 :: RESULT.
!    RESULT is the approximate value of the integral.
!
      implicit none
!     
      integer num
!     
      real*8 :: del(3)
      real*8 :: e
      real*8 :: f
      real*8 :: feints
      real*8 :: g(3)
      integer i
      integer n
      real*8 :: pi(3)
      real*8 :: result
      real*8 :: sum1
      real*8 :: x(num)
      real*8 :: x1
      real*8 :: x2
      real*8 :: x3
      real*8 :: y(num)
!     
      result = 0.0d0
      
      if ( num <= 2 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SIMPNE - Fatal error!'
         write ( *, '(a)' ) '  NUM <= 2.'
         stop
      end if
      
      n = 1
      
      do
         
         x1 = x(n)
         x2 = x(n+1)
         x3 = x(n+2)
         e = x3*x3-x1*x1
         f = x3*x3*x3-x1*x1*x1
         feints = x3-x1
         del(1) = x3-x2
         del(2) = x1-x3
         del(3) = x2-x1
         g(1) = x2+x3
         g(2) = x1+x3
         g(3) = x1+x2
         pi(1) = x2*x3
         pi(2) = x1*x3
         pi(3) = x1*x2
         
         sum1 = 0.0d0
         do i = 1, 3
            sum1 = sum1 + y(n-1+i)*del(i)*(f/3.0d0-g(i)*0.5d0*e+pi(i)
     &           *feints)
         end do
         result = result - sum1 / ( del(1) * del(2) * del(3) )
         
         n = n+2
         
         if ( n + 1 >= num ) then
            exit
         end if
         
      end do
      
      if ( mod(num,2) /= 0 ) then
         return
      end if
      
      n = num-2
      x3 = x(num)
      x2 = x(num-1)
      x1 = x(num-2)
      e = x3*x3-x2*x2
      f = x3*x3*x3-x2*x2*x2
      feints = x3-x2
      del(1) = x3-x2
      del(2) = x1-x3
      del(3) = x2-x1
      g(1) = x2+x3
      g(2) = x1+x3
      g(3) = x1+x2
      pi(1) = x2*x3
      pi(2) = x1*x3
      pi(3) = x1*x2
      
      sum1 = 0.0d0
      do i = 1, 3
         sum1 = sum1 + y(n-1+i) * del(i) * 
     &     ( f / 3.0d0 - g(i) * 0.5d0 * e + pi(i) * feints )
      end do
      
      result = result - sum1 / ( del(1) * del(2) * del(3) )
      
      return
      end
      
