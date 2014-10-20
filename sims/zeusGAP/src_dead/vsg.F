      subroutine vsgiso(cylindrad, zaxis, beta)
      use planet
c     we solve for the approximate vertical structure of a self-gravitating, vertically isothermal disk    
      implicit none  
      real*8, intent(in) :: cylindrad, zaxis
      real*8, intent(out) :: beta
      real*8 :: beta0, bH, paramK, zprime, radius 
      real*8, external :: bigK, get_varray, bigH
      
      common /share/ paramK, radius 
      common /sharebeta/ beta0 
      radius = cylindrad 
c     Initialize a guess for beta0
      beta0 = beta_trial
c     get the K parameter for this radius. and re-scale the vertical axis (make it z')
      
      bH = bigH(cylindrad)
      paramK    = bigK(cylindrad)
      zprime = zaxis/(bH*sqrt(paramK))
     
      call get_beta0(beta0) ! Solve for beta_0
      
      beta = get_varray(zprime)*exp(0.5*(zaxis/bH)**2)
           
c     print*, 'got here', beta
      return
      end
      real*8 function bigK(r)
c     Get the radius-dependent paramter K in the ODE
      use planet
      implicit none
      real*8, parameter :: bigG = 1d0, pi=2d0*acos(0d0)
      real*8 :: r, ciso, omega_k
      real*8, external :: surfden, csq
      ciso = sqrt(csq(r))
      omega_k = r**(-3d0/2d0)
      bigK = ciso*omega_k/(pi*bigG*surfden(r))
      bigK = bigK*sqrt(pi/8d0)
      return
      end function bigK
      real*8 function surfden(r)
!the density is sigma_0 * f(r) *exp(-z^2/2H^2) /( sqrt(2*pi)*H)
!this function calculations the sigma_0 * f(r) \equiv sigma(r) 
      use planet
      use domain
      implicit none
      real*8 :: r, hin
      real*8, external :: bigH, inner_hole
   
      surfden = den0*r**(-sig)*inner_hole(r)
      
      return 
      end function surfden
      subroutine get_beta0(beta_0)
!! beta_0 is an initial condition, but we do not know how to set its value directly.
!! instead, we impose the surface density after SG-modification remains the same as the 
!! standard gaussian profile. beta_0 will be updated to the "correct" value by this routine 
      
      implicit none
      logical :: check
      integer, parameter :: n = 1
      real*8 :: beta_0, paramK, radius
      common /share/ paramK, radius 
      
      call newt(beta_0, n, check)
      
      return
      end subroutine get_beta0
      subroutine funcv(n, x, fvec)
!     return the function to be zeroed by adjusting the beta_0 parameter 
      use planet 
      use domain 
      implicit none
 
      integer :: n
      real*8, parameter :: eps=1d-9, pi=2d0*acos(0d0)
      real*8 :: x, fvec(n), paramK, a, b,
     &     res, beta0, n_H, radius 
      external get_varray
      real*8, external :: bigH 
      common /share/ paramK, radius 
      common /sharebeta/ beta0
      n_H = radius*tan(pi/2d0 - thetamin)/bigH(radius)
      a = 0d0
      b = n_H/sqrt(paramK)
      beta0 = x 
      call simp (get_varray, a, b, eps, res)
      
      fvec(1) = erf(n_H/sqrt(2d0)) - sqrt(2d0*paramK/pi)*res 
      end subroutine funcv
      real*8 function get_varray(ztarget)
!     given a beta(z=0), solve 1st order ODE to get v(z') at some desired 
      use planet 
      implicit none
      integer, parameter :: nvar = 2
      real*8, parameter :: eps=1d-9
      integer :: nok, nbad 
      real*8, intent(in) :: ztarget  
      real*8 :: paramK, y(nvar), t, tout, h1, hmin, beta0, radius 
      external dvdzprime, rkqs1 
      common /share/ paramK, radius
      common /sharebeta/ beta0 
      
      y(1) = dlog(beta0) ;  y(2) = 0d0 
      t = 0d0 ; tout = ztarget 
      
!     ! numerical recipes RK45 solver.
      h1 = zmax/1024d0 ; hmin = 0d0
      call odeint1(y, nvar, t, tout, eps, h1, hmin, nok, nbad, 
     &     dvdzprime, rkqs1)
      
      get_varray = exp(y(1))    ! Want to integrate this w.r.t. z' 
      
      return
      end function get_varray
      subroutine dvdzprime(t, y, yp)
!     evaluate derivative
      implicit none
      integer, parameter :: neq = 2
      real*8 :: t, y(neq), yp(neq)
      real*8 :: paramK, radius 
      common /share/ paramK, radius
      
!     ! two ODE form. this can be used to integrate from exactly z=0. y1 = v, y2 = dv/dz
      yp(1)  = y(2)
      yp(2)  = -paramK - exp(y(1))
      
      return
      end subroutine dvdzprime
      
      subroutine simp ( func, a, b, eps, result )
!*****************************************************************************80
!
!! SIMP approximates the integral of a function by an adaptive Simpson's rule.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    James Lyness,
!    Algorithm 379:
!    SQUANK - Simpson quadrature used adaptively, noise killed,
!    Communications of the ACM,
!    Volume 13, Number 4, April 1970, pages 260-263.
!
!    William McKeeman, Lawrence Tesler,
!    Algorithm 182:
!    Nonrecursive adaptive integration,
!    Communications of the ACM,
!    Volume 6, 1963, page 315.
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      FUNCTION FUNC ( X )
!    which evaluates the function at the point X.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input, real ( kind = 8 ) EPS, the requested error tolerance.
!
!    Output, real ( kind = 8 ) RESULT, the approximation to the integral.
!
      implicit none
      integer ( kind = 4 ), parameter :: maxlev = 30
      
      real ( kind = 8 ) a
      real ( kind = 8 ) a1
      real ( kind = 8 ) absar
      real ( kind = 8 ) b
      real ( kind = 8 ) da
      real ( kind = 8 ) dx(maxlev)
      real ( kind = 8 ) ep
      real ( kind = 8 ) eps
      real ( kind = 8 ) epsp(maxlev)
      real ( kind = 8 ) est
      real ( kind = 8 ) est1
      real ( kind = 8 ) est2(maxlev)
      real ( kind = 8 ) est3(maxlev)
      real ( kind = 8 ) f1
      real ( kind = 8 ) f2(maxlev)
      real ( kind = 8 ) f3(maxlev)
      real ( kind = 8 ) f4(maxlev)
      real ( kind = 8 ) fa
      real ( kind = 8 ) fb
      real ( kind = 8 ) fbp(maxlev)
      real ( kind = 8 ) fm
      real ( kind = 8 ) fmp(maxlev)
      real ( kind = 8 ), external :: func
      integer ( kind = 4 ) l
      integer ( kind = 4 ) lvl
      integer ( kind = 4 ) nrtr(maxlev)
      real ( kind = 8 ) pval(maxlev,3)
      real ( kind = 8 ) result
      real ( kind = 8 ) sum1
      real ( kind = 8 ) sx
      real ( kind = 8 ) x2(maxlev)
      real ( kind = 8 ) x3(maxlev)
      
      result = 0.0D+00
      
      if ( a == b ) then
         return
      end if
      
      ep = eps
      a1 = a
      nrtr(1:maxlev) = 0
      pval(1:maxlev,1:3) = 0.0D+00
      
      lvl = 0
      absar = 0.0D+00
      est = 0.0D+00
      da = b - a1
      
      fa = func ( a1 )
      fm = 4.0D+00 * func ( ( a1 + b ) * 0.5D+00 )
      fb = func ( b )
!     
!     1 = RECUR
!     
 30   continue
      
      lvl = lvl + 1
      dx(lvl) = da / 3.0D+00
      sx = dx(lvl) / 6.0D+00
      f1 = 4.0D+00 * func ( 0.5D+00 * dx(lvl) + a1 )
      x2(lvl) = a1 + dx(lvl)
      f2(lvl) = func ( x2(lvl) )
      x3(lvl) = x2(lvl) + dx(lvl)
      f3(lvl) = func ( x3(lvl) )
      epsp(lvl) = ep
      f4(lvl) = 4.0D+00 * func ( dx(lvl) * 0.5D+00 + x3(lvl) )
      fmp(lvl) = fm
      est1 = sx * (fa+f1+f2(lvl))
      fbp(lvl) = fb
      est2(lvl) = sx * ( f2(lvl) + f3(lvl) + fm )
      est3(lvl) = sx * ( f3(lvl) + f4(lvl) + fb )
      sum1 = est1 + est2(lvl) + est3(lvl)
      absar = absar - abs ( est ) + abs ( est1 ) + abs ( est2(lvl) ) 
     &     + abs ( est3(lvl) )
      
      if ( abs ( est - sum1 ) <= epsp(lvl) * absar ) then
         go to 40
      end if
      
      if ( maxlev <= lvl ) then
         go to 50
      end if
!     
!     2 = UP
!     
 40   continue
      
      if ( 1 < lvl ) then
         lvl = lvl-1
      end if
      
      l = nrtr(lvl)
      
      if ( l == 0 ) then
         go to 50
      end if
      
      pval(lvl,l) = sum1
      
      if ( l == 1 ) then
         go to 60
      end if
      
      if ( l == 2 ) then
         go to 70
      end if
      
      if ( l == 3 ) then
         go to 80
      end if
      
 50   continue
      
      nrtr(lvl) = 1
      est = est1
      fm = f1
      fb = f2(lvl)
      ep = epsp(lvl) / 1.7D+00
      da = dx(lvl)
      go to 30
      
 60   continue
      
      nrtr(lvl) = 2
      fa = f2(lvl)
      fm = fmp(lvl)
      fb = f3(lvl)
      est = est2(lvl)
      a1 = x2(lvl)
      ep = epsp(lvl) / 1.7D+00
      da = dx(lvl)
      go to 30
      
 70   continue
      
      nrtr(lvl) = 3
      fa = f3(lvl)
      fm = f4(lvl)
      fb = fbp(lvl)
      est = est3(lvl)
      a1 = x3(lvl)
      ep = epsp(lvl) / 1.7D+00
      da = dx(lvl)
      go to 30
      
 80   continue
      
      sum1 = pval(lvl,1) + pval(lvl,2) + pval(lvl,3)
      
      if ( 1 < lvl ) then
         go to 40
      end if
      
 90   continue
      
      result = sum1
      
      return
      end subroutine simp
      
      SUBROUTINE odeint1(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,
     &  derivs,rkqs1)
c
c .................................................................
c .   Subroutine odeint 
c .   Runge-Kutta driver with adaptive stepsize control.  Integrate
c .   the starting values ystart(1:nvar) from x1 to x2 with accuracy 
c .   eps, storing intermediate results in the common block /path/.
c .   h1 should be set as a guessed first stepsize, hmin as the 
c .   minimum allowed stepsize (can be zero).
c .   On output, nok and nbad are the number of good and bad (but
c .   retried and fixed) steps taken, and ystart is replaced by 
c .   values at the end of the integration interval.
c .   derivs is the user-supplied subroutine for calculating the
c .   right-hand side derivatives, while rkqs is the name of 
c .   the stepper routine to be used.  
c .   /path/ contains its own information about how often an
c .   intermediate value is to be stored.
c .................................................................
c
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs1
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.0d-30)
      INTEGER i,kmax,kount,nstp
      REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     & yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do i=1,nvar
        y(i)=ystart(i)
      end do
      if (kmax.gt.0) xsav=x-2.*dxsav
      do nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
        end do
        if(kmax.gt.0) then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1) then
              kount=kount+1
              xp(kount)=x
              do i=1,nvar
                yp(i,kount)=y(i)
              end do
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1) .gt. 0.0d0) h=x2-x
        call rkqs1(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h) then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1) .ge. 0.0d0)then
          do i=1,nvar
            ystart(i)=y(i)
          end do
          if(kmax.ne.0) then
            kount=kount+1
            xp(kount)=x
            do i=1,nvar
              yp(i,kount)=y(i)
            end do
          endif
          return
        endif
        if(abs(hnext) .lt. hmin) pause
     &     'stepsize smaller than minimum in odeint'
        h=hnext
        end do
      pause 'too many steps in odeint'
      return
      END
      SUBROUTINE rkqs1(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
c
c ................................................................
c .   Subroutine rkqs is a fifth-order Runge-Kutta step with 
c .   monitoring of local truncation error to ensure accuracy
c .   and adjust stepsize.   
c .
c .   Inputs:
c .   y(1:n) = the dependent variable vector
c .   dydx(1:n) = deriv of y wrt x at starting value of indep var x
c .   htry = attempted step size
c .   eps = required accuracy
c .   yscal(1:n) = scaling vector for the error
c .   
c .   Outputs:
c .   y = new value of y
c .   x = new value of x
c .   hdid = actual stepsize accomplished
c .   hnext = estimated next stepsize   
c .
c .   derivs is the user-supplied subroutine that computes the 
c .    right-hand side derivatives
c ................................................................
c
      INTEGER n,NMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
c     Maximum number of equations
      INTEGER i
      REAL*8 errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     &  PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry
 1    call rkck1(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.0d0
      do i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
      end do
      errmax=errmax/eps
      if(errmax.gt.1.0d0) then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp), 0.1*abs(h)), h)
        xnew=x+h
        if(xnew.eq.x) pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON) then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.0d0*h
        endif
        hdid=h
        x=x+h
        do i=1,n
          y(i)=ytemp(i)
        end do
        return
      endif
      END
      SUBROUTINE rkck1(y,dydx,n,x,h,yout,yerr,derivs)
c
c ...............................................................
c .   Subroutine rkck.
c .   Given values for n variables y and their derivatives dydx
c .   known at x, use the fifth-order Cash-Karp Runge-Kutta method
c .   to advance the solution over an interval h and return the 
c .   incremented variables as yout.  Also return an estimate of
c .   the local truncation error in yout using the embedded 
c .   fourth-order method.  The user supplies the subroutine
c .   derivs (x,y,dydx) which returns derivatives dydx at x.
c ................................................................
c
      INTEGER n,NMAX
      REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
C     USES derivs
      INTEGER i
      REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     &  ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,
     &  B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,
     &  DC5,DC6
      PARAMETER (A2=0.2,A3=0.3,A4=0.6,A5=1.0,A6=0.875,B21=0.2,
     &  B31=3.0/40.0,B32=9.0/40.0,B41=0.3,B42=-0.9,B43=1.2,
     &  B51=-11.0/54.0,B52=2.5,B53=-70.0/27.0,B54=35.0/27.0,
     &  B61=1631.0/55296.0,B62=175.0/512.0,B63=575.0/13824.0,
     &  B64=44275.0/110592.0,B65=253.0/4096.0,C1=37.0/378.0,
     &  C3=250.0/621.0,C4=125.0/594.0,C6=512.0/1771.0,
     &  DC1=C1-2825.0/27648.0,DC3=C3-18575.0/48384.0,
     &  DC4=C4-13525.0/55296.0,DC5=-277.0/14336.0,DC6=C6-0.25)
      do i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
      end do 
      call derivs(x+A2*h,ytemp,ak2)
      do i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      end do
      call derivs(x+A3*h,ytemp,ak3)
      do i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      end do
      call derivs(x+A4*h,ytemp,ak4)
      do i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
      end do
      call derivs(x+A5*h,ytemp,ak5)
      do i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     &     B65*ak5(i))
      end do
      call derivs(x+A6*h,ytemp,ak6)
      do i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
      end do
c     estimate error as difference between fourth and fifth order methods
      do i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+
     &   DC6*ak6(i))
      end do
      return
      END
c     NEWTON SCHEME
c
      SUBROUTINE newt(x,n,check)
      INTEGER n,nn,NP,MAXITS
      LOGICAL check
      REAL*8 x(n),fvec,TOLF,TOLMIN,TOLX,STPMX
      PARAMETER (NP=40,MAXITS=200,TOLF=1.d-4,TOLMIN=1.d-6,TOLX=1.d-7,
     *STPMX=100.d0)
      COMMON /newtv/ fvec(NP),nn
      SAVE /newtv/
CU    USES fdjac,fmin,lnsrch,lubksb,ludcmp
      INTEGER i,its,j,indx(NP)
      REAL*8 d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP),
     *xold(NP),fmin
      EXTERNAL fmin
      nn=n
      f=fmin(x)
      test=0.
      do 11 i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
11    continue
      if(test.lt..01*TOLF)then
        check=.false.
        return
      endif
      sum=0.
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum),float(n))
      do 21 its=1,MAXITS
        call fdjac(n,x,fvec,NP,fjac)
        do 14 i=1,n
          sum=0.
          do 13 j=1,n
            sum=sum+fjac(j,i)*fvec(j)
13        continue
          g(i)=sum
14      continue
        do 15 i=1,n
          xold(i)=x(i)
15      continue
        fold=f
        do 16 i=1,n
          p(i)=-fvec(i)
16      continue
        call ludcmp(fjac,n,NP,indx,d)
        call lubksb(fjac,n,NP,indx,p)
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
        test=0.
        do 17 i=1,n
          if(abs(fvec(i)).gt.test)test=abs(fvec(i))
17      continue
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then
          test=0.
          den=max(f,.5*n)
          do 18 i=1,n
            temp=abs(g(i))*max(abs(x(i)),1.)/den
            if(temp.gt.test)test=temp
18        continue
          if(test.lt.TOLMIN)then
            check=.true.
          else
            check=.false.
          endif
          return
        endif
        test=0.
        do 19 i=1,n
          temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
          if(temp.gt.test)test=temp
19      continue
        if(test.lt.TOLX)return
21    continue
      pause 'MAXITS exceeded in newt'
      END
      SUBROUTINE fdjac(n,x,fvec,np,df)
      INTEGER n,np,NMAX
      REAL*8 df(np,np),fvec(n),x(n),EPS
      PARAMETER (NMAX=40,EPS=1.e-4)
CU    USES funcv
      INTEGER i,j
      REAL*8 h,temp,f(NMAX)
      do 12 j=1,n
        temp=x(j)
        h=EPS*abs(temp)
        if(h.eq.0.)h=EPS
        x(j)=temp+h
        h=x(j)-temp
        call funcv(n,x,f)
        x(j)=temp
        do 11 i=1,n
          df(i,j)=(f(i)-fvec(i))/h
11      continue
12    continue
      return
      END
      FUNCTION fmin(x)
      INTEGER n,NP
      REAL*8 fmin,x(*),fvec
      PARAMETER (NP=40)
      COMMON /newtv/ fvec(NP),n
      SAVE /newtv/
CU    USES funcv
      INTEGER i
      REAL*8 sum
      call funcv(n,x,fvec)
      sum=0.
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin=0.5*sum
      return
      END
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
      INTEGER n
      LOGICAL check
      REAL*8 f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
      PARAMETER (ALF=1.e-4,TOLX=1.e-7)
      EXTERNAL func
CU    USES func
      INTEGER i
      REAL*8 a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
     *test,tmplam
      check=.false.
      sum=0.
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=func(x)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.)then
            tmplam=-slope/(2.*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.)then
              tmplam=-slope/(2.*b)
            else
              disc=b*b-3.*a*slope
              if(disc.lt.0.) pause 'roundoff problem in lnsrch'
              tmplam=(-b+sqrt(disc))/(3.*a)
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1*alam)
      goto 1
      END
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
