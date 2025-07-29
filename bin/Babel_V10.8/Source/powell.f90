module powell_min

  implicit none

  public :: powell
  private :: PO_linmin, PO_f1dim, PO_mnbrak, PO_brent

  REAL(kind(0.d0)), DIMENSION(:), allocatable, save, private ::pcom,xicom
  INTEGER, save, private :: ncom


  !----------------------
contains

  ! taken from Numerical Recipes in Fortran 77 and adapted for Fortran 90
  ! http://www.ulib.org/webRoot/Books/Numerical_Recipes/


  SUBROUTINE powell(fcost,p,xi,n,np,ftol,iter,fret, print_fCost)

    implicit none

    REAL(kind(0.d0)), external :: fcost
    INTEGER, intent(in) :: np   ! Maximal size of the vector p
    INTEGER, intent(in) :: n    ! Real size of the vector p
    REAL(kind(0.d0)), DIMENSION(np), intent(inout) :: p
    REAL(kind(0.d0)), DIMENSION(np,np), intent(inout) :: xi
    REAL(kind(0.d0)), intent(in) :: ftol
    INTEGER, intent(out) :: iter
    REAL(kind(0.d0)), intent(out) :: fret
    EXTERNAL :: print_fCost

    INTEGER, PARAMETER :: ITMAX=200

    ! USES fcost,PO_linmin
    
    !   Minimization of a function fcost of n variables.  Input consists of an 
    !   initial starting point p(1:n); an initial matrix xi(1:n,1:n) with 
    !   physical dimensions np by np, and whose columns contain the initial 
    !   set of directions (usually the n unit vectors); and ftol, the 
    !   fractional tolerance in the function value such that failure to 
    !   decrease by more than this amount on one iteration signals doneness. 
    !   On output, p is set to the best point found, xi is the then-current 
    !   direction set, fret is the returned function value at p, and iter 
    !   is the number of iterations taken. The routine PO_linmin is used. 

    ! PARAMETERS: Maximum value of n, maximum allowed iterations, and a 
    !   small number.

    INTEGER :: i,ibig,j
    REAL(kind(0.d0)) :: del,fp,fptt,t
    REAL(kind(0.d0)), DIMENSION(np) :: pt,ptt,xit
  
    ! Table allocation
    ALLOCATE(pcom(1:np)) ; ALLOCATE(xicom(1:np))
    ncom=n

    fret=fcost(p,n)
    pt(1:n)=p(1:n)  ! save the initial point
    iter=0

1   iter=iter+1

    !--------------------------
    ! Print current values 
    write(6,'(a,i0,a,g20.12)') 'iter = ', iter, ' - fcost= ', fret
    DO i=1, n
       WRITE(6,'(a,i0,a,g20.12)') '  x(',i,') = ', p(i)
    END DO
    CALL Print_fCost(6)
    !--------------------------
    fp=fret
    ibig=0
    del=0.d0        ! will be the biggest function decrease
    do i=1,n        ! in each iteration, loop over all directions in the set
       xit(1:n)=xi(1:n,i)            ! copy te direction
       fptt=fret
       call PO_linmin(fcost, p, xit, n, np, fret)     ! minimize along it
       if(abs(fptt-fret).gt.del)then ! and record it if it is the largest 
          del=abs(fptt-fret)         !   decrease so far
          ibig=i
       endif
    enddo
    if(2.d0*abs(fp-fret).le.ftol*(abs(fp)+abs(fret))) then
            DeAllocate(pcom) ; DeAllocate(xicom)
            return!Termination criter.
    END IF
    if(iter.eq.ITMAX) STOP 'powell exceeding maximum iterations'
    do j=1,n                         ! construct the extrapolated point and the
       ptt(j)=2.d0*p(j)-pt(j)        !   average direction moved.
       xit(j)=p(j)-pt(j)             !   save tge old starting point  
       pt(j)=p(j)
    enddo
    fptt=fcost(ptt,n)                   ! function value at extrapolated point
    if(fptt.ge.fp)goto 1             ! one reason not to use new direction  
    t=2.d0*(fp-2.d0*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
    if(t.ge.0.d0)goto 1              ! other reason 
    call PO_linmin(fcost, p, xit, n, np, fret)     ! move to the minimum of the new direction
    do j=1,n                         ! and save the new direction
       xi(j,ibig)=xi(j,n)
       xi(j,n)=xit(j)
    enddo
    goto 1        ! back for another iteration          
    
  END SUBROUTINE powell
  !--------------------------------------------------------------------------
  SUBROUTINE PO_linmin(fcost, p, xi, n, np, fret)
    
    implicit none
    
    REAL(kind(0.d0)), external :: fcost
    INTEGER, intent(in) :: n, np
    REAL(kind(0.d0)), DIMENSION(np), intent(inout) :: p
    REAL(kind(0.d0)), DIMENSION(np), intent(inout) :: xi
    REAL(kind(0.d0)), intent(out) :: fret

    REAL(kind(0.d0)), PARAMETER :: TOL=1.d-4
    
    ! USES PO_brent,PO_f1dim,PO_mnbrak
    
    !   Given an n-dimensional point p(1:n) and an n-dimensional direction 
    !   xi(1:n), moves and resets p to where the function fcost(p,n) takes on 
    !   a minimum along the direction xi from p, and replaces xi by the 
    !   actual vector displacement that p was moved. Also returns as fret 
    !   the value of fcost at the returned location p. This is actually all 
    !   accomplished by calling the routines PO_mnbrak and PO_brent. 
    
    INTEGER :: j
    REAL(kind(0.d0)) :: ax,bx,fa,fb,fx,xmin,xx
    
    pcom(1:n)=p(1:n)
    xicom(1:n)=xi(1:n)
    
    ax=0.d0               ! initial guess for brackets
    xx=1.d0
    
    call PO_mnbrak(fcost,n,ax,xx,bx,fa,fx,fb)
    fret=PO_brent(fcost,n,ax,xx,bx,TOL,xmin)
    
    do j=1,n              ! construct the vector result to return
       xi(j)=xmin*xi(j)
       p(j)=p(j)+xi(j)
    enddo
    return
    
  END SUBROUTINE PO_linmin
  
  !----------------------------------------------------------------------------
  
  FUNCTION PO_f1dim(fcost, x, np)

    implicit none
    
    REAL(kind(0.d0)), external :: fcost
    REAL(kind(0.d0)), intent(in) :: x
    INTEGER, intent(in) :: np
    REAL(kind(0.d0)) :: PO_f1dim
    
    !   Used by PO_linmin as the function passed to PO_mnbrak and PO_brent.
    
    REAL(kind(0.d0)), DIMENSION(np) :: xt
    
    xt(1:ncom)=pcom(1:ncom)+x*xicom(1:ncom)
    PO_f1dim=fcost(xt, np)
    return
    
  END FUNCTION PO_f1dim
  
  !---------------------------------------------------------------------------
  FUNCTION PO_brent(fcost, np, ax,bx,cx,tol,xmin)
    
    implicit none
    
    REAL(kind(0.d0)), external :: fcost
    INTEGER, intent(in) :: np
    REAL(kind(0.d0)), intent(in) ::  ax,bx,cx,tol
    REAL(kind(0.d0)), intent(out) ::  xmin
    REAL(kind(0.d0)) ::  PO_brent

    INTEGER, PARAMETER :: ITMAX=100
    REAL(kind(0.d0)), PARAMETER :: CGOLD=.3819660d0,ZEPS=1.0d-10
    
    !    Given a function f, and given a bracketing triplet of abscissas 
    !    ax, bx, cx (such that bx is between ax and cx, and f(bx) is less 
    !    than both f(ax) and f(cx)), this routine isolates the minimum to a 
    !    fractional precision of about tol using Brent's method. The abscissa 
    !    of the minimum is returned as xmin, and the minimum function value 
    !    is returned as PO_brent, the returned function value. 
    
    ! Parameters: Maximum allowed number of iterations; 
    !    golden ratio; and a small number that protects against trying 
    !    to achieve fractional accuracy for a minimum that happens to be 
    !    exactly zero.
    
    INTEGER :: iter
    REAL(kind(0.d0)) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.d0
    fx=PO_f1dim(fcost, x, np)
    fv=fx
    fw=fx
    do iter=1,ITMAX
       xm=0.5d0*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.d0*tol1
       if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
       if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5d0*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x) 
          goto 2
       endif
1      if(x.ge.xm) then
          e=a-x
       else
          e=b-x
       endif
       d=CGOLD*e
2      if(abs(d).ge.tol1) then
          u=x+d
       else
          u=x+sign(tol1,d)
       endif
       fu=PO_f1dim(fcost, u, np)
       if(fu.le.fx) then
          if(u.ge.x) then
             a=x
          else
             b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
       else
          if(u.lt.x) then
             a=u
          else
             b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
             v=w
             fv=fw
             w=u
             fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
             v=u
             fv=fu
          endif
       endif
    end do
    STOP 'PO_brent exceed maximum iterations'
3   xmin=x
    PO_brent=fx
    return
    
  END FUNCTION PO_brent
  
  !---------------------------------------------------------------------------
  
  SUBROUTINE PO_mnbrak(fcost, np, ax,bx,cx,fa,fb,fc)
    
    implicit none
    
    REAL(kind(0.d0)), external :: fcost
    INTEGER, intent(in) :: np
    REAL(kind(0.d0)), intent(inout) :: ax,bx
    REAL(kind(0.d0)), intent(out) :: cx,fa,fb,fc

    REAL(kind(0.d0)), PARAMETER :: GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20

    !  Given a function func, and given distinct initial points ax and bx, 
    !  this routine searches in the downhill direction (defined by the 
    !  function as evaluated at the initial points) and returns new points 
    !  ax, bx, cx that bracket a minimum of the function. Also returned are 
    !  the function values at the three points, fa, fb, andfc. 
    
    !  Parameters: GOLD is the default ratio by which successive intervals 
    !    are magnified; 
    !    GLIMIT is the maximum magnification allowed for a parabolic-fit step.
    
    REAL(kind(0.d0)) :: dum,fu,q,r,u,ulim
    
    fa=PO_f1dim(fcost, ax, np)
    fb=PO_f1dim(fcost, bx, np)
    if(fb.gt.fa)then
       dum=ax
       ax=bx
       bx=dum
       dum=fb
       fb=fa
       fa=dum
    endif
    cx=bx+GOLD*(bx-ax)
    fc=PO_f1dim(fcost, cx, np)
1   if(fb.ge.fc)then
       r=(bx-ax)*(fb-fc)
       q=(bx-cx)*(fb-fa)
       u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
       ulim=bx+GLIMIT*(cx-bx)
       if((bx-u)*(u-cx).gt.0.d0)then
          fu=PO_f1dim(fcost, u, np)
          if(fu.lt.fc)then
             ax=bx
             fa=fb
             bx=u
             fb=fu
             return
          else if(fu.gt.fb)then
             cx=u
             fc=fu
             return
          endif
          u=cx+GOLD*(cx-bx)
          fu=PO_f1dim(fcost, u, np)
       else if((cx-u)*(u-ulim).gt.0.d0)then
          fu=PO_f1dim(fcost, u, np)
          if(fu.lt.fc)then
             bx=cx
             cx=u
             u=cx+GOLD*(cx-bx)
             fb=fc
             fc=fu
             fu=PO_f1dim(fcost, u, np)
          endif
       else if((u-ulim)*(ulim-cx).ge.0.d0)then
          u=ulim
          fu=PO_f1dim(fcost, u, np)
       else
          u=cx+GOLD*(cx-bx)
          fu=PO_f1dim(fcost, u, np)
       endif
       ax=bx
       bx=cx
       cx=u
       fa=fb
       fb=fc
       fc=fu
       goto 1
    endif
    return
  END SUBROUTINE PO_mnbrak

end module powell_min
