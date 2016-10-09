program main
use module_input
implicit none 
real,dimension(imin:imax)     :: a,dx,r,z,a_true,r_true
real,dimension(imin-1:imax+1) :: x			 
real,dimension(jmin-1:jmax+1) :: aR,aL,a6,dela,abar,fL,fR
real,dimension(jmin-2:jmax+1) :: a_12
real :: u=1. ,dt  
integer :: i,j



x  = (/((i*length/(SIZE(x)-1)),i=0,SIZE(x)-1,1)/)   ! Make domain uniformly spaced
dx =  x(imin:imax)-x(imin-1:imax-1)   				! Spacing of cells
r  = .5*(x(imin:imax)+x(imin-1:imax-1))         	! Cell centers
dt = Courant*MINVAL(dx)/abs(u)

call OpenFiles()
call InitCond(a,r)
call Initcond(a_true,r)
write(2,*) r(jmin:jmax)

do i=1,itterMax
	write(1,*) a(jmin:jmax)
    write(3,*) a_true(jmin:jmax)
    
    if (Method.EQ.'PPM') then
        call PPM_interpol(a,dx,aL,aR,dela,a6,a_12)
        z = (u*dt/dx)
        fL(jmin-1:jmax) = aR(jmin-1:jmax) - .5*z(jmin-1:jmax)*(dela(jmin-1:jmax)-(1.-2.*z(jmin-1:jmax)/3.)*a6(jmin-1:jmax))
        z = (-u*dt/dx)
        fR(jmin-1:jmax) = aL(jmin:jmax+1) + .5*z(jmin:jmax+1)*(dela(jmin:jmax+1)+(1.-2.*z(jmin:jmax+1)/3.)*a6(jmin:jmax+1))
    else if (Method.EQ.'MUSCL') then
        call MUSCL_interpol(a,dx,aL,aR,dela)
        z = (u*dt)
        fL(jmin-1:jmax) = a(jmin-1:jmax)+dela(jmin-1:jmax)*(x(jmin-1:jmax)-r(jmin-1:jmax)-.5*z(jmin-1:jmax))/dx(jmin-1:jmax)
        z = (-u*dt)
        fR(jmin-1:jmax) = a(jmin:jmax+1)+dela(jmin:jmax+1)*(x(jmin-1:jmax)-r(jmin:jmax+1)+.5*z(jmin-1:jmax))/dx(jmin:jmax+1)
    else if (Method.EQ.'Godunov') then
        fL(jmin-1:jmax) = a(jmin-1:jmax)
        fR(jmin-1:jmax) = a(jmin:jmax+1)
    end if 
	if (u.LT.delta) then
		abar(jmin-1:jmax) = fR(jmin-1:jmax)
	else 
		abar(jmin-1:jmax) = fL(jmin-1:jmax)
	end if 
	a(jmin:jmax) = a(jmin:jmax) + (u*dt/dx(jmin:jmax))*(abar(jmin-1:jmax-1)-abar(jmin:jmax))
    r_true = (/(MODULO(r(j)-u*i*dt,r(jmax)-r(jmin)),j=imin,imax,1)/)  ! Makes true solution periodic
    call InitCond(a_true,r_true)
    call Boundaries(a,u)
end do 

end program

subroutine Boundaries(a,u)
use module_input
implicit none
real,dimension(imin:imax),intent(inout) :: a
real,intent(in) :: u
integer :: i
if (u.GT.delta) a(jmin) = a(jmax)
if (u.LE.delta) a(jmax) = a(jmin)
end subroutine

subroutine PPM_interpol(a,dx,aL,aR,dela,a6,a_12)
use module_input
implicit none
real,intent(in),dimension(imin:imax) :: a,dx
real,intent(inout),dimension(jmin-1:jmax+1) :: aR,aL,a6,dela
real,intent(inout),dimension(jmin-2:jmax+1) :: a_12
real,dimension(jmin-2:jmax+1) :: one,two,three,four,five  ! Cells in which we need the a_j+1/2
real,dimension(jmin-2:jmax+1) :: del
integer :: js,je
js = jmin-2   ! index of first cell interpolated, that is a_j+1/2
je = jmax+1   ! inded of last  cell interpolated, that is a_j+1/2
one(js:je)   = dx(js:je)/(dx(js:je)+dx(js+1:je+1))
two(js:je)   = 1./(dx(js-1:je-1)+dx(js:je)+dx(js+1:je+1)+dx(js+2:je+2))
three(js:je) = (2.*dx(js+1:je+1)*dx(js:je))/(dx(js:je)+dx(js+1:je+1))
four(js:je)  = (dx(js-1:je-1)+dx(js:je))/(2.*dx(js:je)+dx(js+1:je+1))
five(js:je)  = (dx(js+2:je+2)+dx(js+1:je+1))/(2.*dx(js+1:je+1)+dx(js:je))

CALL PPM_avgSlope(a,dx,del)
a_12(js:je) = a(js:je) + one(js:je)*(a(js+1:je+1)-a(js:je)) &
            + two(js:je)*(three(js:je)*(four(js:je)-five(js:je))*(a(js+1:je+1)-a(js:je)) &
            - dx(js:je)*four(js:je)*del(js+1:je+1)+dx(js+1:je+1)*five(js:je)*del(js:je))  
aL(jmin-1:jmax+1) = a_12(jmin-2:jmax)   
aR(jmin-1:jmax+1) = a_12(jmin-1:jmax+1)   
call PPM_monotone(a,aL,aR)
dela = aR-aL
a6 = 6.*(a(jmin-1:jmax+1)-.5*(aL+aR))
end subroutine

subroutine PPM_monotone(a,aR,aL)
!-------------------------------
! Ensures interpolation is monotone using CW Eq. 1.10
! Input:
!   a  - Current solution
! Input/Output:
!   aR - Approx. of a_(j+1/2) from cell to right
!   aL - Approx. of a_(j-1/2) from cell to left
!-----------------------------
use module_input
implicit none
real,intent(inout),dimension(jmin-1:jmax+1) :: aR,aL
real,intent(in),dimension(imin:imax)       :: a
real,dimension(jmin-1:jmax+1) ::cond1,cond2,cond3
integer :: i
cond1 = (aR-a(jmin-1:jmax+1))*(a(jmin-1:jmax+1)-aL)
cond2 = (aR-aL)*(a(jmin-1:jmax+1)-.5*(aL+aR))
cond3 = ((aR-aL)**(2))/6.
do i=jmin-1,jmax+1  
    if (cond2(i).GT.cond3(i)) then
        aL(i) = 3.*a(i)-2.*aR(i)
    end if 
    if (-cond3(i).GT.cond2(i)) then
        aR(i) = 3.*a(i)-2.*aL(i)
    end if 
    if (cond1(i).LE.0.) then
        aL(i) = a(i)
        aR(i) = a(i)
    end if 
end do
end subroutine

subroutine PPM_avgSlope(a,dx,delmj)
!------------------------------
! Calculate the average slope in the jth cell of
! the parabola with zone average a j,j+1,j-1 from
! CW Eq. 1.7 and Eq. 1.8
! Input:
!   a - Current Solution 
!   dx- Cell widths 
! Input/Output:
!   delj - Average slove 
!-----------------------------
use module_input
implicit none
real,intent(in),dimension(imin:imax)    :: a,dx
real,intent(inout),dimension(jmin-2:jmax+1) :: delmj 
real,dimension(jmin-2:jmax+1) :: one,two,three,diff1,diff2,delj
integer :: st0,stm1,stp1,en0,enm1,enp1,i
real :: sgn
! These value are the start (st#) and the end (en#) 
! of the parts of a and dx that we evaluate in Eq. 1.8
! We will assign value to the iMIN-1-iMAX+2 elements of 
! the delmj vector.
st0 = jmin-2 ! First cell I calculate a_j+1/2 in
en0 = jmax+1 ! Last  cell I calculate a_j+1/2 in
stm1= st0 -1
enm1= en0 -1
stp1= st0 +1
enp1= en0 + 1
one   = dx(st0:en0)/(dx(stm1:enm1)+dx(st0:en0)+dx(stp1:enp1))
two   = (2.*dx(stm1:enm1)+dx(st0:en0))/(dx(stp1:enp1)+dx(st0:en0))
three = (dx(st0:en0)+2.*dx(stp1:enp1))/(dx(stm1:enm1)+dx(st0:en0)) 
delj  = one*(two*(a(stp1:enp1)-a(st0:en0))+three*(a(st0:en0)-a(stm1:enm1)))
diff1 = a(stp1:enp1)-a(st0:en0)
diff2 = a(st0:en0)  -a(stm1:enm1)
do i=st0,en0
  if (diff1(i)*diff2(i).GT. 0.) then
    if (delj(i).LT.0.) then
      sgn = -1.
    else 
      sgn = 1.
    end if 
    delmj(i) = SIGN( MIN(ABS(delj(i)),2.*ABS(diff1(i)),2.*ABS(diff2(i))),sgn)
  else 
    delmj(i) = 0. 
  end if 
end do
end subroutine

subroutine MUSCL_interpol(a,dx,aL,aR,dela)
use module_input
implicit none
real,intent(in),dimension(imin:imax) :: a,dx
real,intent(inout),dimension(jmin-1:jmax+1) :: aL,aR,dela
integer :: i
do i=jmin-1,jmax+1
    if ((a(i+1)-a(i))*(a(i)-a(i-1)).GT.delta) then
        dela(i) = SIGN(MIN(2.*ABS(a(i)-a(i-1)),2.*ABS(a(i+1)-a(i)),.5*abs(a(i+1)-a(i-1))),a(i+1)-a(i-1))
    else 
        dela(i) = 0 
    end if 
end do 
aL = a(jmin-1:jmax+1) - .5*dela
aR = a(jmin-1:jmax+1) + .5*dela
end subroutine

subroutine OpenFiles()
implicit none
open(unit=1,file='Output/a.txt')
open(unit=2,file='Output/r.txt')
open(unit=3,file='Output/a_true.txt')
end subroutine

