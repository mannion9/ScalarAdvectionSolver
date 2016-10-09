module module_input
implicit none
!----------------------- USER INPUTS ----------------------------------------- ! 
integer,parameter :: N 	      = 100    ,& ! Number of cells + 1
					 itterMax = 1000   ,& ! Maximum itterations
					 writeStep= 10        ! Time steps difference between write outs
real,parameter    :: Courant  = .9     ,&
					 length   = 1.0    ,&
					 tmax     = 0.3       
character (LEN=*),parameter ::  Method = 'PPM'  ! ('Godunov','MUSCL','PPM')
!----------------------- Program Constants ----------------------------------------- !                                
real,parameter    :: delta    = 1.e-30      ! A small number
integer,parameter :: jmin     = 0        ,& ! Index of left  most LHS real  cell center
					 jmax     = N        ,& ! Index of right most RHS real  cell center
					 imin     = -3       ,& ! Index of left  most LHS ghost cell center
					 imax     = jmax+4      ! Index of right most RHS ghost cell center
      
contains   
subroutine InitCond(a,x)
implicit none
real,intent(inout),dimension(imin:imax) :: a
real,intent(in),dimension(imin:imax)    :: x
integer :: i
! Square Wave
do i=imin,imax
	if (x(i).LE.length*.25) then
		a(i) = delta
	else if (x(i).GT.length*.25 .AND. x(i).LE.length*.75) then
		a(i) = 1.
	else 
		a(i) = delta
	end if
end do
! Constant state
!a = 1.
! Sine Wave
!a  = (/(SIN(2*3.14*x(i)),i=imin,imax,1)/)
! Gaussian
a = (/(EXP(-((x(i)-.5*length)**(2)/(length/100))),i=imin,imax,1)/)
end subroutine
end module
