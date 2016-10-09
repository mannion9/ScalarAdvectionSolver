This code will solve the linear advection equation given by finite volume considerations.
The linear advection equation to be solved it.

	a_t + c * a_x = 0

where a subscript indicates a partial derivative, c is the speed of propogation, and a is
the scalar function a = a(x,t).

The files contained in this package are

src		 - Contains source code
Plot.py 	 - Plotting tool
makefile 	 - Makes files

In the src/ you will find four files

main.f90         - Main program driver
module_input.f90 - Contains all user inputs and paramters
Output		 - Directory in which data is written to
makefile	 - Makes source code

The dependencies are

gfortran (or other fortran compiler, specificed in /src/makefile)
python 3.0
matplotlib

The code is build so that the user should only interact with the module_input.f90 file.
This file contains two sections, the initialization and the initial condition.

The initialization contains two sections. The first is the user inputs and the second is
program constants, which should never be changed. The user inputs are

N 	 - Number of grid points
itterMax - Maximum number of time steps
writeStep- Number of time steps between write out dump
Courant  - CLF condition < 1.0
length   - Length of domain (begining at zero)
tmax     - Time cut off
Solver   - Solution method

The program constants are 

jmin     - Array index of left  most real  cell center
jmax     - Array index of right most real  cell center
imin     - Array index of left  most ghost cell center
imax  	 - Array index of right most ghost cell center

The user inserts the chosen initial conditiion function the subroutine InitCond. This
subroutine has two inputs, a the solution, and r the cell centers. These should not 
be changed. The user has freedom to define the function as is desired. A series of 
example initial conditions are contained with the subroutine.
 

 
