CF 	= gfortran
objects = module_input.o main.o

execname: $(objects)
	$(CF) $(objects)
module_input.mod : module_input.o module_input.f90
	$(CF) -c module_input.f90
module_input.o: module_input.f90
	$(CF) -c module_input.f90
main.o: module_input.mod main.f90
	$(CF) -c main.f90
clean:
	rm $(objects) module_input.mod 
	rm Output/*.txt

