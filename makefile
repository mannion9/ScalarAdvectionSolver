# Goes into my /src directory, and calls the make file
all:
	cd src; make 
	cd src; ./a.out
	cp -r src/Output/ Output/
clean: 
	rm -rf Output
	cd src; make clean

