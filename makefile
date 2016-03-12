.cpp.o:	
	gcc -c $< -o $@ -lboost_iostreams -lboost_system -lboost_filesystem -lutil -fopenmp

# Abhaengigkeiten fuer das executable
fluidDemo: fluidDemo.o
	gcc $^ -o $@ -lboost_iostreams -lboost_system -lboost_filesystem -lutil -fopenmp


#Abhaengigkeiten fuer die object files
fluidDemo.o : fluidDemo.cpp fluidFuncs.cpp defines.h gnuplot-iostream.h


.PHONY: run
run: fluidDemo
	./fluidDemo


.PHONY: clean
clean:
	rm -f *.o fluidDemo *~
