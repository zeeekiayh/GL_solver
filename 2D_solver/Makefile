PROGRAM = gl_fdm

PROGRAMc = $(PROGRAM).cpp
PROGRAMo = $(PROGRAM).o

SUBR = SC_classes_base.cpp SC_classes_derived.cpp he3bulk.cpp linear_eq_solver.cpp

SUBRo = $(SUBR:.cpp=.o)

LIBS = 

#FLAGS = -Wall -O3 -fopenmp

#mac osx needs pointing where to include eigen from; hence extra flag -I
FLAGS = -Wall -O3 #-fopenmp #-I/usr/local/include

#C = /usr/local/mpich-1.2.7/C/bin
#C = /usr/bin
compiler= g++  
#compiler = /opt/local/bin/g++-mp-9 # Mac Pro

.SUFFIXES : .cpp 
.cpp.o :
	$(compiler) $(FLAGS) -c $<

$(PROGRAM) : $(PROGRAMo) $(SUBRo)
	$(compiler) $(FLAGS) -o $(PROGRAM) $(PROGRAMo) $(SUBRo) $(LIBS)

#################################
#################################

clean_this: 
	rm     $(PROGRAMo) SC_classes_base.o SC_classes_derived.o 

clean:
	rm *.o;

plot5:
	python3 plot_solution.py output_OP5.txt
plot3:
	python3 plot_solution.py output_OP3.txt
