# Makefile 

# Precision
PFLAG = -Ddouble_precision
#PFLAG = -Dsingle_precision

# Compiler
CPP=g++

# Optimization flags
OFLAG=-O3

# Debug flags
DFLAG=-g3 -pg -Ddebug_mode

# Static flag
SFLAG=-static

# Warning flags
WFLAG=-Wall -Wno-deprecated

# Thread flags
TFLAG=-lpthread -fopenmp

# Optimization
FLAGS=$(OFLAG) $(SFLAG) $(TFLAG) $(WFLAG) $(PFLAG)

# Debugging
#FLAGS=$(DFLAG) $(SFLAG) $(TFLAG) $(WFLAG) $(PFLAG)


OBJECTS=\
Vector.o \
Matrix.o \
Cell.o \
Timer.o \
Parameters.o \
RegularMesh.o \
AdaptiveMesh.o \
main.o

HEADERS=\
Global.h \
Vector.h \
Cell.h \
Timer.h \
Parameters.h \
Matrix.h \
RegularMesh.h \
AdaptiveMesh.h

# Compile source files
.cpp.o:
	$(CPP) -c $(FLAGS) $<

# Link objects
lara: $(OBJECTS) $(HEADERS)
	$(CPP) -o lara $(FLAGS) $(OBJECTS) -I /usr/include

clean:	
	rm -f *.o
	rm -f *~
	rm -f lara


