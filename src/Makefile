## need module load mkl gcc48

#INCL = -I/usr/include/gcc/darwin/3.3-fast \
	 -I/System/Library/Frameworks/vecLib.framework/Headers/

INCL = -I$(MKL_HOME)/include

LIB = -L$(MKL_HOME)/lib/intel64

CFLAGS =    -O3  -funroll-loops -fstrict-aliasing -mtune=native -m64 -ffast-math #-Wno-long-double #-DNDEBUG 

LIBS =  -lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lpthread -lm  -lstdc++

CC       = g++

MPE_C =       mpe.cpp protein.cpp BD.cpp pdb.cpp mcoeff.cpp gradcoeff.cpp \
	      xform.cpp shcoeff.cpp tricoeff.cpp torqcoeff.cpp rotcoeff.cpp \
		transcoeff.cpp

MPE_OBJS      =       $(MPE_C:.cpp=.o)

mpe:  $(MPE_OBJS) Makefile
	rm -f ../bin/mpe
	$(CC)  $(CFLAGS) $(MPE_OBJS) $(LIB) $(LIBS) -o mpe
	mv mpe ../bin/

.cpp.o:
	$(CC) $(CFLAGS) $(INCL) -c $<

clean:
	rm *.o

depend:
	makedepend $(INCL) $(MPE_C) $(FIT_C) $(PARAMFIT_C) $(TESTFIT_C)

# DO NOT DELETE THIS LINE -- make depend depends on it.

mpe.o: mpe.h util.h mcoeff.h gradcoeff.h tricoeff.h
mpe.o: torqcoeff.h xform.h shcoeff.h rotcoeff.h transcoeff.h protein.h pdb.h   
protein.o: protein.h pdb.h util.h 
protein.o: torqcoeff.h xform.h shcoeff.h rotcoeff.h transcoeff.h  mpe.h mcoeff.h gradcoeff.h tricoeff.h
BD.o: BD.h util.h 
BD.o: pdb.h protein.h mpe.h mcoeff.h gradcoeff.h
BD.o: tricoeff.h torqcoeff.h xform.h shcoeff.h rotcoeff.h transcoeff.h
pdb.o: pdb.h util.h 
pdb.o: protein.h mpe.h mcoeff.h gradcoeff.h
pdb.o: tricoeff.h torqcoeff.h xform.h shcoeff.h rotcoeff.h transcoeff.h
mcoeff.o: mcoeff.h util.h 
mcoeff.o: shcoeff.h
gradcoeff.o: gradcoeff.h tricoeff.h mcoeff.h util.h 
xform.o: xform.h mcoeff.h util.h 
xform.o: gradcoeff.h tricoeff.h torqcoeff.h
xform.o: shcoeff.h rotcoeff.h transcoeff.h
shcoeff.o: shcoeff.h util.h 
tricoeff.o: tricoeff.h mcoeff.h util.h 
torqcoeff.o: torqcoeff.h tricoeff.h mcoeff.h util.h 
torqcoeff.o: shcoeff.h
rotcoeff.o: rotcoeff.h mcoeff.h util.h 
rotcoeff.o: tricoeff.h shcoeff.h
transcoeff.o: shcoeff.h util.h 
transcoeff.o: mcoeff.h transcoeff.h
mcoeff.o: mcoeff.h util.h 
mcoeff.o: shcoeff.h
gradcoeff.o: gradcoeff.h tricoeff.h mcoeff.h util.h 
shcoeff.o: shcoeff.h util.h 
shcoeff.o: mcoeff.h
tricoeff.o: tricoeff.h mcoeff.h util.h 
transcoeff.o: shcoeff.h util.h 
transcoeff.o: mcoeff.h transcoeff.h
mcoeff.o: mcoeff.h util.h 
mcoeff.o: shcoeff.h
gradcoeff.o: gradcoeff.h tricoeff.h mcoeff.h util.h 
shcoeff.o: shcoeff.h util.h 
shcoeff.o: mcoeff.h
tricoeff.o: tricoeff.h mcoeff.h util.h 
